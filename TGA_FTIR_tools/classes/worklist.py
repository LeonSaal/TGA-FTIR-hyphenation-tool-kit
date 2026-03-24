import logging
from pathlib import Path
import pickle
import re
from dataclasses import dataclass, field
from typing import List, Literal, Optional, Union
import matplotlib.pyplot as plt
import pandas as pd
from itertools import chain, zip_longest
from ..classes import Sample, Baseline
from ..config import PATHS, DEFAULTS, MERGE_CELLS
from ..fitting import get_presets, robustness
from ..input_output import samplelog, time
from ..plotting import plot_results, plot_robustness, plots
from concurrent.futures import ProcessPoolExecutor
from inspect import signature
from typing import get_args
from types import NoneType
import numpy as np

logger = logging.getLogger(__name__)
import os

# for concurrent.futures
# def fit_sample(args):
#     sample, reference, presets, mod_samples, kwargs = args
#     return sample.fit(
#         reference,
#         presets=presets,
#         mod_sample=mod_samples,
#         save=False,
#         **kwargs
#     )


@dataclass
class Worklist:
    names: Union[List[Union[Sample, str]], str, Sample] = field(default_factory=list)
    name: Optional[str] = "worklist"
    profile: str = DEFAULTS["profile"]
    aliases: list = field(default_factory=list)
    samples: list = field(default_factory=list)
    _results: dict = field(default_factory=lambda: {"fit": pd.DataFrame(), "robustness": pd.DataFrame()})

    def __post_init__(self):
        len_names = 1 if isinstance(self.names, str) else len(self.names)
        match self.aliases:
            case str():
                self.aliases = [self.aliases] * len_names 

        match self.samples:
            case str():
                self.samples = [self.samples] * len_names 

        match self.names:
            case str():
                self.names = [Sample(self.names, profile=self.profile, alias = self.aliases[0], sample = self.samples[0])]
            case list():
                self.names = [
                    (
                        name
                        if isinstance(name, Sample)
                        else Sample(name, profile=self.profile, alias = alias if alias else name, sample=sample)
                    )
                    for name, alias, sample in zip_longest(self.names, self.aliases, self.samples)
                ]
            case Sample():
                self.names = [self.names]
            case _:
                logger.warning(
                    f"{self.names!r} has wrong input type ({type(self.names)}). Supply one of: str, Sample, list[str|Sample]."
                )

    def __add__(self, other):
        if isinstance(other, Sample):
            return Worklist(self.names + [other], name=f"{self.name}+{other.name}")
        elif isinstance(other, Worklist):
            return Worklist(
                self.names + other.names, name=f"{self.name}+{other.name}"
            )
        else:
            logger.warning("Can only add Sample or other Worklist")

    def __repr__(self) -> str:
        return (
            f"Worklist(name='{self.name}', samples= \n"
            + "\n".join(
                [f"[{i}] {sample.__repr__()}" for i, sample in enumerate(self.names)]
            )
            + "\n)"
        )

    def __getitem__(self, i:Union[int, slice, str, list]):
        if type(i) == int:
            return self.names[i]
        elif type(i) == slice:
            return Worklist(
                names=self.names[i], name=f"{self.name}_({i.start}-{i.stop})"
            )
        elif type(i) == str:
            if i in (d := {sample.name: sample for sample in self.names}):
                return d[i]
        elif type(i) == list and i is not []:
            samples = []
            for elem in i:
                samples.append(self.__getitem__(elem))
            if len(samples) == 1:
                return samples[0]
            return Worklist(names=samples, name=f"{self.name} {repr(i)}")

    def __iter__(self):
        yield from self.names

    def get(self, pattern: str, attr: str = "name"):
        return Worklist(
            names=[
                sample
                for sample in self.names
                if re.search(pattern, sample.__dict__[attr])
            ],
            name=pattern,
        )

    def append(self, other) -> NoneType:
        if isinstance(other, Worklist):
            self.names.extend(other.names)
        if isinstance(other, Sample):
            self.names.append(other)

    def fit(
        self, reference: str, save:bool=True, mod_samples:bool=True, presets:dict[pd.DataFrame]=None, **kwargs
    ) -> pd.DataFrame:

        # load default presets
        if not presets:
            presets = get_presets(reference)

        if presets is None:
            logger.error(
                "Presets neither were supplied, nor could be loaded from defaults."
            )
            return

        # make subdirectory to save data
        if save:
            path = PATHS["fitting"] / f"{time()}_{reference}_{self.name}"
            os.makedirs(path)
            os.chdir(path)
        else:
            path = None

        # cycling through samples
        fits = []
        for sample in self.names:
            # writing data to output DataFrames
            if reference not in sample.results["fit"] or not mod_samples:
                fits.append(sample.fit(
                    reference,
                    presets=presets,
                    mod_sample=mod_samples,
                    **kwargs,
                    save=save,
                    save_dir = path
                ))      
        res = self.results["fit"] if mod_samples else pd.concat(fits).pint.convert_object_dtype()

        if save:
            with pd.ExcelWriter("summary.xlsx") as writer:
                res.to_excel(writer)

        os.chdir(PATHS["home"])
        return res
    @property
    def results(self):
        out = {"fit": None, "robustness": None}
        oldres = list(
            chain(
                *[
                    [
                        res
                        for res in s.results["fit"].values()
                        if isinstance(res, pd.DataFrame)
                    ]
                    for s in self
                ]
            )
        )

        if oldres:
            out["fit"] = pd.concat(oldres)

        # conform robustness df to long format
        if self._results["robustness"].empty:
            out["robustness"] = {}
        else:
            idx_cols = self._results["robustness"].index.names
            new_col = "varied"
            out["robustness"] = (
                self._results["robustness"]
                .astype(np.float64)
                .reset_index()
                .melt(
                    id_vars=idx_cols,
                    ignore_index=False,
                    var_name=new_col,
                    value_name="mmol_per_mg",
                )
            ).set_index(idx_cols+[new_col])
        return out

    def robustness(self, reference: str, plot:bool=True, save:bool=True, var_T:float=10, var_rel:float=0.3, **kwargs) -> pd.DataFrame:
        self._results["robustness"], summary = robustness(self, reference, **kwargs)

        if save:
            path = PATHS["robustness"] / "_".join([time(), reference, f"{var_T=}_{var_rel=}"])
            os.makedirs(path)
            os.chdir(path)
            logger.info(f"Plots and results are saved.'{path=}'.")
            with pd.ExcelWriter("robustness_in_mmol_per_mg.xlsx") as writer:
                summary.to_excel(writer, sheet_name="summary", merge_cells=MERGE_CELLS)
                self._results["robustness"].to_excel(writer, sheet_name="data", merge_cells=MERGE_CELLS)

        if plot:
            self.plot("robustness", save=True, save_dir=Path("."))
        os.chdir(PATHS["home"])
        
        return self._results["robustness"], summary

    def plot(self, plot:Literal["fit", "robustness", "TG", "EGA", "DTG", "heat_flow"]=None, ax:Union[None, plt.Axes]=None, save:bool=False, save_dir:Union[str, None]=None, reference_name:str=None, **kwargs) -> NoneType:
       
        if plot in ["fit", "robustness"]:
            results = self.results[plot]
            if results is None:
                logger.warning(f"No results to plot. Run .{plot}().")
                return
            results = (
                results.pint.convert_object_dtype().astype(np.float64).reset_index()
            )
            avail_gases = results.gas.unique()
            avail_references = results.reference.unique()

            # validate inputs
            references = kwargs.get("references", avail_references)
            references = [ref for ref in references if ref in avail_references]
            gases = kwargs.get("gases", avail_gases)
            gases = [g for g in gases if g in avail_gases]
            if not references or not gases:
                logger.error(
                    f"There are no valid references or gases to plot.\nValid references are {avail_references!r}\nValid gases are {avail_gases!r}"
                )
                return
            # subset dataset
            data = results.query("gas in @gases & reference in @references")

            # validate plot arguments
            sig = signature(plot_results).parameters
            args = {
                key: (
                    kwargs.get(key)
                    if kwargs.get(key) in get_args(sig[key].annotation)
                    else sig[key].default
                )
                for key in sig.keys()
                if key != "data"
            }
            fig, ax = plot_results(data, **args)
        else:
            if not ax:
                fig, ax = plt.subplots()
            plots(self.names, plot, ax, **kwargs)

        if save:
            path_plot = save_dir if save_dir else PATHS["plots"] / plot
            if not path_plot.exists():
                path_plot.mkdir(parents=True)

            path_pic = path_plot / "_".join([time(), self.name])
            fig.savefig(path_pic)
        return ax

    def __len__(self) -> int:
        return len(self.names)

    def corr(self, baselines, corrs:dict, **kwargs) -> NoneType:
        if type(baselines) == list:
            if (ls := len(self)) != (lb := len(baselines)):
                logger.error(
                    f"Lengths of samples ({ls}) and baselines ({lb}) do not match."
                )
                return
        else:
            baselines = list(baselines)*len(self)
        
        for sample, baseline in zip(self.names, baselines):
            match baseline:
                case str():
                    baseline = Baseline(baseline)
                case Sample():
                    pass
                case Baseline():
                    pass
                case _:
                    logger.error(f"{baseline!r} is of wrong type! Allowed are str (sample name), Sample and Baseline. Skipping.")
                    continue
            sample.corr(baseline, corrs, **kwargs)

    def pop(self, i: int) -> NoneType:
        return self.names.pop(i)

    def save(
        self,
        how: Literal["samplelog", "pickle"] = "samplelog",
        **kwargs,
    ) -> NoneType:
        if how == "samplelog":
            samplelog(self.info, how="samplelog")

        elif how == "pickle":
            self.to_pickle()

    def from_pickle(name: str) -> NoneType:
        if (p:= PATHS["output"] / f"{name}.pkl").exists():
            with open(p , "rb") as inp:
                obj = pickle.load(inp)
            if isinstance(obj, Sample):
                return obj
        else:
            logger.error(f"{p.as_posix()!r} does not exist!")

    def to_pickle(self):
        if not PATHS["output"].exists():
            os.makedirs(PATHS["output"])
        with open(PATHS["output"] / f"{self.name}.pkl", "wb") as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
    
    @property
    def info(self) -> pd.DataFrame:
        return pd.concat([sample.info.to_row() for sample in self.names])
    
    @property
    def profiles(self) -> list:
        return [s.profile for s in self.names]
        
    def calibrate(self, **kwargs):
        from ..calibration import calibrate
        return calibrate(worklist = self, mode="recalibrate", profile=self.profile, **kwargs)
        
    def from_samplelog(sheet_name:Union[str, int, list,None]=0):
        worklist = samplelog(sheet_name=sheet_name)
        samples = worklist.index.to_list()
        profiles = worklist.profile
        aliases = worklist.alias

        return Worklist(samples, name = sheet_name, profiles=profiles, aliases=aliases)
