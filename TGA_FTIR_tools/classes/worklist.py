import logging
import pickle
import re
from dataclasses import dataclass, field
from typing import List, Literal, Optional, Union
import matplotlib.pyplot as plt
import pandas as pd
from itertools import chain
from ..classes import Sample
from ..config import PATHS, DEFAULTS
from ..fitting import get_presets, robustness
from ..input_output import samplelog, time
from ..plotting import bar_plot_results, plot_robustness, plots
from concurrent.futures import ProcessPoolExecutor

logger = logging.getLogger(__name__)
import os

# for concurrent.futures
def fit_sample(args):
    sample, reference, presets, mod_samples, kwargs = args
    return sample.fit(
        reference,
        presets=presets,
        mod_sample=mod_samples,
        save=False,
        **kwargs
    )

@dataclass
class Worklist:
    samples: Union[List[Union[Sample,str]],str,Sample] = field(default_factory=list)
    name: Optional[str] = "worklist"
    profile: str = DEFAULTS["profile"]
    _results: dict = field(default_factory=lambda: {"fit": {}, "robustness": {}})

    def __post_init__(self):
        match self.samples:
            case str():
                self.samples = [Sample(self.samples, profile=self.profile)]
            case list():
                self.samples = [sample if isinstance(sample, Sample) else Sample(sample, profile=self.profile) for sample in self.samples]
            case Sample():
                self.samples = [self.samples]
            case _:
                logger.warning(f"{self.samples!r} has wrong input type ({type(self.samples)}). Supply one of: str, Sample, list[str|Sample].")
        if any([sample.profile!=self.profile for sample in self.samples]):
            logger.warning(f"Some of the supplied Samples have another profile than {self.profile!r}")

    def __add__(self, other):
        if isinstance(other, Sample):
            return Worklist(self.samples + [other], name=f"{self.name}+{other.name}")
        elif isinstance(other, Worklist):
            return Worklist(
                self.samples + other.samples, name=f"{self.name}+{other.name}"
            )
        else:
            logger.warning("Can only add Sample or other Worklist")

    def __repr__(self) -> str:
        return (
            f"Worklist(name='{self.name}', samples= \n"
            + "\n".join(
                [f"[{i}] {sample.__repr__()}" for i, sample in enumerate(self.samples)]
            )
            + "\n)"
        )

    def __getitem__(self, i):
        if type(i) == int:
            return self.samples[i]
        elif type(i) == slice:
            return Worklist(
                samples=self.samples[i], name=f"{self.name}_({i.start}-{i.stop})"
            )
        elif type(i) == str:
            if i in (d := {sample.name: sample for sample in self.samples}):
                return d[i]
        elif type(i) == list and i is not []:
            samples = []
            for elem in i:
                samples.append(self.__getitem__(elem))
            if len(samples) == 1:
                return samples[0]
            return Worklist(samples=samples, name=f"{self.name} {repr(i)}")

    def __iter__(self):
        yield from self.samples

    def get(self, key: str, attr: str = "name"):
        return Worklist(
            samples=[
                sample
                for sample in self.samples
                if re.search(key, sample.__dict__[attr])
            ],
            name=key,
        )

    def append(self, other) -> None:
        if isinstance(other, Worklist):
            self.samples.extend(other.samples)
        if isinstance(other, Sample):
            self.samples.append(other)

    def fit(
        self, reference: str, save=True, mod_samples=True, presets=None, **kwargs
    ) -> pd.DataFrame:

        # load default presets
        if not presets:
            presets = get_presets(reference)

        if presets is None:
            logger.error("Presets neither were supplied, nor could be loaded from defaults.")
            return

        # make subdirectory to save data
        if save:
            path = PATHS["fitting"] / f"{time()}{reference}_{self.name}"
            os.makedirs(path)
            os.chdir(path)

        # cycling through samples
        samples_to_fit = [s for s in self if reference not in s.results["fit"]]
        args_list = [(s, reference, presets, mod_samples, kwargs) for s in samples_to_fit]

        with ProcessPoolExecutor() as e:
            futures = [e.submit(fit_sample, args) for args in args_list]

        os.chdir(PATHS["home"])
        return [f.result() for f in futures]  # self.results["fit"]

    @property
    def results(self):
        out = {"fit": None, "robustness": None}
        oldres = list(chain(*[[res for res in s.results["fit"].values() if isinstance(res, pd.DataFrame)] for s in self]))

        if oldres:
            out["fit"] = pd.concat(oldres)
        out["robustness"] = self._results["robustness"]

        return out

    def robustness(self, reference: str, plot=True, **kwargs) -> pd.DataFrame:
        self._results["robustness"] = robustness(self, reference, **kwargs)
        if plot:
            self.plot("robustness")
        return self._results["robustness"]

    def plot(self, plot=None, ax=None, save=False,**kwargs) -> None:
        if not ax:
            fig, ax = plt.subplots()

        if plot == "robustness" and "robustness" in self._results:
            logger.warning("Under Maintenance")
            # plot_robustness(self._results["robustness"][0])
        elif plot == "results":
            logger.warning("Under Maintenance")
            # bar_plot_results(self, **kwargs)
        else:
            plots(self.samples, plot,ax, **kwargs)

        if save:
            path_plot = PATHS["plots"] / plot
            if not path_plot.exists() :
                path_plot.mkdir(parents=True)

            path_pic = path_plot / "_".join([time(), self.name])
            fig.savefig(path_pic)



    def __len__(self) -> int:
        return len(self.samples)

    def corr(self, references=None, plot=False, **kwargs) -> None:
        if type(references) == list:
            if (ls := len(self)) != (lr := len(references)):
                logger.error(
                    f"Lengths of samples ({ls}) and references ({lr}) do not match."
                )
                return
        elif type(references) == str:
            references = [references for _ in self.samples]

        elif not references:
            references = [None for _ in self.samples]

        for sample, baseline in zip(self.samples, references):
            sample.corr(baseline, plot=plot, **kwargs)

    def pop(self, i: int) -> None:
        self.samples.pop(i)

    def save(
        self,
        fname: str = None,
        how: Literal["samplelog", "pickle"] = "samplelog",
        **kwargs,
    ) -> None:
        if how == "samplelog":
            # collect sample info
            infos = {sample.name: sample.info for sample in self.samples}
            info = {name: {key: str(value) for key, value in info.__dict__.items() if value is not None} for name, info in infos.items()}
            
            #make df and save
            data = pd.DataFrame.from_dict(info, orient="index")
            samplelog(data, how = "samplelog")

        elif how == "pickle":
            path_output = PATHS["output"]
            if not path_output.exists():
                os.makedirs(path_output)
            if not fname:
                fname = self.name
            path = path_output / f"{fname}.wkl"
            with open(path, "wb") as output:
                pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

            # for sample in self.samples:
        #     sample.save(**kwargs)

    def load(self, fname: str) -> None:
        with open(PATHS["output"] / f"{fname}.wkl", "rb") as inp:
            obj = pickle.load(inp)
        for key in obj.__dict__:
            self.__dict__[key] = obj.__dict__[key]
