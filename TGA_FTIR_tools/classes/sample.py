import copy
import json
import logging
import os
import pickle
import re
from dataclasses import InitVar, dataclass, field
from types import NoneType
from typing import Any, Literal, Mapping, Optional, List
from webbrowser import get
import matplotlib.pyplot as plt
from typing import get_args
from inspect import signature
import inspect

import numpy as np
import pandas as pd
from traitlets import signature_has_traits
import pint_pandas

import pint
from pint import DimensionalityError
ureg = pint.get_application_registry()
import scipy as sp
from scipy.signal import savgol_filter

from ..config import DEFAULTS, PATHS, SAVGOL, update_config
from ..input_output import (FTIR, TGA, corrections, general, mass_step,
                            read_data, samplelog, read_profile_json)
from ..plotting import plot_dweight, plot_mass_steps, plot_calibration_single, plot_residuals_single, plot_calibration_combined, plot_results
from ..utils import select_import_profile, check_profile_exists
from .info import SampleInfo

WINDOW_LENGTH_REL = SAVGOL.getfloat("window_length_rel")
POLYORDER = int(SAVGOL.getfloat("polyorder"))


logger = logging.getLogger(__name__)


@dataclass
class Sample:
    name: str = field(default=None)
    ega: Optional[pd.DataFrame] = field(default=None)
    tga: Optional[pd.DataFrame] = field(default=None)
    linreg: Optional[pd.DataFrame] = field(default=None)
    alias: str = field(default=None)
    reference: str = field(default=None)
    sample: str = field(default=None) 
    run: str = field(default=None)
    baseline = None
    _info:dict =field(default_factory=dict) 
    results: dict = field(default_factory=lambda: {"fit": {}, "robustness": {}})
    profile: str = DEFAULTS["profile"]

    def __post_init__(self, **kwargs):
        if self.name:
            self.init(**kwargs)

    def init(self, **kwargs):
        # check and get import profile
        self.check_profile()
        self.profile_data = read_profile_json(self.profile)
        logger.info(f"Initializing {self.name!r} using profile {self.profile!r}.")

        # load data
        try:
            self.__dict__.update(read_data(self.name, profile=self.profile))
        except Exception as e:
            logger.error(f"Could not read any data for {self.name!r}. Have you specified the right {self.profile=}")
            return
        self._info = SampleInfo(self.name, info=self._info, profile=self.profile)

        self.post_init_tga()
        self.post_init_ega(**kwargs)

        if self.ega is None and self.tga is None:
            logger.error(f"Failed to initialize '{self.name}'.")
            return
        else:
            logger.info(f"'{self.name}' successfully initialiazed.")

        # assigning alias, sample, run from samplelog
        log = samplelog(create=False)
        self.alias = log.get("alias", {}).get(self.name, self.name) if not self.alias else str(self.alias)
        self.sample = log.get("sample", {}).get(self.name, self.alias) if not self.alias else str(self.alias)
        self.run = log.get("run", {}).get(self.name, 0) if not self.run else self.run
        self.raw = copy.deepcopy(self)

        if corrs:=self.info.get("corrections"):
            self.baseline = Baseline().from_sample(self, corrs=corrs)
            self.__dict__.update(corrections.corr_baseline(self.raw, self.baseline))
            self.reference = self.baseline.name

    def post_init_tga(self):
        if self.tga is None:
            logger.error(f"No TGA data found for {self.name!r}.")
            self._info = SampleInfo(name=self.name, alias=self.alias, profile=self.profile)
            return
        else:
            logger.info("TGA data found.")

        logger.debug("Checking required columns in TGA data.")
        
        # deriving TG info
        if not self._info["initial_mass"] and "sample_mass" in self.tga.columns:
            self._info["initial_mass"]=self.tga["sample_mass"].iloc[0]

        # calculate sample mass with mass loss or percentage mass and initial mass
        try:
            if "sample_mass" not in self.tga.columns and "mass_loss" in self.tga.columns:
                if "initial_mass" in self.info:
                    self.tga["sample_mass"] =  self.tga["mass_loss"] + self.info["initial_mass"]

            elif "sample_mass" in self.tga.columns and  self.tga.sample_mass.dtype == "pint[%]":
                if "initial_mass" in self.info:
                    self.tga["sample_mass"] = self.info["initial_mass"] * (self.tga["sample_mass"].astype("pint[dimensionless]"))
        except DimensionalityError as e:
            logger.error(f"Failed. {e}")
        self._info.steps_idx.update({})

        # check required columns
        if missing:=self.missing_tga_columns():
            logger.error(f"Required column(s) {missing} not found in TGA data of {self.name!r}")
        else:   
            try:
                window_length = int(self.tga.index.size * WINDOW_LENGTH_REL/2)*2 + 1
                dtg = savgol_filter(
                    self.tga["sample_mass"].values._data, window_length, POLYORDER, deriv=1
                )
                dtype =  self.tga["sample_mass"].dtype
                self.tga = self.tga.assign(dtg = -pd.Series(dtg).astype(dtype))
            except ValueError as e:
                print(e)
            except AttributeError as e:
                print(e)

    def post_init_ega(self, **kwargs):
        if self.ega is None:
            logger.error(f"No EGA data found for {self.name!r}.")
            return
        
        self._info["gases"] = self.ega.columns[1:].to_list()
        
        # load calibration
        self.calibrate(mode="load", profile=self.profile)
        
        gases = self._info["gases"]
        if self.linreg is not None:
            calibrated_gases = self.linreg.index 
            gases_annotated = [gas+ ("*"if  gas in calibrated_gases else "")
                        for gas in gases
                    ]
            cali_msg = " (* and calibrated) "
        else:
            gases_annotated = gases
            cali_msg = ""
        logger.info(
            "EGA data found{} for {}.".format(cali_msg,", ".join(gases_annotated))
        )

        if "sample_temp" not in self.ega and "time" in self.ega:
            self.ega.time = self.ega.time.astype(self.tga.time.dtype)
            self.ega = self.ega.set_index("time").reindex(self.tga.time, method="nearest").reset_index()
            try:
                self.ega = pd.merge(
                    self.tga.filter(
                        ["time", "sample_temp", "reference_temp", "sample_mass"], axis=1
                    ),
                    self.ega,
                    how="left",
                    on="time",
                )#.dropna(axis=0)
            except Exception as e:
                logger.debug(f"Unable to merge temperature data. {e}")
        try:
            self.dry_weight(plot=False, **kwargs)
        except Exception as e:
            logger.warning(f"Failed to derive dry weight: {e}")
                    
    # initialize object from pickle file
    def from_pickle(name):
        if (p:= PATHS["output"] / f"{name}.pkl").exists():
            with open(p , "rb") as inp:
                obj = pickle.load(inp)
            if isinstance(obj, Sample):
                return obj
        else:
            logger.error(f"{p.as_posix()!r} does not exist!")

    def __repr__(self):
        attr_names = ["name", "alias", "profile","sample", "reference", "run"]
        attr_vals = ", ".join(
            [f"{key}={repr(self.__dict__[key])}" for key in attr_names]
        )
        return f"Sample({attr_vals})"

    def __add__(self, other):
        from .worklist import Worklist
        if isinstance(other, Sample):
            return Worklist([self, other], name=f"{self.name}+{other.name}")
        elif isinstance(other, Worklist):
            return other + self
        else:
            logger.debug("Can only add to Sample- or Worklist-type")

    def get(self, name):
        return self.__dict__.get(name)

    @property
    def info(self):
        self._info.alias = self.alias
        # deriving IR info
        try:
            self._info.update(FTIR.FTIR_info(self))
        except Exception as e:
            logger.debug(f"Unable to derive EGA-info. {e}")
            #raise e

        return self._info

    @property
    def reference_mass(self):
        ref_mass_name=self.info.reference_mass_name
        step_data = self.step_data()
        return step_data[step_data.step == ref_mass_name].sample_mass.values[0]

    def step_data(self, ref_mass_name=None):
        idxs = self.info.steps_idx
        step_data = self.tga.iloc[list(idxs.values()),:].assign(step= idxs.keys()).sort_index()
        step_data["mass_loss"] = step_data.sample_mass.diff()
        
        # calc relative mass loss 
        if not ref_mass_name:
            ref_mass_name=self.info.reference_mass_name
        if ref_mass_name in step_data.step:
            ref_mass = step_data[step_data.step==ref_mass_name].sample_mass
            step_data["mass_loss_rel"] = step_data.mass_loss / ref_mass

        return step_data

    def check_profile(self):
        # check if profile is supplied
        if not check_profile_exists(self.profile):
            self.profile = select_import_profile()
            update_config(section="defaults",key="profile", value=self.profile)
            logger.info(f"Updated profile in {PATHS["ini"].name!r} to {self.profile!r}.")
        

    def missing_tga_columns(self):
        "check if required columns are present in TGA data"
        required_columns = [
            "sample_mass",
            "sample_temp",
            "time"
        ]
        missing = [col for col in required_columns if col not in self.tga.columns]
        return missing if missing else None

    def corr(self, baseline, corrs: dict, **kwargs):
        """_summary_

        Args:
            baseline (Baseline | Sample): Initialized Baseline used for correction.
            corrs (dict): Dict of structure {"ega": {signal: fun}, "tga": ...}
            fun is function with signature fun(x, y, **kwargs) -> z that takes in numeric vectors of same length and returns one vector of same length. 
            x and y are the respective signal column for the sample to be corrected and the associated baseline data respectively.
        """        
        if not self.baseline:
            self.baseline = Baseline(name="from baseline")

        for data, traces in corrs.items():
            to_corr = self.__dict__.get(data)
            use_to_corr = baseline.__dict__.get(data)
            if not isinstance(to_corr, pd.DataFrame) and not isinstance(use_to_corr, pd.DataFrame):
                logger.warning(f"{data!r} not in sample or baseline or not of the correct type (required: pandas.DataFrame). Skipping.")
                continue
            
            to_update = pd.DataFrame()
            for signal, fun in traces.items():
                if signal not in to_corr or signal not in use_to_corr:
                    logger.warning(f"{signal!r} not in sample or baseline. Skipping.")
                    continue
                to_update[signal] = fun(to_corr[signal], use_to_corr[signal], **kwargs)

            self.__dict__[data].update(to_update)
            index_cols = self.__dict__[data].filter(items=['time', 'sample_temp', 'reference_temp', 'sample_mass'])
            self.baseline.__dict__[data] = pd.concat([index_cols, to_update], axis=1)
            self.baseline.info["gases"] = to_update.columns.to_list()
            

    def get_value(self, *values, which="sample_mass", at="sample_temp"):
        "extract values from TG data at e.g. certain temperatures"
        new_idx = pd.Series(values, dtype=np.float64, name="sample_temp")
        tmp = self.tga[[at, which]].copy().sort_values(at)
        tmp[at] = tmp[at].astype(np.float64)
        tmp = pd.merge_asof(new_idx, tmp, on=at)
        tmp[at] = tmp[at].astype(self.tga[at].dtype)
        return tmp
    

    def dry_weight(self, plot=False, **kwargs):
        "determine dry point and mass of sample"

        try:
            TGA.dry_weight(self, **kwargs)
            logger.info(
                "'.info' was updated. To store these in Samplelog.xlsx run '.save()'"
            )
            if plot:
                plot_dweight(self, **kwargs)

        except Exception as e:
            logger.error(f"Failed to derive dry weight. {e}")

    def mass_step(
        self, ax=None, plot=True, min_rel_height=0.2, width_T=np.array([0, np.inf]),**kwargs
    ):
        #calculate width in terms of indices from supplied temperature width
        step_height, rel_step_height, step_starts_idx, step_ends_idx, peaks_idx = mass_step(
            self,
            min_rel_height=min_rel_height,
            width_T=width_T,
            rel_height=.7,
            **kwargs,
        )
        steps = {f"Step {i}": idx for i, idx in enumerate(step_ends_idx, start=1)}
        self._info.steps_idx.update(steps)
        if plot:
            self.plot(
                "mass_steps",
                ax=ax,
                steps=self.step_data().sample_temp.astype(np.float64),
                y_axis="orig",
            )
        return step_height, rel_step_height, step_starts_idx, step_ends_idx, peaks_idx

    def plot(self, plot=Literal["TG","mass_steps","heat_flow","EGA", "cumsum","EGA_to_DTG","fit", "calibration", "results"], ax=None, save=False,reference=None, **kwargs):
        from ..plotting import FTIR_to_DTG, plot_fit, plot_FTIR, plot_TGA
        "plotting TG and or IR data"
        options = [
            "TG",
            "mass_steps",
            "heat_flow",
            "EGA",
            
            "cumsum",
            "EGA_to_DTG",
            "fit",
            "calibration",
            "results"
        ]
        if plot not in options:
            logger.warning(f"{plot} not in supported {options=}.")
            return
        needs_tg = ["TG", "heat_flow", "EGA_to_DTG", "mass_steps"]
        needs_ega = ["EGA",  "cumsum", "EGA_to_DTG"]
        needs_cali = [ "EGA_to_DTG","calibration"]
        
        # check requisites
        if (self.tga is None) and (plot in needs_tg):
            logger.warning("Option unavailable without TGA data.")
            return
        if self.ega is None and (plot in needs_ega):
            logger.warning("Option unavailable without EGA data.")
            return
        if self.linreg is None and (plot in needs_cali):
            logger.warning("Option unavailable without calibration.")
            return

        if not isinstance(ax, plt.Axes) and plot not in ["EGA_to_DTG", "calibration", "fit", "results"]:
            fig, ax = plt.subplots()

        match plot:
            case "EGA":
                plot_FTIR(self, ax, **kwargs)

            case "cumsum":
                temp = copy.deepcopy(self)
                temp.ega.update(self.ega.filter(self._info["gases"], axis=1).cumsum())
                plot_FTIR(temp, ax, **kwargs)

            case "TG":
                plot_TGA(self,"sample_mass", ax, **kwargs)

            case"heat_flow":
                if "heat_flow" in self.tga.columns:
                    plot_TGA(self, plot, ax, **kwargs)
                else:
                    logger.warning("No heat flow data available!")

            case"mass_steps":
                if "steps" not in kwargs:
                    self.mass_step(plot=False)
                    kwargs["steps"] = self.step_data().sample_temp.astype(np.float64)
                plot_mass_steps(self,ax,**kwargs)

            case "EGA_to_DTG":
                FTIR_to_DTG(self, save=save, **kwargs)

            case "fit":
                if self.results["fit"] == {}:
                    logger.warning(
                        'No fitting results available for plotting. Run ".fit()" first.'
                    )
                    return
                if reference not in self.results["fit"]:
                    logger.warning(
                        f"No fit performed with {reference}. Available options are {[self.results['fit'].keys()]}"
                    )
                    return
                else:
                    if len(self.results["fit"].keys()) == 1:
                        reference = list(self.results["fit"].keys())[0]
                    else:
                        avail_refs = [self.results['fit'].keys()]
                        warn_msg = f"Multiple fitting results available. Specify with keyword 'reference'= one of {avail_refs}"
                        logger.warning(warn_msg)
                        return
                plot_fit(self, reference, **kwargs)
                    
            case "calibration":
                if (gas := kwargs.get("gas")) and kwargs.get("gas") in self._info["gases"]:
                    x = self.xcali.get(gas)
                    y = self.ycali.get(gas)
                    linreg = self.linreg.loc[gas, :]

                    fig, ax = plt.subplots(1,2)
                    plot_calibration_single(x, y, linreg, ax[0])
                    plot_residuals_single(x, y, linreg, ax[1])
                else:
                    plot_calibration_combined(self.xcali, self.ycali, self.linreg, self.linreg.index)
            case "results":
                restype = kwargs.get("type", "fit")
                results = self.results[restype]
                if results:
                    results = pd.concat([
                        res
                        for res in results.values()
                        if isinstance(res, pd.DataFrame)
                    ])

                results = results.pint.convert_object_dtype().astype(np.float64).reset_index()
                avail_gases = results.gas.unique()
                avail_references = results.reference.unique()
                
                # validate inputs
                references = kwargs.get("references", avail_references)
                references = [ref for ref in references if ref in avail_references]
                gases = kwargs.get("gases", avail_gases)
                gases = [g for g in gases if g in avail_gases]
                if not references or not gases:
                    logger.error(f"There are no valid references or gases to plot.\nValid references are {avail_references!r}\nValid gases are {avail_gases!r}")
                    return
                # subset dataset
                data = results.query("gas in @gases & reference in @references")
                
                # validate plot arguments
                sig = signature(plot_results).parameters
                args = {key: kwargs.get(key) if kwargs.get(key) in get_args(sig[key].annotation)  else sig[key].default for key in sig.keys() if key!="data"} 
                fig, ax = plot_results(data, **args)

        if save:
            path_plots = PATHS["plots"]/ plot
            if not path_plots.exists():
                path_plots.mkdir(parents=True)

            path_pic = path_plots / f"{self.name}_{'_'.join(kwargs.keys())}"
            ax.get_figure().savefig(path_pic)
        return ax

    def fit(
        self,
        reference_name,
        T_max=None,
        T_max_tol=50,
        save=True,
        plot=True,
        presets=None,
        mod_sample=True,
        overwrite=False,
        **kwargs,
    ):
        "deconvolution of IR data"
        from ..fitting import fitting, get_presets
        if self.ega is None:
            logger.error("No EGA data to fit available.")
            return

        if reference_name in self.results["fit"]:
            if not overwrite:
                logger.warning(f"Fitting with{reference_name!r} already in results! Pass overwrite=True to redo fit.")
                results = self.results["fit"][reference_name]
            else:
                logger.warning(f"Fitting with{reference_name!r} was already in results! Overwriting.")

        if reference_name not in self.results["fit"] or overwrite:
            # setting upper limit for data
            if T_max is None:
                T_max = max(self.tga["sample_temp"]).magnitude
            elif T_max > (T_max_data := max(self.tga["sample_temp"])):
                logger.warning(f"{T_max=} exceeds maximum temperature of data ({T_max_data}).")
                T_max = max(self.tga["sample_temp"]).magnitude
                logger.info(f'"T_max" has been set to {T_max_data}.')

            # load presets for deconvolution
            if presets is None:
                presets = get_presets(reference_name)

            if presets is None:
                return

            for gas in presets:
                presets[gas] = presets[gas].drop(
                    presets[gas].index[presets[gas].loc[:, "center_0"] > T_max + T_max_tol]
                )

            # setting up output directory
            if save:
                path = (
                    PATHS["fitting"] / f'{general.time()}{reference_name}_{self._info["name"]}'
                )
                os.makedirs(path)
                os.chdir(path)

            # fitting
            logger.info(
                f'Fitting of "{self.name}" according to "{reference_name}" in Fitting_parameters.xlsx is in progress ...'
            )
            temp = copy.deepcopy(self)
            temp_unit = temp.tga.sample_temp.dtypes.units
            temp.tga = temp.tga[temp.tga["sample_temp"] < ureg.Quantity(T_max, temp_unit)]
            temp.ega = temp.ega[temp.ega["sample_temp"] < ureg.Quantity(T_max, temp_unit)]
            peaks = fitting(temp, presets, save=save, **kwargs)
            if save:
                logger.info(f"Plots and results are saved.\n'{path=}'.")
                os.chdir(PATHS["home"])

            logger.info("Fitting finished!")

            results = pd.concat(
                [peaks],
                keys=[(reference_name, self.sample, self.alias, self.run)],
                names=["reference", "sample", "alias", "run"],
            )
            if mod_sample:
                self.results["fit"].update({reference_name: results})

        # plotting
        if plot:
            self.plot("fit", reference=reference_name, **kwargs)

        return results

    def robustness(self, ref, **kwargs):
        from ..fitting import robustness
        data, summary = robustness(self, ref, **kwargs)
        self.results["robustness"].update({ref: {"data":data,"summary":summary}})
        return data, summary

    def save(self, how: Literal["samplelog", "excel", "pickle"]="samplelog", **kwargs):
        "save object or its contents as pickle file or excel"
        options = ["samplelog", "excel", "pickle"]
        if how not in options:
            logger.warning(f"{how=} not in {options=}")
            return

        # update samplelog
        samplelog(self.info.to_row(), create=True, **kwargs)
        match how:
            case "pickle":
                self.to_pickle()
            case "excel":
                self.to_excel() 
            case _:
                pass

    def to_pickle(self):
        if not PATHS["output"].exists():
            os.makedirs(PATHS["output"])
        with open(PATHS["output"] / f"{self.name}.pkl", "wb") as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def to_excel(self):
        if not PATHS["output"].exists():
            os.makedirs(PATHS["output"])

        path = PATHS["output"] / f'{self._info["name"]}.xlsx'
        with pd.ExcelWriter(path) as writer:
            try:
                pd.DataFrame.from_dict(
                    self._info.__dict__, orient="index"
                ).to_excel(writer, sheet_name="info")
            except PermissionError:
                logger.warning(
                    f"Unable to write on {path=} as the file is opened by another program."
                )
            for key in ["tga", 'ega']:
                try:
                    if self.__dict__[key] is not None:
                        self.__dict__[key].to_excel(writer, sheet_name=key)
                except PermissionError:
                    logger.warning(
                        f"Unable to write on {path=} as the file is opened by another program."
                    )
    def calibrate(self, profile=None,**kwargs):
        "calibrate object"
        from ..calibration import calibrate
        if not profile:
            profile=self.profile

        self.linreg, self.stats, self.xcali, self.ycali = calibrate(profile=profile, **kwargs)

    def __len__(self) -> int:
        return 1

    def __iter__(self):
        yield self



class Baseline(Sample):
    def __init_subclass__(cls) -> None:
        return super().__init_subclass__()

    def __repr__(self):
        attrs = ["name"]
        return f'Baseline({", ".join([f"{key}={repr(value)}" for key, value in self.__dict__.items() if key in attrs])})'

    def from_sample(self, sample:Sample, corrs: dict):
        self.name = "from data"
        self._info = SampleInfo(self.name, reference = None, alias = "Baseline", profile="from_sample", initial_mass = ureg.Quantity(0, "mg"), final_mass = ureg.Quantity(0, "mg")) 
        for name, methods in corrs.items():
            if not (sample.get(name) is not None and isinstance(sample.get(name), pd.DataFrame)):
                logger.warning(f"{name} is no correctable attribute of {sample.name}.")
                continue
            else:
                data = sample.__dict__[name].copy()
            all_signals = data.columns
            
            for pattern, method in methods.items():
                signals = [signal for signal in all_signals if re.match(pattern, signal)]

                if not signals:
                    logger.warning(f"{pattern!r} did not match any signal!")
                    continue

                kwargs = method.get("kwargs", {})
                match method.get("method"):
                    case "als":
                        fn = lambda x: corrections.baseline_als(x.astype(np.float64), **kwargs)
                    case "arpls":
                        fn = lambda x: corrections.baseline_arpls(x.astype(np.float64), **kwargs)
                    case "min":
                        fn = lambda x: np.ones_like(x) * x.min()
                    case "const":
                        fn= lambda x: corrections.baseline_const(x.astype(np.float64), **kwargs)
                    case "debug":
                        fn = lambda x: print(x)
                    case "none":
                        fn = lambda x: np.zeros_like(x)
                    case _:
                        fn = lambda x: np.zeros_like(x)

                data.update(data[signals].apply(fn))
            index_cols = sample.__dict__[name].filter(items=['time', 'sample_temp', 'reference_temp', 'sample_mass'])
            self.__dict__[name] = pd.concat([index_cols, data], axis=1)
        return self