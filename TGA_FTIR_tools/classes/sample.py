import copy
import json
import logging
import os
import pickle
import re
from dataclasses import InitVar, dataclass, field
from types import NoneType
from typing import Any, Literal, Mapping, Optional, List
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import pint_pandas

import pint
from pint import DimensionalityError
ureg = pint.get_application_registry()
import scipy as sp
from scipy.signal import savgol_filter

from ..config import DEFAULTS, PATHS, SAVGOL, update_config
from ..input_output import (FTIR, TGA, corrections, general, mass_step,
                            read_data, samplelog, read_profile_json)
from ..plotting import plot_dweight, plot_mass_steps, plot_calibration_single, plot_residuals_single, plot_calibration_combined
from ..utils import select_import_profile, check_profile_exists
from .info import SampleInfo

WINDOW_LENGTH_REL = SAVGOL.getfloat("window_length_rel")
POLYORDER = int(SAVGOL.getfloat("polyorder"))


logger = logging.getLogger(__name__)


@dataclass
class Sample:
    name: str
    ega: Optional[pd.DataFrame] = field(default=None)
    tga: Optional[pd.DataFrame] = field(default=None)
    linreg: Optional[pd.DataFrame] = field(default=None)
    alias: str = field(default=None)
    reference: str = field(default=None)
    _info = None 
    baseline = None
    results: dict = field(default_factory=lambda: {"fit": {}, "robustness": {}})
    mode: Literal["construct", "pickle"] = "construct"
    profile: str = DEFAULTS["profile"]

    def __post_init__(self, **kwargs):
        if self.mode == "construct":
            logger.info(f"Initializing {self.name!r} with {self.profile!r}")

            self.check_profile()
            self.profile_data = read_profile_json(self.profile)
            logger.debug(f"Using profile {self.profile!r}.")

            # load data
            self.__dict__.update(read_data(self.name, profile=self.profile))
            self._info = SampleInfo(self.name, info=self._info)

            if self.tga is not None:
                logger.info("TGA data found.")
                logger.debug("Checking required columns in TGA data.")
                
                # deriving TG info
                # try:
            
                # except Exception as e:
                #     logger.info(f"Failed to derive TG info. Using default values. {e}")
                #     self._info = TGA.default_info(self.name, self.tga)

                if not self._info["initial_mass"] and "sample_mass" in self.tga.columns:
                    self._info["initial_mass"]=self.tga["sample_mass"].iloc[0]

                # calculate sample mass with mass loss and initial mass
                try:
                    if "sample_mass" not in self.tga.columns and "mass_loss" in self.tga.columns:
                        if "initial_mass" in self._info:
                            self.tga["sample_mass"] =  self.tga["mass_loss"] + self._info["initial_mass"]
                except DimensionalityError:
                    logger.error("Failed")
                self._info.steps_idx.update({})

                # check required columns
                if missing:=self.missing_tga_columns():
                    logger.error(f"Required column(s) {missing} not found in TGA data of {self.name!r}")
                else:   
                    try:
                        window_length = int(self.tga.index.size * WINDOW_LENGTH_REL)
                        dtg = savgol_filter(
                            self.tga["sample_mass"].values._data, window_length if window_length%2 ==0 else window_length+1 , POLYORDER, deriv=1
                        )
                        dtype =  self.tga["sample_mass"].dtype
                        self.tga = self.tga.assign(dtg = -pd.Series(dtg).astype(dtype))
                    except ValueError as e:
                        print(e)


            else:
                logger.error(f"No TGA data found for {self.name!r}.")
                self._info = SampleInfo(name=self.name, alias=self.alias)

            if self.ega is not None:
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

                # deriving IR info
                try:
                    self._info.update(FTIR.FTIR_info(self))
                except Exception as e:
                    logger.error(f"Unable to derive EGA-info. {e}")
                    #raise e
                
                if "sample_temp" not in self.ega and "time" in self.ega:
                    self.ega.time = self.ega.time.astype(self.tga.time.dtype)
                    self.ega = self.ega.set_index("time").reindex(self.tga.time, method="nearest").reset_index()
                    try:
                        self.ega = pd.merge(
                            self.tga.filter(
                                ["time", "sample_temp", "reference_temp"], axis=1
                            ),
                            self.ega,
                            how="left",
                            on="time",
                        )#.dropna(axis=0)
                    except Exception as e:
                        logger.info(f"Unable to merge temperature data. {e}")
            try:
                self.dry_weight(plot=False, **kwargs)
            except Exception as e:
                logger.error(f"Failed to derive dry weight: {e}")

            if self.ega is None and self.tga is None:
                logger.error(f"Failed to initialize '{self.name}'.")
            else:
                logger.info(f"'{self.name}' successfully initialiazed.")

            # assigning alias
            if not self.alias:
                log  = samplelog(create=False)
                if self.name in log.index:
                    self.alias = log.loc[self.name, "alias"]
                else:
                    self.alias = self.name
            else:
                self.alias = str(self.alias)
            self._info.alias = self.alias

            if match := re.match(r"(?P<sample>.+)_(?P<run>\d{,3})$", self.name):
                self.sample = match.group("sample")
                self.run = match.group("run")
            elif match := re.match(r"(?P<sample>.+)_(?P<run>\d{,3})$", self.alias):
                self.sample = match.group("sample")
                self.run = match.group("run")
            else:
                self.sample = self.alias
                self.run = 0

        # initialize object from pickle file
        if self.mode == "pickle":
            with open(PATHS["output"] / f"{self.name}.pkl", "rb") as inp:
                obj = pickle.load(inp)
            for key in obj.__dict__:
                self.__dict__[key] = obj.__dict__[key]
        self.raw = copy.deepcopy(self)

    def __repr__(self):
        attr_names = ["name", "alias", "sample", "reference", "run"]
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
            logger.warning("Can only add to Sample- or Worklist-type")

    @property
    def info(self):
        self._info.alias = self.alias
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

    def corr(
        self,
        reference: NoneType | str = None,
        plot: Mapping | bool = False,
        update=False,
        dry_args={},
        ega_args={},
    ):
        "correction of TG and IR data"
        if self.ega is None and self.tga is None:
            logger.error("There is no data to correct.")
            return

        if plot == True:
            plot = {}
            plot["ega"] = True
            plot["tga"] = True
            plot["dry_weight"] = True
        if plot == False:
            plot = {}
            plot["ega"] = False
            plot["tga"] = False
            plot["dry_weight"] = False

        if self.reference and not update:
            logger.warning(
                "Sample has already been corrected! Re-initialise object for correction or specify 'update'=True."
            )
            return
        # try to load reference from samplelog if none is supplied
        if not reference:
            if self.reference:
                reference = self.reference
            else:
                if self.name in (log := samplelog(create=False).reference).index:
                    reference = log.loc[self.name]
                else:
                    logger.warning(
                        "No reference found in Samplelog. Please supply 'reference = '"
                    )
                    return

        # correction of data
        self.baseline = Baseline(reference)
        self.reference = reference
        if self.tga is not None:
            try:
                self.tga = corrections.corr_TGA_Baseline(
                    self.raw.tga, self.baseline, plot=plot["tga"]
                )
                self.tga["dtg"] = -savgol_filter(
                    self.tga.sample_mass, WINDOW_LENGTH_REL, POLYORDER, deriv=1
                )
            except PermissionError:
                logger.error("Failed to correct TG data.")

            # filling Sample.info
            try:
                if self._info['reference_mass_name'] == "initial_mass":
                    # By default dry_weight() asumes how_dry = 'H2O'. If during initialization how_dry = None, this is catched here.
                    # kwargs = dict(kwargs, how_dry=None)
                    # However, there is no distinction between the other how_dry options (e.g. float), that still have to be passed to corr() again!
                    pass

                self.dry_weight(plot=plot["dry_weight"], **dry_args)
                logger.info(
                    f'".info" of {self._info["name"]} was updated. To store these in Samplelog.xlsx run ".save()"'
                )
                success = True
            except PermissionError:
                logger.error("Failed to derive TG info.")

        if (self.ega is not None) and (self.baseline.ega is not None):
            try:
                self.ega.update(
                    corrections.corr_FTIR(
                        self.raw, self.baseline, plot=plot['ega'], **ega_args
                    )
                )
            except PermissionError:
                logger.error("Failed to correct IR data.")

            try:
                self._info.update(FTIR.FTIR_info(self))
                if not success:
                    logger.info(
                        "'.info' was updated. To store these in Samplelog.xlsx run '.save()'"
                    )
            except PermissionError:
                logger.error("Failed to derive IR info.")

    def get_value(self, *values, which="sample_mass", at="sample_temp"):
        "extract values from TG data at e.g. certain temperatures"

        out = pd.DataFrame(index=[which], columns=pd.Index(values, name=at))
        for value in values:
            out.loc[which, value] = self.tga[which][self.tga[at] >= value].values[0]

        return out

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
            raise e
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
                steps=self.step_data().sample_temp,
                y_axis="orig",
            )
        return step_height, rel_step_height, step_starts_idx, step_ends_idx, peaks_idx

    def plot(self, plot=Literal["TG","mass_steps","heat_flow","EGA", "DIR","cumsum","IR_to_DTG","fit", "calibration"], ax=None, save=False,reference=None, **kwargs):
        from ..plotting import FTIR_to_DTG, plot_fit, plot_FTIR, plot_TGA

        "plotting TG and or IR data"
        options = [
            "TG",
            "mass_steps",
            "heat_flow",
            "EGA",
            "DIR",
            "cumsum",
            "IR_to_DTG",
            "fit",
            "calibration"
        ]
        if plot not in options:
            logger.warning(f"{plot} not in supported {options=}.")

        needs_tg = ["TG", "heat_flow", "IR_to_DTG", "mass_steps"]
        needs_ega = ["EGA", "DIR", "cumsum", "IR_to_DTG"]
        needs_cali = [ "IR_to_DTG","calibration"]
        
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

        if not isinstance(ax, plt.Axes) and plot not in ["IR_to_DTG", "calibration", "fit"]:
            fig, ax = plt.subplots()

        match plot:
            case "EGA":
                plot_FTIR(self, ax, **kwargs)

            case "DIR":
                temp = copy.deepcopy(self)
                temp.ega.update(
                    self.ega.filter(self._info["gases"], axis=1).diff().ewm(span=10).mean()
                )
                plot_FTIR(temp, ax, **kwargs)

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
                    kwargs["steps"] = self._info.step_temp
                plot_mass_steps(self,ax,**kwargs)

            case "IR_to_DTG":
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

                    fig, axs = plt.subplots(1,2)
                    plot_calibration_single(x, y, linreg, axs[0])
                    plot_residuals_single(x, y, linreg, axs[1])
                else:
                    plot_calibration_combined(self.xcali, self.ycali, self.linreg, self._info["gases"])

        if save:
            path_plots = PATHS["plots"]/ plot
            if not path_plots.exists():
                path_plots.mkdir(parents=True)

            path_pic = path_plots / f"{self.name}_{'_'.join(kwargs.keys())}"
            ax.get_figure().savefig(path_pic)
        return ax

    def fit(
        self,
        reference,
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

        if reference in self.results["fit"]:
            if not overwrite:
                logger.warning(f"Fitting with{reference!r} already in results! Pass overwrite=True to redo fit.")
                results = self.results["fit"][reference]
            else:
                logger.warning(f"Fitting with{reference!r} was already in results! Overwriting.")

        if reference not in self.results["fit"] or overwrite:
            # setting upper limit for data
            if T_max is None:
                T_max = max(self.tga["sample_temp"]).magnitude
            elif T_max > (T_max_data := max(self.tga["sample_temp"])):
                logger.warning(f"{T_max=} exceeds maximum temperature of data ({T_max_data}).")
                T_max = max(self.tga["sample_temp"]).magnitude
                logger.info(f'"T_max" has been set to {T_max_data}.')

            # load presets for deconvolution
            if presets is None:
                presets = get_presets(reference)

            if presets is None:
                return

            for gas in presets:
                presets[gas] = presets[gas].drop(
                    presets[gas].index[presets[gas].loc[:, "center_0"] > T_max + T_max_tol]
                )

            # setting up output directory
            if save:
                path = (
                    PATHS["fitting"] / f'{general.time()}{reference}_{self._info["name"]}'
                )
                os.makedirs(path)
                os.chdir(path)

            # fitting
            logger.info(
                f'Fitting of "{self.name}" according to "{reference}" in Fitting_parameters.xlsx is in progress ...'
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
                keys=[(reference, self.sample, self.alias, self.run)],
                names=["reference", "sample", "alias", "run"],
            )
            if mod_sample:
                self.results["fit"].update({reference: results})

        # plotting
        if plot:
            self.plot("fit", reference=reference, **kwargs)

        return results

    def robustness(self, ref, **kwargs):
        from ..fitting import robustness
        data, summary = robustness(self, ref, **kwargs)
        self.results["robustness"].update({ref: {"data":data,"summary":summary}})
        return data, summary

    def save(self, how: Literal["samplelog", "excel", "pickle"], **kwargs):
        "save object or its contents as pickle file or excel"
        options = ["samplelog", "excel", "pickle"]
        if how not in options:
            logger.warning(f"{how=} not in {options=}")
            return

        # update samplelog
        info = {key: str(value) for key, value in info.items() if value is not None}
        data = pd.DataFrame.from_dict(info, orient="index").T
        data.set_index("name", inplace=True)
        samplelog(data, create=True, **kwargs)
        path_output = PATHS["output"]
        if how == "samplelog":
            return

        if not path_output.exists():
            os.makedirs(path_output)

        # save object
        if how == "pickle":
            with open(path_output / f"{self.name}.pkl", "wb") as output:
                pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        elif how == "excel":
            path = path_output / f'{self._info["name"]}.xlsx'
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
