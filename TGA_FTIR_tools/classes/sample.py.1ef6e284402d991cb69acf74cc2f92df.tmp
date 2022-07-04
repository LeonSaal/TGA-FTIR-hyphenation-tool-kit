import pandas as pd
from scipy.signal import savgol_filter
import copy
import os
import pickle
from ..config import SAVGOL, PATHS, COUPLING
from ..input_output import corrections, TGA, FTIR, general, samplelog, mass_step
import numpy as np
import scipy.stats as sp

from ..plotting import plot_fit, plot_dweight
import logging
import re


from dataclasses import dataclass, field
from typing import Optional

WINDOW_LENGTH = int(SAVGOL.getfloat("window_length"))
POLYORDER = int(SAVGOL.getfloat("POLYORDER"))


logger = logging.getLogger(__name__)


@dataclass
class Sample:
    name: str
    ir: Optional[pd.DataFrame] = field(default=None)
    tga: Optional[pd.DataFrame] = field(default=None)
    linreg: Optional[pd.DataFrame] = field(default=None)
    info: Optional[dict] = field(default_factory=dict)
    profile: str = COUPLING["profile"]
    alias: str = field(default=None)
    results: dict = field(default_factory=dict)

    def __post_init__(self, mode: str = "construct", **kwargs):
        if mode == "construct":
            logger.info(f'Initializing "{self.name}"')
            # load TG data
            self.tga = TGA.read_TGA(self.name, profile=self.profile)
            if self.tga is not None:
                logger.info("TGA data found.")
                self.tga["dtg"] = -savgol_filter(
                    self.tga["sample_mass"], WINDOW_LENGTH, POLYORDER, deriv=1
                )
                # deriving TG info
                try:
                    self.info = TGA.TGA_info(self.name, self.tga, profile=self.profile)
                except:
                    logger.info("Failed to derive TG info. Using default values.")
                    self.info = TGA.default_info(self.name, self.tga)
                try:
                    self.dry_weight(self, plot=False, **kwargs)
                except:
                    pass

            # load IR data
            self.ir = FTIR.read_FTIR(self.name)
            if self.ir is not None:
                from ..calibration import calibrate

                self.info["gases"] = set(self.ir.columns[1:].to_list())
                # load calibration
                self.linreg, self.stats = calibrate(mode="load")
                logger.info(
                    "IR data found{} for gases {}.".format(
                        " (* and calibrated) " if self.linreg is not None else "",
                        ", ".join(
                            [
                                gas
                                + (
                                    "*"
                                    if self.linreg is not None
                                    and gas in self.linreg.index
                                    else ""
                                )
                                for gas in self.info["gases"]
                            ]
                        ),
                    )
                )
                # derivinf IR info
                try:
                    self.info.update(FTIR.FTIR_info(self))
                except:
                    pass
                try:
                    self.ir["time"] = (
                        self.ir["time"] + 60 * self.info["background_delay"]
                    )
                    self.ir = pd.merge(
                        self.tga.filter(
                            ["time", "sample_temp", "reference_temp"], axis=1
                        ),
                        self.ir,
                        how="left",
                        on="time",
                    ).dropna(axis=0)
                except:
                    logger.info("No IR data found.")

                if self.ir is None and self.tga is None:
                    logger.error(f"Failed to initialize '{self.name}'.")
                else:
                    logger.info(f"'{self.name}' successfully initialiazed.")

            # assigning alias
            if not self.alias:
                if "alias" not in (log := samplelog(create=False)):
                    self.alias = log.loc[self.name, "alias"]
                else:
                    self.alias = self.name
            else:
                self.alias = str(self.alias)

            if match := re.match(r"(?P<sample>.+)_(?P<run>\d{3})$", self.name):
                self.sample = match.group("sample")
                self.run = match.group("run")
            elif match := re.match(r"(?P<sample>.+)_(?P<run>\d{3})$", self.alias):
                self.sample = match.group("sample")
                self.run = match.group("run")
            else:
                self.sample = self.alias
                self.run = 0

        # initialize object from pickle file
        if mode == "pickle":
            with open(os.path.join(PATHS["output"], self.name + ".pkl"), "rb") as inp:
                obj = pickle.load(inp)
            for key in obj.__dict__:
                self.__dict__[key] = obj.__dict__[key]

    def __repr__(self):
        attrs = ["name", "alias", "sample", "run"]
        return f'Sample({", ".join([f"{key}={repr(value)}" for key, value in self.__dict__.items() if key in attrs])})'

    def corr(self, reference=None, plot=False, **kwargs):
        "correction of TG and IR data"

        if "reference" in self.info:
            logger.warning(
                "Sample has already been corrected! Re-initialise object for correction."
            )
            return
        # try to load reference from samplelog if none is supplied
        if not reference:
            try:
                reference = samplelog(create=False).loc[self.info["name"], "reference"]
            except KeyError:
                logger.warning(
                    "No reference found in Samplelog. Please supply 'reference = '"
                )
                return

        # correction of data
        self.info["reference"] = reference
        if self.tga is not None:
            try:
                self.tga = corrections.corr_TGA(self.tga, reference, plot=plot)
                self.tga["dtg"] = -savgol_filter(
                    self.tga["sample_mass"], WINDOW_LENGTH, POLYORDER, deriv=1
                )
            except PermissionError:
                logger.error("Failed to correct TG data.")

            # filling TG_IR.info
            try:
                if self.info["reference_mass"] == "initial_mass":
                    # By default dry_weight() asumes how_dry = 'H2O'. If during initialization how_dry = None, this is catched here.
                    # kwargs = dict(kwargs, how_dry=None)
                    # However, there is no distinction between the other how_dry options (e.g. float), that still have to be passed to corr() again!
                    pass

                self.dry_weight(plot=plot, **kwargs)
                logger.info(
                    f'".info" of {self.info["name"]} was updated. To store these in Samplelog.xlsx run ".save()"'
                )
                success = True
            except PermissionError:
                logger.error("Failed to derive TG info.")

        if self.ir is not None:
            try:
                self.ir.update(corrections.corr_FTIR(self, reference, plot=plot))
            except:
                logger.error("Failed to correct IR data.")

            try:
                self.info.update(FTIR.FTIR_info(self))
                if not success:
                    logger.info(
                        "'.info' was updated. To store these in Samplelog.xlsx run '.save()'"
                    )
            except:
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

        except:
            logger.error("Failed to derive TG info.")

    def mass_step(
        self, plot=True, height=0, width=100, prominence=0, rel_height=0.9, **kwargs
    ):
        w_T = width / sp.mode(np.diff(self.tga.reference_temp)).mode[0]
        return mass_step(
            self.tga,
            plot=plot,
            height=height,
            width=w_T,
            prominence=prominence,
            rel_height=rel_height,
            **kwargs,
        )

    def plot(self, plot, **kwargs):
        from ..plotting import plot_FTIR, plot_TGA, FTIR_to_DTG, plot_fit

        "plotting TG and or IR data"

        options = ["TG", "heat_flow", "IR", "DIR", "cumsum", "IR_to_DTG", "fit"]
        if plot not in options:
            logger.warn(f"{plot} not in supported {options=}.")

        if self.ir is None and (plot in ["IR", "DIR", "cumsum", "IR_to_DTG"]):
            logger.warn("Option unavailable without IR data.")
            return
        else:
            if plot == "IR":
                plot_FTIR(self, **kwargs)

            if plot == "DIR":
                temp = copy.deepcopy(self)
                temp.ir.update(
                    self.ir.filter(self.info["gases"], axis=1)
                    .diff()
                    .ewm(span=10)
                    .mean()
                )
                plot_FTIR(temp, **kwargs)

            if plot == "cumsum":
                temp = copy.deepcopy(self)
                temp.ir.update(self.ir.filter(self.info["gases"], axis=1).cumsum())
                plot_FTIR(temp, **kwargs)

        if (self.tga is None) and (plot in ["TG", "heat_flow", "IR_to_DTG"]):
            logger.warn("Option unavailable without TGA data.")
            return
        else:
            if plot == "TG":
                plot_TGA(self, "sample_mass", **kwargs)

            if plot == "heat_flow":
                if "heat_flow" in self.tga.columns:
                    plot_TGA(self, plot, **kwargs)
                else:
                    logger.warn("No heat flow data available!")

        if (self.linreg is not None) and (plot == "IR_to_DTG"):
            FTIR_to_DTG(self, **kwargs)
        elif (self.linreg is None) and (plot == "IR_to_DTG"):
            logger.warn("Option unavailable without calibration!")
            return

        if plot == "fit":
            if "fit" in self.results:
                plot_fit(self, **kwargs)
            else:
                logger.warn(
                    'No fitting results available for plotting. Run ".fit()" first.'
                )
                
    def fit(
        self,
        reference,
        T_max=None,
        T_max_tol=50,
        save=True,
        plot=True,
        presets=None,
        mod_sample=True,
        **kwargs,
    ):
        "deconvolution of IR data"
        from ..fitting import fitting, get_presets

        # setting upper limit for data
        if T_max is None:
            T_max = max(self.tga["sample_temp"])
        elif T_max > (T_max_data := max(self.tga["sample_temp"])):
            logger.warn(f"{T_max=} exceeds maximum temperature of data ({T_max_data}).")
            T_max = max(self.tga["sample_temp"])
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
            path = os.path.join(
                PATHS["fitting"], general.time() + reference + "_" + self.info["name"],
            ).replace(os.sep, os.altsep)
            os.makedirs(path)
            os.chdir(path)

        # fitting
        logger.info(
            f'Fitting of "{self.name}" according to "{reference}" in Fitting_parameters.xlsx is in progress ...'
        )
        temp = copy.deepcopy(self)
        temp.tga = temp.tga[temp.tga["sample_temp"] < T_max]
        temp.ir = temp.ir[temp.ir["sample_temp"] < T_max]
        peaks = fitting(temp, presets, save=save, **kwargs)
        if save:
            logger.info(f"Plots and results are saved.\n'{path=}'.")
            os.chdir(PATHS["home"])

        logger.info("Fitting finished!")
        if mod_sample:
            self.results["fit"] = peaks

        # plotting
        if plot:
            self.plot("fit", **kwargs)

        return pd.concat(
            [peaks], keys=[(self.sample, self.run)], names=["sample", "run"]
        )
    def robustness():
        pass

    def save(self, how="samplelog", **kwargs):
        "save object or its contents as pickle file or excel"

        # update samplelog
        samplelog(self.info.__dict__, create=True, **kwargs)
        path_output = PATHS["output"]
        if how == "samplelog":
            return

        if os.path.exists(path_output) == False:
            os.makedirs(path_output)
        # save object
        if how == "pickle":
            with open(
                os.path.join(path_output, self.info["name"] + ".pkl"), "wb"
            ) as output:
                pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        elif how == "excel":
            path = os.path.join(path_output, self.info["name"] + ".xlsx")
            with pd.ExcelWriter(path) as writer:
                try:
                    pd.DataFrame.from_dict(self.info.__dict__, orient="index").to_excel(
                        writer, sheet_name="info"
                    )
                except PermissionError:
                    logger.warn(
                        f"Unable to write on {path=} as thre file is opened by another program."
                    )
                for key in ["tga", "ir"]:
                    try:
                        self.__dict__[key].to_excel(writer, sheet_name=key)
                    except PermissionError:
                        logger.warn(
                            f"Unable to write on {path=} as thre file is opened by another program."
                        )

    def calibrate(self, **kwargs):
        "calibrate object"
        from ..calibration import calibrate

        self.linreg, self.stats = calibrate(**kwargs)
