import copy
import logging
import os
import re
import shutil as sh

import numpy as np
import pandas as pd
import scipy as sp

from ..config import BOUNDS, PATHS, config
from ..input_output.general import time
from ..input_output.corrections import baseline_als
from ..utils import multi_gauss
from .fitdata import FitData
import pint_pandas
import pint
ureg = pint.get_application_registry()

logger = logging.getLogger(__name__)

url_fitting_param = "https://raw.githubusercontent.com/BAMresearch/TGA-data.ir-hyphenation-tool-kit/9382aaea97048e507bdc56715f971a9dec25be6d/Fitting_parameter.xlsx"





def fitting(
    sample, presets, func=multi_gauss, y_axis="orig", save=True, predef_tol=0.01,
):
    from ..calibration import calibration_stats
    "deconvolve IR data of sample with func and multiple presets"
    # temp_presets = copy.deepcopy(presets)

    # initializing output DataFrame
    # extract links between groups and sort gases accordingly
    # thresholds for fit parameters

    data = FitData(sample, copy.deepcopy(presets))

    # cycling through gases
    for gas in data.gases:
        # correction of water-signal drift
        if gas == "H2O":
            data.ega[gas] -= baseline_als(data.ega[gas].astype(np.float64))

        # molar desorption
        tot_area = np.sum(data.ega[gas])

        if gas == "H2O":
            # total area of water is calculated above dry-point, if this exists (therefore, try)
            if "dry_temp" in data.info:
                minimum_temp = data.info["dry_temp"]
            else:
                minimum_temp = data.ega[gas].min()
            tot_area = np.sum(data.ega[gas][data.ega.sample_temp > minimum_temp])

        if y_axis == "rel" and data.linreg is not None:
            tot_mol = (tot_area - data.linreg.intercept[gas]) / data.linreg.slope[gas]
            data.ega.update(data.ega[gas] / tot_area * tot_mol)

        # predefining guesses and bounds
        params_0, params_min, params_max = data.get_bounds(gas, sample, predef_tol)

        # actual fitting
        x, y = data.ega.sample_temp.to_numpy(dtype=np.float64), data.ega[gas].to_numpy(dtype=np.float64)
        try:
            popt, _ = sp.optimize.curve_fit(
                func, x, y, p0=params_0, bounds=(params_min, params_max)
            )
        except ValueError:
            logger.error(f"Failed to fit {gas} signal")
            break

        # return values
        data.update_peaks(popt, gas, (x, y))

    if data.linreg is not None:
        data.calculate_molar_mass()

    # calculate summarized groups
    #data.summarize()

    #if sample.stats is not None:
    #    check_LODQ_frame(data.peaks, sample.stats)

    # save results to excel
    if save:
        data.save(y_axis)

    return data.peaks

def get_presets(reference, file = config["fitting_params"]):
    "load deconvolution presets from excel file"
    # load raw data from file
    presets = dict()

    if not os.path.exists(file):
        if PATHS["fitting_params"].exists():
            sh.copy(PATHS["fitting_params"], file)
        else:
            logger.warning("Unable to get default settings.")

    references = pd.read_excel(
        file, index_col=0, header=[0, 1], sheet_name=None,
    )

    if (
        reference not in (options := references["center_0"].index.to_list())
        and reference is not None
    ):
        logger.warning(f"{reference=} is an invalid option. {options=} ")
        return
    elif reference is None:
        return options

    gases = list(set(references["center_0"].columns.get_level_values(1)))
    # organizing data in dict, sorted by gases and filling in missing values with [fitting] parameters of settings.ini

    for gas in gases:
        index = [
            group
            for group, group_gas in references["center_0"].columns
            if group_gas == gas
        ]
        data = pd.DataFrame(index=index)
        for key in references:
            data[key] = pd.DataFrame(
                references[key]
                .loc[reference, :][references[key].columns.get_level_values(1) == gas]
                .T.values,
                index=index,
                columns=[key],
            )  # .dropna(axis=1)
        presets[gas] = data.dropna(axis=0, how="all")

        params = [
            "height_0",
            "hwhm_0",
            "center_min",
            "hwhm_min",
            "height_min",
            "center_max",
            "hwhm_max",
            "height_max",
            "link",
        ]
        vals = [
            BOUNDS.getfloat("height_0"),
            BOUNDS.getfloat("hwhm_0"),
            pd.Series(presets[gas].loc[:, "center_0"] - BOUNDS.getfloat("tol_center")),
            BOUNDS.getfloat("hwhm_min"),
            BOUNDS.getfloat("height_min"),
            pd.Series(presets[gas].loc[:, "center_0"] + BOUNDS.getfloat("tol_center")),
            BOUNDS.getfloat("hwhm_max"),
            BOUNDS.getfloat("height_max"),
            "0",
        ]
        infill = dict(zip(params, vals))
        presets[gas] = presets[gas].fillna(infill).dropna()
        if presets[gas].empty:
            del presets[gas]

    return presets


def check_LODQ(val: float, gas: str, *, stats: pd.DataFrame):
    for limit in ["x_LOD", "x_LOQ"]:
        if val < stats.loc[gas, limit]:
            return f"< {limit[2:]}"

    return "> LOQ"


def check_LODQ_frame(peaks: pd.DataFrame, stats: pd.DataFrame) -> None:
    total_cond = peaks.index.get_level_values("group") == "total"
    total = peaks[total_cond].dropna(subset=["area"])
    peaks.loc[("total", slice(None)), "limits"] = total.apply(
        lambda x: check_LODQ(x.mmol, x.name[1], stats=stats), axis=1
    )
    for gas, _ in peaks.groupby("gas"):
        peaks["limits"] = peaks.loc[("total", gas), "limits"]

