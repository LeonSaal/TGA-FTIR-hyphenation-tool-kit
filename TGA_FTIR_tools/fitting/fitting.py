import copy
import logging
import os
import re
import shutil as sh

import numpy as np
import pandas as pd
import scipy as sp

from ..config import BOUNDS, COUPLING, MERGE_CELLS, PATHS, config
from ..input_output.general import time
from ..utils import gaussian, multi_gauss
from .fitdata import FitData

logger = logging.getLogger(__name__)

url_fitting_param = "https://raw.githubusercontent.com/BAMresearch/TGA-data.ir-hyphenation-tool-kit/9382aaea97048e507bdc56715f971a9dec25be6d/Fitting_parameter.xlsx"


def baseline_als(
    y, lam=1e6, p=0.01, niter=10
):  # https://stackoverflow.com/questions/29156532/python-baseline-correction-library
    L = len(y)
    D = sp.sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sp.sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sp.sparse.linalg.spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def fitting(
    sample, presets, func=multi_gauss, y_axis="orig", save=True, predef_tol=0.01,
):
    from ..calibration import stats
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
            data.ir[gas] -= baseline_als(data.ir[gas])

        # molar desorption
        tot_area = np.sum(data.ir[gas])

        if gas == "H2O":
            # total area of water is calculated above dry-point, if this exists (therefore, try)
            if "dry_temp" in data.info:
                minimum_temp = data.info["dry_temp"]
            else:
                minimum_temp = data.ir[gas].min()
            tot_area = np.sum(data.ir[gas][data.ir.sample_temp > minimum_temp])

        if y_axis == "rel" and data.linreg is not None:
            tot_mol = (tot_area - data.linreg.intercept[gas]) / data.linreg.slope[gas]
            data.ir.update(data.ir[gas] / tot_area * tot_mol)

        # predefining guesses and bounds
        params_0, params_min, params_max = data.get_bounds(gas, sample, predef_tol)

        # actual fitting
        x, y = data.ir.sample_temp, data.ir[gas]
        try:
            popt, _ = sp.optimize.curve_fit(
                func, x, data.ir[gas], p0=params_0, bounds=(params_min, params_max)
            )
        except ValueError:
            logger.error(f"Failed to fit {gas} signal")
            break

        # return values
        data.update_peaks(popt, gas, (x, y))

    if data.linreg is not None:
        data.calculate_molar_mass()

    # calculate summarized groups
    data.summarize()

    if stats is not None:
        check_LODQ_frame(data.peaks, stats)

    # save results to excel
    if save:
        data.save(y_axis)

    return data.peaks


def fits(worklist, reference, save=True, presets=None, mod_sample=True, **kwargs):
    "perform decovolution on multiple sample objects"
    from ..calibration import stats

    # load default presets
    if not presets:
        presets = get_presets(reference)

    if presets is None:
        return

    # make subdirectory to save data
    if save:
        sample_names = "".join(
            list(set([str(sample.info["name"]) for sample in worklist]))
        )
        sample_names = "".join(
            [x if (x.isalnum() or x in "._- ") else "" for x in str(sample_names)]
        )  # to catch invalide sample names
        path = os.path.join(
            PATHS["fitting"], time() + reference + "_" + "_" + sample_names
        ).replace(os.sep, os.altsep)
        # check path length and if necessary shorten file name by list of samples, regarding expacted .png files to be saved to this directory
        longest_name_length = 0
        for sample in worklist:
            obj_name_length = len(sample.info["name"])
            if obj_name_length > longest_name_length:
                longest_name_length = obj_name_length
        if (len(path) + longest_name_length) > 258:
            sample_names = sample_names[
                : (len(sample_names) - ((len(path) + longest_name_length) - (258 - 8)))
            ]  # -8 for _gas and .png
            path = os.path.join(
                PATHS["fitting"], time() + reference + "_" + "_" + sample_names
            ).replace(os.sep, os.altsep)
        os.makedirs(path)
        os.chdir(path)

    results = {}
    # cycling through samples
    for sample in worklist:
        # writing data to output DataFrames
        results[sample.alias] = sample.fit(
            reference, presets=presets, mod_sample=mod_sample, **kwargs, save=False
        )

    # return results
    collection = pd.concat([item for item in results.values()])

    statistics = []
    dm = COUPLING.getfloat("mass_resolution")

    # group by samples
    names = ["run", "group"]
    for sample, groups in collection.groupby("sample"):
        sample_stats = []
        # group by desorbed group
        for (group, gas), values in groups.drop("limits", axis=1).groupby(
            ["group", "gas"]
        ):
            if gas in ["stddev"]:
                continue
            # calculate statistical values
            mean = pd.concat(
                [values.mean().to_frame(gas).T], keys=[("mean", group)], names=names,
            )
            dev = pd.concat(
                [values.std().to_frame(gas).T], keys=[("stddev", group)], names=names,
            )
            mean.index.rename(["run", "group", "gas"], inplace=True)
            dev.index.rename(["run", "group", "gas"], inplace=True)

            # calculate error propagation
            if stats is not None and gas in stats.index:
                lod = stats["x_LOD"][gas]
                mmol_p_mg = values["mmol_per_mg"]
                mmol = values["mmol"]
                mg = mmol / mmol_p_mg

                dmmol_mg_i = (
                    np.power(np.power(lod / mmol, 2) + np.power(dm / mg, 2), 0.5)
                    * mmol_p_mg
                )
                dmmol = np.power(np.sum(np.power(dmmol_mg_i, 2)), 0.5)

                err_dev = pd.concat(
                    [pd.DataFrame([dmmol], columns=["mmol_per_mg"], index=[gas])],
                    keys=[("err_dev", group)],
                    names=["run", "group"],
                )
                err_dev.index.rename(["run", "group", "gas"], inplace=True)
                dev = pd.concat([dev, err_dev])

            rel_dev = dev / mean.values
            rel_dev.rename(
                {old: f"rel_{old}" for old in rel_dev.index.get_level_values("run")},
                inplace=True,
            )
            sample_stats.append(
                pd.concat([pd.concat([mean, dev, rel_dev])], keys=[sample])
            )

        stats_df = pd.concat(sample_stats)
        statistics.append(stats_df)

    results = pd.concat([collection, pd.concat(statistics)])
    results.sort_index(level=["sample", "group"], inplace=True)

    # check results for LOD and LOQ and add columns 'rel_dev', limits

    # exporting data
    if save:
        logger.info(f"Fitting finished! Results are saved in {path=}.")
        with pd.ExcelWriter("summary.xlsx") as writer:
            results.to_excel(writer, sheet_name="summary", merge_cells=MERGE_CELLS)
        os.chdir(PATHS["home"])
    if mod_sample:
        worklist.results["fit"] = results
    return results.drop("total", level="group")


def get_presets(reference):
    "load deconvolution presets from excel file"
    # load raw data from file
    presets = dict()

    if not os.path.exists(config["fitting_params"]):
        if os.path.exists(PATHS["fitting_params"]):
            sh.copy(PATHS["fitting_params"], config["fitting_params"])
        else:
            logger.warn("Unable to get default settings.")

    references = pd.read_excel(
        config["fitting_params"], index_col=0, header=[0, 1], sheet_name=None,
    )

    if (
        reference not in (options := references["center_0"].index.to_list())
        and reference is not None
    ):
        logger.warn(f"{reference=} is an invalid option. {options=} ")
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

