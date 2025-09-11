import logging
import os
from pathlib import Path

from matplotlib.image import interpolations_names
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from .utils import validate_mf
from sklearn import linear_model
import re
from .config import MERGE_CELLS, PATHS, SEP, UNITS
from .plotting import plot_integration, plot_calibration_single, plot_calibration_combined, plot_residuals_single
import pint
from .input_output.general import time
from molmass import Formula
ureg = pint.get_application_registry()

logger = logging.getLogger(__name__)

REGR_COLS = ["molecular_formula","molmass", "slope", "intercept", "r_value", "p_value", "std_error"]
FIGSIZE = np.array(plt.rcParams["figure.figsize"])

def integrate_peaks(
    ega_data, step_starts_idx, step_ends_idx, peaks_idx, corr_baseline=None, plot=False, ax=None, gases=None
):
    "integrating IR signal in between given bounds"
    gases = list(gases)
    integrals = pd.DataFrame(index=range(len(step_starts_idx)), columns=gases)

    # integration
    baselines = {}
    for gas in gases:
        baselines[gas] = {}
        for i in range(len(step_ends_idx)):
            subset = ega_data[gas].iloc[step_starts_idx[i]:step_ends_idx[i]]

            # baseline correction
            if corr_baseline == "linear":
                baseline = np.linspace(
                    subset.iloc[0], subset.iloc[len(subset) - 1], len(subset)
                )
            elif corr_baseline == "const":
                baseline = min(subset.iloc[0], subset.iloc[len(subset) - 1]) * np.ones(
                    len(subset)
                )  # np.linspace(min(subset),min(subset),len(subset))
            elif corr_baseline == None:
                baseline = np.zeros(len(subset))

            integral = sp.integrate.simpson(subset - baseline)
            integrals.loc[i, gas] = integral
            baselines[gas][i] = baseline

    if plot:
        plot_integration(ega_data, baselines, peaks_idx, step_starts_idx, step_ends_idx, gases, ax)

    return integrals


def eval_lin(x, slope, intercept):
    "evaluating linear equation with slope and intercept at ponts x"
    return slope * x + intercept


def calibration_stats(x_cali, y_cali, linreg, alpha=0.95, beta=None, m=1, k=3):
    "calculating calibration stats according to DIN 32645"
    gases = linreg.index

    n = len(x_cali)
    f = n - 2

    if not beta:
        beta = alpha

    stats_rows = []
    for gas in gases:
        x_unit = x_cali[gas].dtype.units
        y_unit = y_cali[gas].dtype.units
        b = ureg.Quantity(linreg["slope"][gas], y_unit / x_unit)
        a = ureg.Quantity(linreg["intercept"][gas], y_unit)

        s_yx = np.sqrt(np.sum(np.power(b * x_cali[gas] + a - y_cali[gas], 2)) / (n - 2))
        s_x0 = s_yx / b
        x_ = np.mean(x_cali[gas])
        Q_x = np.sum(np.power(x_cali[gas] - x_, 2))

        x_NG = (
            s_x0 * sp.stats.t.ppf(alpha, f) * np.sqrt(1 / m + 1 / n + (x_ * x_) / Q_x)
        )
        # x_EG=x_NG+s_x0*sp.stats.t.ppf(beta,f)*np.sqrt(1/m+1/n+(x_*x_)/Q_x)
        x_BG = k * x_NG

        stats_row = pd.DataFrame(
            [[s_yx / y_unit, s_x0, x_NG, x_BG]],
            index=[gas],
            columns=["s_yx", "s_x0", "x_LOD", "x_LOQ"]
        )
        stats_rows.append(stats_row)

    return pd.concat(stats_rows).pint.convert_object_dtype()


def calibrate(worklist=None, molecular_formulas = {},plot=False, mode="load", method="max", profile=None, width_T=np.array([0, np.inf]), min_rel_height = .2, corr_baseline="linear", **fig_args):
    methods = ['max', "iter", "co_oxi", "co_oxi_iter", 'mlr']
    if method not in methods:
        logger.warning(f'{method=} not in {methods=}.')

    if mode == "load":
        # check if calibration folder is present
        if not PATHS["calibration"].exists():
            PATHS["calibration"].mkdir()
            logger.warning(
                "No calibration data found. To obtain quantitative EGA data run .calibrate(mode='recalibrate')!"
            )
            return None, None, None, None
        os.chdir(PATHS["calibration"])

        # try to load saved calibration
        try:
            cali = pd.read_excel(PATHS["calibration"] / "cali.xlsx", sheet_name=None, index_col=[0,1], header=[0,1])
            #cali["data"] = pd.read_excel("cali.xlsx", sheet_name="data", index_col=[0, 1])
            # cali["data"].reset_index(inplace=True)
            # cali["data"].set_index(cali["data"].columns[:3].to_list(), inplace=True)
            cali = {name: data.xs(profile).pint.quantify() for name, data in cali.items()}

            gases = cali["linreg"].index

        except Exception as e:
            logger.warning(
                f"No calibration data found. To obtain quantitative IR data run .calibrate(mode='recalibrate')! {e}"
            )
            os.chdir(PATHS["home"])
            return None, None, None, None

    # new calibration
    elif mode == "recalibrate":
        if not len(profiles:=set([profile for profile in worklist.profiles]))==1:
            logger.error(f"The worklist contains multiple profiles: {profiles!r}")
            return
        start_time=time()
        # check if calibration folder is present
        if not PATHS["calibration"].exists():
            # make directory
            PATHS["calibration"].mkdir()
        os.chdir(PATHS["calibration"])
        output_path = Path(f"{start_time}_calibration_{profile}")
        output_path.mkdir()

        # setting up output DataFrames
        
        names = ["linreg","stats","x_mass","x_mol","y","data"]
        cali = {name: pd.DataFrame() for name in names}

        if not worklist:
            logger.error(
                "Supply Worklist-object with calibration measurements with 'worklist=' and rerun this command."
            )
            os.chdir(PATHS["home"])
            return None, None, None, None

        # calculating mass steps and integrating FTIR_data signals for all samples
        logger.info("Calibrating...")
        figdim = len(worklist), 2
        if plot:
            fig, axs = plt.subplots(*figdim,sharex="col", gridspec_kw={"hspace":.5, "wspace":.5}, figsize=FIGSIZE*figdim[::-1])

        data = []
        for i, sample_data in enumerate(worklist):
            # calculating mass steps and integrating FTIR_data signals
            # calculate width in terms of indices from supplied temperature width
            #width_idxs = width_T / sp.stats.mode(np.diff(sample_data.tga.sample_temp.values._data)).mode
            steps, _, step_starts_idx, step_ends_idx, peaks_idx = sample_data.mass_step(
                plot=plot,
                ax=axs[i, 0],
                min_rel_height=min_rel_height, #only allow positive peaks
                width_T=width_T,
            )
            integrals = integrate_peaks(
                sample_data.ega,
                step_starts_idx,
                step_ends_idx,
                peaks_idx,
                plot=plot,
                ax=axs[i, 1],
                corr_baseline=corr_baseline,
                gases=sample_data.info.gases,
            )

            integrals.insert(
                loc=0, column="mass loss", value=steps,
            )
            data.append(pd.concat({sample_data.name: integrals}, names=["samples", "step"]))
            
            axs[i,0].set_title(sample_data.alias)
            if i!=len(worklist)-1:
                axs[i, 0].set_xlabel("")
                axs[i, 1].set_xlabel("")
        logger.info("Finished integrating data.")
        cali["data"] = pd.concat(data)
        cali["data"].pint.dequantify().to_excel(output_path / f"integration_data.xlsx")
        if plot:
            fig.align_xlabels()
            plt.show()
            fig.savefig(output_path / f"integration.png")
        
        step_dtype = cali["data"]["mass loss"].dtype
        # assigning gases to mass steps
        for sample in cali["data"].index.levels[0]:
            release_steps = []
            for i, step in enumerate(
                cali["data"].loc[sample, "mass loss"]
            ):
                integrals = cali["data"].loc[sample].drop(
                    ["mass loss"], axis=1
                )
                
                # check gas that has its highest area in the particular step
                norm = integrals.divide(integrals.max(axis=0).values, axis=1).loc[i]
                gases = norm.loc[norm == 1].index.to_list()
             
                assignments = {gas: molecular_formulas.get(gas, gas) for gas in gases}
                release_steps.append(assignments.values())

                # check duplicates
                if (l:=len(set(assignments.values()))) > 1:
                    dupl_assignments = [trace if trace==mf else f"{trace} ({mf})" for trace, mf in assignments.items()]
                    logger.warning(f"{sample!r}: {l} EGA-traces ('{"', '".join(dupl_assignments)}') were assigned to mass step {i+1}. Treating them as independant.")
                
                # fill x_mass and y
                for gas in gases:
                    gas_dtype = cali["data"][gas].dtype
                    if gas not in cali["x_mass"].columns:
                        cali["x_mass"][gas] = pd.Series(dtype=step_dtype)
                        cali["y"][gas] = pd.Series(dtype=gas_dtype)
                    if sample not in cali["x_mass"].index:
                        cali["x_mass"] = pd.concat([cali["x_mass"], pd.DataFrame(index=[sample])])
                        cali["y"] = pd.concat([cali["y"], pd.DataFrame(index=[sample])])
                    
                    cali["x_mass"].loc[sample, gas] = step
                    cali["y"].loc[sample, gas] = integrals.loc[i, gas]

        # calculate molar calibration of mass values
        logger.info(f"Performing calibration based on {method=}")
        cali["x_mass"] = cali["x_mass"].pint.convert_object_dtype()
        cali["y"] = cali["y"].pint.convert_object_dtype()
       
        match method:
            case "max":
                calibrate_max(cali, molecular_formulas)
            case "iter":
                logger.error(f"{method=!r} is currently not available.")
                #calibrate_iter(cali, molecular_formulas)
            case "co_oxi":
                logger.error(f"{method=!r} is currently not available.")
                #calibrate_co_oxi(cali, molecular_formulas)
            case "co_oxi_iter":
                logger.error(f"{method=!r} is currently not available.")
                #calibrate_co_oxi_iter(cali, molecular_formulas)
            case "mlr":
                logger.error(f"{method=!r} is currently not available.")
                #calibrate_mlr(cali, molecular_formulas)

        cali["linreg"] = cali["linreg"].pint.convert_object_dtype()
        cali["stats"] = calibration_stats(cali["x_mol"], cali["y"], cali["linreg"])        
        logger.info("Calibration finished.")

        # saving of dfs
        # reset step index for easier loading
        cali["data"].reset_index(level="step", inplace=True)
        cali = {name:data.pint.convert_object_dtype() for name, data in cali.items()}
        try:
            path = PATHS["calibration"]/ "cali.xlsx"
            path_exists=path.exists()
            with pd.ExcelWriter(path, mode="a" if path_exists else "w", engine="openpyxl", if_sheet_exists="replace" if path_exists else None) as writer:
                for name, data in cali.items():
                    # transform before saving
                    data = data.pint.dequantify()
                    data = data.rename({"":"dimensionless"}, level=1, axis=1)
                    data_new = pd.concat({profile:data}, names=["profile"])

                    # get old contents of file and update rows with new data
                    if path_exists:
                        data_old = pd.read_excel(writer, sheet_name=name, header=[0,1], index_col=[0,1])
                        data_new = pd.concat([data_old, data_new]).drop_duplicates(keep="last")
                    data_new.to_excel(writer, sheet_name=name,merge_cells=MERGE_CELLS)

                logger.info(f"Calibration completed, data is stored under {path.as_posix()!r}")
        except PermissionError:
            logger.error(f"Could not write on {path.as_posix()!r}. Close file and rerun command.")

    # plotting
    if plot:
        gases = cali["x_mass"].columns
        plot_gases = [gas for gas in gases if gas in cali["linreg"].index]
        figdim = len(plot_gases),2
        fig, axs = plt.subplots(*figdim ,gridspec_kw = {"hspace":.5, "wspace":.25}, figsize=FIGSIZE*figdim[::-1],**fig_args)
        for i,gas in enumerate(plot_gases):
            x = cali["x_mol"][gas]
            y = cali["y"][gas]
            plot_calibration_single(x, y, cali["linreg"].loc[gas,:], axs[i, 0])
            plot_residuals_single(x, y, cali["linreg"].loc[gas,:], axs[i, 1])

            # set labels, titles only at borders
            if i == 0:
                axs[i, 0].set_title("Regression")
                axs[i, 1].set_title("Residuals")
            if i == len(plot_gases)-1:
                axs[i, 0].set_xlabel(UNITS["molar_amount"])
                axs[i, 1].set_xlabel(f"$\\hat{{y}}_i$ {SEP} {UNITS['int_ega']}")
        fig.savefig(output_path / f"regression.png")

        x = cali["x_mol"]
        y = cali["y"].pint.convert_object_dtype()
        fig, axs = plot_calibration_combined(x, y, cali["linreg"], plot_gases)
        fig.savefig(output_path / f"regression_combined.png")

    # return home
    os.chdir(PATHS["home"])
    return cali["linreg"], cali["stats"], cali["x_mol"], cali["y"]

def get_regression_df():
    pass

def calibrate_max(cali: dict, molecular_formulas:dict) -> dict:
    cali["x_mol"] = pd.DataFrame(index = cali["x_mass"].index, columns = cali["x_mass"].columns)
    invalid_mfs = []
    for gas in cali["x_mass"].columns:
        molecular_formula = gas if gas not in molecular_formulas else molecular_formulas[gas]
        if validate_mf(molecular_formula):
            molar_mass = ureg.Quantity(Formula(molecular_formula).mass, "g/mol")
        else:
            invalid_mfs.append(molecular_formula)
            cali["x_mol"].drop(gas, axis=1, inplace=True)
            continue
        molar_amount = (cali["x_mass"][gas] / molar_mass).pint.to_base_units()
        cali["x_mol"][gas] = molar_amount
        x = cali["x_mol"][gas].values._data#.dropna(axis=0)#.astype(float)
        y = cali["y"][gas].values._data#.dropna(axis=0)#.astype(float)
        regr_res = list(sp.stats.linregress(x, y))
        regression = pd.DataFrame(
            [[molecular_formula, molar_mass]+regr_res], index=[gas], columns=REGR_COLS
        )

        if gas not in cali["linreg"].index:
            cali["linreg"] = pd.concat([cali["linreg"], regression], verify_integrity=True)
        else:
            cali["linreg"].loc[[gas]] = regression

    if invalid_mfs:
        invalid_mfs_msg = ", ".join([f"{mf!r}: 'Formula {i}'" for i, mf in enumerate(invalid_mfs)])
        logger.error(
                    f"Could not determine molar mass for {invalid_mfs!r}. Skipping. Supply valid molecular formula with 'molecular_formulas = {{{invalid_mfs_msg}}}'")

def calibrate_iter(cali: dict, molecular_formulas:dict):
    for gas in cali["x_mass"].columns:
        # slope, intercept, r_value, p_value, std_err=sp.stats.linregress(x_cali[gas].dropna(axis=0).astype(float),y_cali[gas].dropna(axis=0).astype(float))
        x = cali["x_mass"][gas].dropna(axis=0).astype(float)
        y = cali["y"][gas].dropna(axis=0).astype(float)
        regression = pd.DataFrame(
            [sp.stats.linregress(x, y)], index=[gas], columns=REGR_COLS
        )

        if gas not in cali["linreg"].index:
            cali["linreg"] = pd.concat([cali["linreg"], regression], verify_integrity=True)
        else:
            cali["linreg"].loc[[gas]] = regression

    n_iter = 10
    X_cali = pd.DataFrame()  # attention: X_cali is not x_cali
    temp_linreg = pd.DataFrame(index=release_steps, columns=REGR_COLS)
    for i in range(n_iter):
        for step in cali["data"].index.levels[1]:
            gas = release_steps[step]
            X_cali[gas] = cali["data"].loc[
                (slice(None), step), "mass loss",
            ].droplevel(1)
            for other in set(release_steps) - set([gas]):
                corr = (
                    cali["data"].loc[(slice(None), step), other].droplevel(1)
                    - cali["linreg"]["intercept"][other]
                ) / cali["linreg"]["slope"][other]
                X_cali[gas] = np.subtract(X_cali[gas], corr * (corr > 0))

            x = X_cali[gas].values.astype(float)
            y = cali["data"].loc[(slice(None), step), gas].values.astype(float)
            temp_linreg.update(
                pd.DataFrame(
                    np.array([sp.stats.linregress(x, y)]),
                    columns=REGR_COLS,
                    index=[gas],
                )
            )

        cali["linreg"].update(temp_linreg)

    for step in cali["data"].index.levels[1]:
        gas = release_steps[step]
        for other in set(release_steps) - set([gas]):
            cali["x_mass"][gas] = cali["x_mass"][gas].subtract(
                (
                    cali["data"].loc[(slice(None), step), other].droplevel(1)
                    - cali["linreg"]["intercept"][other]
                )
                / cali["linreg"]["slope"][other]
            )
def calibrate_co_oxi(cali: dict, molecular_formulas:dict):
    co_corr = (
        (cali["data"].loc[(slice(None), 1), "CO2"] - cali["linreg"].loc["CO2", "intercept"])
        / cali["linreg"].loc["CO2", "slope"]
    ).values
    for i in range(
        len(co_corr)
    ):  # apply correction step wise to check each value
        j = cali["x_mass"].index[i]
        if (
            co_corr[i] > 0.0
        ):  # if the correction would be negative, this is due to the intercept and should not be applied!
            cali["x_mass"].loc[j, "CO"] = cali["x_mass"].loc[j, "CO"] - co_corr[i]

    x = cali["x_mass"]["CO"].dropna(axis=0).astype(float)
    y = cali["y"]["CO"].dropna(axis=0).astype(float)
    regression = pd.DataFrame(
        [sp.stats.linregress(x, y)], index=["CO"], columns=REGR_COLS
    )
    cali["linreg"].loc[["CO"]] = regression

def calibrate_co_oxi_iter(cali: dict, molecular_formulas:dict):
    n_iter = 10
    X_cali = cali["x_mass"].copy(deep=True)  # attention: X_cali is not x_cali
    for i in range(n_iter):
        co_corr = (
            (cali["data"].loc[(slice(None), 1), "CO2"] - cali["linreg"].loc["CO2", "intercept"])
            / cali["linreg"].loc["CO2", "slope"]
        ).values
        for j in range(
            len(co_corr)
        ):  # apply correction step wise to check each value
            k = X_cali.index[j]
            if (
                co_corr[j] > 0.0
            ):  # if the correction would be negative, this is due to the intercept and should not be applied!
                X_cali.loc[k, "CO"] = cali["x_mass"].loc[k, "CO"] - co_corr[j]
        x = X_cali["CO"].dropna(axis=0).astype(float)
        y = cali["y"]["CO"].dropna(axis=0).astype(float)
        regression = pd.DataFrame(
            [sp.stats.linregress(x, y)], index=["CO"], columns=REGR_COLS
        )
        cali["linreg"].loc[["CO"]] = regression

        co2_corr = (
            (cali["data"].loc[(slice(None), 2), "CO"] - cali["linreg"].loc["CO", "intercept"])
            / cali["linreg"].loc["CO", "slope"]
        ).values
        for j in range(
            len(co2_corr)
        ):  # apply correction step wise to check each value
            k = X_cali.index[j]
            if (
                co2_corr[j] > 0.0
            ):  # if the correction would be negative, this is due to the intercept and should not be applied!
                X_cali.loc[k, "CO2"] = cali["x_mass"].loc[k, "CO2"] - co2_corr[j]
        x = X_cali["CO2"].dropna(axis=0).astype(float)
        y = cali["y"]["CO2"].dropna(axis=0).astype(float)
        regression = pd.DataFrame(
            [sp.stats.linregress(x, y)], index=["CO2"], columns=REGR_COLS
        )
        cali["linreg"].loc[["CO2"]] = regression

    cali["x_mass"] = X_cali.copy(deep=True)


def calibrate_mlr(cali: dict, molecular_formulas:dict):
    # Y_cali and X_cali is different to methods above, here it is more like the whole data dataframe
    # For better perfomance, a kind of x_cali and y_cali should be reconstructed, not to show confusing calibrations plots,
    # because the result ist quite good. But not needed for now...
    cali["data"] = cali["data"].clip(lower=0)  # this is not necessary? (but reasonable)
    X_cali = cali["data"].loc[
        (slice(None), slice(None)), "mass loss",
    ]  # attention: Y_cali is not y_cali
    Y_cali = cali["data"].drop(["mass loss"], axis=1).dropna(
        axis=1
    )  # attention: X_cali is not x_cali, dropna() to exclude gases not available for every measurement
    mlr = (
        linear_model.LinearRegression()
    )  # RANSACRegressor()#fit_intercept=0.0)
    mlr.fit(Y_cali, X_cali)

    # fill linreg manually
    for i, gas in enumerate(Y_cali.columns):
        cali["linreg"].loc[[gas]] = pd.DataFrame(
            [
                [
                    1 / mlr.coef_[i] * Formula(gas).mass,
                    0,
                    np.nan,
                    np.nan,
                    np.nan,
                ]
            ],
            index=[gas],
            columns=REGR_COLS,
        )

    # create new x_cali and y_cali
    gases = cali["linreg"].index
    cali["x_mass"] = Y_cali.copy()
    for step in cali["data"].index.levels[1]:
        step_sums = [0] * len(cali["data"].index.levels[0])
        for gas in gases:
            step_sums += Y_cali.loc[(slice(None), step), gas].values / (
                cali["linreg"].loc[gas, "slope"] * cali["linreg"]["molmass"][gas]
            )
        for gas in gases:
            step_gas = Y_cali.loc[(slice(None), step), gas].values / (
                cali["linreg"].loc[gas, "slope"] * cali["linreg"]["molmass"][gas]
            )
            cali["x_mass"].loc[(slice(None), step), gas] = (
                X_cali[(slice(None), step)].values * step_gas / step_sums
            ) / cali["linreg"]["molmass"][gas]
    cali["y"] = Y_cali

    # from new x_cali and y_cali a linear regression is calculated
    # comment this lines out to get the result of mlr only
    for gas in gases:
        x = cali["x_mass"][gas].dropna(axis=0).astype(float)
        y = cali["y"][gas].dropna(axis=0).astype(float)
        regression = pd.DataFrame(
            [sp.stats.linregress(x, y)], index=[gas], columns=REGR_COLS
        )

        if gas not in cali["linreg"].index:
            cali["linreg"] = pd.concat([cali["linreg"], regression], verify_integrity=True)
        else:
            cali["linreg"].loc[[gas]] = regression