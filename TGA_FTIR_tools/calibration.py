import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from chempy import Substance
from sklearn import linear_model

from .classes import Sample
from .config import COUPLING, MERGE_CELLS, PATHS, SAVGOL, SEP, UNITS
from .plotting import get_label

logger = logging.getLogger(__name__)


def integrate_peaks(
    FTIR_data, step_start_idx, step_end_idx, peaks_idx, corr_baseline=None, plot=False, ax=None, gases=None
):
    "integrating IR signal in between given bounds"
    gases = list(gases)
    integrals = pd.DataFrame(index=range(len(step_start_idx)), columns=gases)

    # plotting
    if plot:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        x = FTIR_data["time"]
        y = FTIR_data[gases[0]]
        step_starts = x[step_start_idx]
        step_ends = x[step_end_idx]
        peaks = x[peaks_idx]
        # setup figure and plot first gas
        graph = []
        #fig.subplots_adjust(right=0.8)
        graph.append(ax)
        graph[0].set_xlabel(f'{get_label("time")} {SEP} { UNITS["time"]}')
        graph[0].set_ylabel(f'{get_label(gases[0])} {SEP} {UNITS["ir"]}')
        graph[0].yaxis.label.set_color(colors[0])
        graph[0].plot(x, y)
        graph[0].set_ylim(0 - (max(y) / 20), max(y))
        for step_start, step_end in zip(step_starts, step_ends):
            graph[0].axvspan(step_start, step_end, alpha=.5)
        graph[0].vlines(peaks, 0, max(y), linestyle="dotted")

        # append secondary, third... y-axis on right side
        for i, gas in enumerate(gases[1:]):
            y = FTIR_data[gas]

            graph.append(graph[0].twinx())
            graph[i + 1].spines["right"].set_position(("axes", 1 + i * 0.1))
            graph[i + 1].plot(x, y, color=colors[i + 1])

            graph[i + 1].set_ylabel(f"{get_label(gas)} {SEP} {UNITS['ir']}")
            graph[i + 1].yaxis.label.set_color(colors[i + 1])
            graph[i + 1].set_ylim(0 - (max(y) / 20), max(y))

    # integration
    for gas in gases:
        for i in range(len(step_end_idx)):
            subset = FTIR_data[gas].iloc[step_start_idx[i]:step_end_idx[i]]
            

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

            if plot:
                x = (
                    FTIR_data["time"].iloc[step_start_idx[i]:step_end_idx[i]]
                )
                graph[gases.index(gas)].plot(
                    x, baseline, color=colors[gases.index(gas)], linestyle="dashed"
                )

    return integrals


def eval_lin(x, slope, intercept):
    "evaluating linear equation with slope and intercept at ponts x"
    return slope * x + intercept


def calibration_stats(x_cali, y_cali, linreg, alpha=0.95, beta=None, m=1, k=3):
    "calculating calibration stats according to DIN 32645"
    gases = linreg.index

    n = len(x_cali)
    if n == 2:
        return pd.DataFrame()

    f = n - 2
    if not beta:
        beta = alpha

    stats = pd.DataFrame()
    for gas in gases:
        b = linreg["slope"][gas]
        a = linreg["intercept"][gas]

        s_yx = np.sqrt(np.sum(np.power(b * x_cali[gas] + a - y_cali[gas], 2)) / (n - 2))
        s_x0 = s_yx / b
        x_ = np.mean(x_cali[gas])
        Q_x = np.sum(np.power(x_cali[gas] - x_, 2))

        x_NG = (
            s_x0 * sp.stats.t.ppf(alpha, f) * np.sqrt(1 / m + 1 / n + (x_ * x_) / Q_x)
        )

        # x_EG=x_NG+s_x0*sp.stats.t.ppf(beta,f)*np.sqrt(1/m+1/n+(x_*x_)/Q_x)

        x_BG = k * x_NG

        stats = pd.concat(
            [
                stats,
                pd.DataFrame(
                    [[s_yx, s_x0, x_NG, x_BG]],
                    index=[gas],
                    columns=["s_yx", "s_x0", "x_LOD", "x_LOQ"],
                ),
            ]
        )
    return stats


def calibrate(worklist=None, molecular_formulas = {},plot=False, mode="load", method="max", profile=None, corr_baseline=None, rel_height=.95, **fig_args):
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
            return None, None
        os.chdir(PATHS["calibration"])

        # try to load saved calibration
        try:
            cali = pd.read_excel(PATHS["calibration"] / "cali.xlsx", sheet_name=None, index_col=[0,1])
            cali["data"] = pd.read_excel("cali.xlsx", sheet_name="data", index_col=[0, 1])
            cali["data"].reset_index(inplace=True)
            cali["data"].set_index(cali["data"].columns[:3].to_list(), inplace=True)
            cali = {name: data.xs(profile, level="profile") for name, data in cali.items()}

            gases = cali["linreg"].index

        except Exception as e:
            logger.warning(
                f"No calibration data found. To obtain quantitative IR data run .calibrate(mode='recalibrate')! {e}"
            )
            os.chdir(PATHS["home"])
            return None, None

    # new calibration
    elif mode == "recalibrate":
        # check if calibration folder is present
        if not PATHS["calibration"].exists():
            # make directory
            PATHS["calibration"].mkdir()
        os.chdir(PATHS["calibration"])

        # setting up output DataFrames
        
        names = ["linreg","stats",f"x in {UNITS['molar_amount']}",f"y in {UNITS['int_ir']}","data"]
        cali = {name: pd.DataFrame() for name in names}

        if not worklist:
            logger.error(
                "Supply Worklist-object with calibration measurements with 'worklist=' and rerun this command."
            )
            os.chdir(PATHS["home"])
            return None, None

        # calculating mass steps and integrating FTIR_data signals for all samples
        logger.info("Calibrating...")

        fig, axs = plt.subplots(len(worklist), 2,sharex="col", gridspec_kw={"hspace":.5, "wspace":.5}, **fig_args)

        for i, sample_data in enumerate(worklist):

            # calculating mass steps and integrating FTIR_data signals
            steps, _, step_starts_idx, step_ends_idx, peaks_idx = sample_data.mass_step(
                plot=plot,
                ax=axs[i, 0],
                height=0.00025,
                width=50,
                prominence=0.00025,
                rel_height=rel_height,
            )
            integrals = integrate_peaks(
                sample_data.ir,
                step_starts_idx,
                step_ends_idx,
                peaks_idx,
                plot=plot,
                ax=axs[i, 1],
                corr_baseline=corr_baseline,
                gases=sample_data.info.gases,
            )

            integrals.insert(
                loc=0, column=f"mass loss in {UNITS['sample_mass']}", value=steps,
            )
            cali["data"] = pd.concat(
                [cali["data"], pd.concat({sample_data.name: integrals}, names=["samples", "step"])]
            )
            axs[i,0].set_title(sample_data.alias)
            if i!=len(worklist)-1:
                axs[i, 0].set_xlabel("")
                axs[i, 1].set_xlabel("")
        logger.info("Finished integrating data.")

        if plot:
            fig.align_xlabels()
            plt.show()

        # assigning gases to mass steps
        for sample in cali["data"].index.levels[0]:
            release_steps = []
            for i, step in enumerate(
                cali["data"].loc[sample, f"mass loss in {UNITS['sample_mass']}"]
            ):
                integrals = cali["data"].loc[sample].drop(
                    [f"mass loss in {UNITS['sample_mass']}"], axis=1
                )
                
                norm = integrals.divide(integrals.max(axis=0).values, axis=1).loc[i]
                gases = norm.loc[norm == 1].index.to_list()
             
                assignments = {gas: molecular_formulas.get(gas, gas) for gas in gases}
                release_steps.append(assignments.values())

                # check duplicates
                if (l:=len(set(assignments.values()))) > 1:
                    dupl_assignments = [trace if trace==mf else f"{trace} ({mf})" for trace, mf in assignments.items()]
                    logger.warning(f"{sample!r}: {l} EGA-traces ('{"', '".join(dupl_assignments)}') were assigned to mass step {i+1}. Treating them as independant.")
                
                for gas in gases:
                    if gas not in cali[f"x in {UNITS['molar_amount']}"].columns:
                        cali[f"x in {UNITS['molar_amount']}"][gas] = np.nan
                        cali[f"y in {UNITS['int_ir']}"][gas] = np.nan
                    if sample not in cali[f"x in {UNITS['molar_amount']}"].index:
                        cali[f"x in {UNITS['molar_amount']}"] = pd.concat([cali[f"x in {UNITS['molar_amount']}"], pd.DataFrame(index=[sample])])
                        cali[f"y in {UNITS['int_ir']}"] = pd.concat([cali[f"y in {UNITS['int_ir']}"], pd.DataFrame(index=[sample])])
                    cali[f"x in {UNITS['molar_amount']}"].loc[sample, gas] = step
                    cali[f"y in {UNITS['int_ir']}"].loc[sample, gas] = integrals.loc[i, gas]

        cols = ["molmass", "slope", "intercept", "r_value", "p_value", "std_error"]
        gases = cali[f"x in {UNITS['molar_amount']}"].columns

        logger.info(f"Performing calibration based on {method=}")
        if method == "iter":
            for gas in gases:
                # slope, intercept, r_value, p_value, std_err=sp.stats.linregress(x_cali[gas].dropna(axis=0).astype(float),y_cali[gas].dropna(axis=0).astype(float))
                x = cali[f"x in {UNITS['molar_amount']}"][gas].dropna(axis=0).astype(float)
                y = cali[f"y in {UNITS['int_ir']}"][gas].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=[gas], columns=cols
                )

                if gas not in cali["linreg"].index:
                    cali["linreg"] = pd.concat([cali["linreg"], regression], verify_integrity=True)
                else:
                    cali["linreg"].loc[[gas]] = regression

            n_iter = 10
            X_cali = pd.DataFrame()  # attention: X_cali is not x_cali
            temp_linreg = pd.DataFrame(index=release_steps, columns=cols)
            for i in range(n_iter):
                for step in cali["data"].index.levels[1]:
                    gas = release_steps[step]
                    X_cali[gas] = cali["data"].loc[
                        (slice(None), step), f"mass loss in {UNITS['sample_mass']}",
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
                            columns=cols,
                            index=[gas],
                        )
                    )

                cali["linreg"].update(temp_linreg)

            for step in cali["data"].index.levels[1]:
                gas = release_steps[step]
                for other in set(release_steps) - set([gas]):
                    cali[f"x in {UNITS['molar_amount']}"][gas] = cali[f"x in {UNITS['molar_amount']}"][gas].subtract(
                        (
                            cali["data"].loc[(slice(None), step), other].droplevel(1)
                            - cali["linreg"]["intercept"][other]
                        )
                        / cali["linreg"]["slope"][other]
                    )

        # regression (method 'max')
        invalid_mfs = []
        for gas in gases:
            molecular_formula = gas if gas not in molecular_formulas else molecular_formulas[gas]
            try:
                molar_mass = Substance.from_formula(molecular_formula).mass
            except ValueError:
                invalid_mfs.append(molecular_formula)
                continue
            cali[f"x in {UNITS['molar_amount']}"].update(cali[f"x in {UNITS['molar_amount']}"][gas] / molar_mass)
            x = cali[f"x in {UNITS['molar_amount']}"][gas].dropna(axis=0).astype(float)
            y = cali[f"y in {UNITS['int_ir']}"][gas].dropna(axis=0).astype(float)
            regression = pd.DataFrame(
                [[molar_mass]+list(sp.stats.linregress(x, y))], index=[gas], columns=cols
            )

            if gas not in cali["linreg"].index:
                cali["linreg"] = pd.concat([cali["linreg"], regression], verify_integrity=True)
            else:
                cali["linreg"].loc[[gas]] = regression

        if invalid_mfs:
            invalid_mfs_msg = ", ".join([f"{mf!r}: 'Formula {i}'" for i, mf in enumerate(invalid_mfs)])
            logger.error(
                        f"Could not determine molar mass for {invalid_mfs!r}. Skipping. Supply valid molecular formula with 'molecular_formulas = {{{invalid_mfs_msg}}}'")
                    
        if method == "co_oxi":
            co_corr = (
                (cali["data"].loc[(slice(None), 1), "CO2"] - cali["linreg"].loc["CO2", "intercept"])
                / cali["linreg"].loc["CO2", "slope"]
            ).values
            for i in range(
                len(co_corr)
            ):  # apply correction step wise to check each value
                j = cali[f"x in {UNITS['molar_amount']}"].index[i]
                if (
                    co_corr[i] > 0.0
                ):  # if the correction would be negative, this is due to the intercept and should not be applied!
                    cali[f"x in {UNITS['molar_amount']}"].loc[j, "CO"] = cali[f"x in {UNITS['molar_amount']}"].loc[j, "CO"] - co_corr[i]

            x = cali[f"x in {UNITS['molar_amount']}"]["CO"].dropna(axis=0).astype(float)
            y = cali[f"y in {UNITS['int_ir']}"]["CO"].dropna(axis=0).astype(float)
            regression = pd.DataFrame(
                [sp.stats.linregress(x, y)], index=["CO"], columns=cols
            )
            cali["linreg"].loc[["CO"]] = regression

        if method == "co_oxi_iter":
            n_iter = 10
            X_cali = cali[f"x in {UNITS['molar_amount']}"].copy(deep=True)  # attention: X_cali is not x_cali
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
                        X_cali.loc[k, "CO"] = cali[f"x in {UNITS['molar_amount']}"].loc[k, "CO"] - co_corr[j]
                x = X_cali["CO"].dropna(axis=0).astype(float)
                y = cali[f"y in {UNITS['int_ir']}"]["CO"].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=["CO"], columns=cols
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
                        X_cali.loc[k, "CO2"] = cali[f"x in {UNITS['molar_amount']}"].loc[k, "CO2"] - co2_corr[j]
                x = X_cali["CO2"].dropna(axis=0).astype(float)
                y = cali[f"y in {UNITS['int_ir']}"]["CO2"].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=["CO2"], columns=cols
                )
                cali["linreg"].loc[["CO2"]] = regression

            cali[f"x in {UNITS['molar_amount']}"] = X_cali.copy(deep=True)

        if method == "mlr":
            # Y_cali and X_cali is different to methods above, here it is more like the whole data dataframe
            # For better perfomance, a kind of x_cali and y_cali should be reconstructed, not to show confusing calibrations plots,
            # because the result ist quite good. But not needed for now...
            cali["data"] = cali["data"].clip(lower=0)  # this is not necessary? (but reasonable)
            X_cali = cali["data"].loc[
                (slice(None), slice(None)), f"mass loss in {UNITS['sample_mass']}",
            ]  # attention: Y_cali is not y_cali
            Y_cali = cali["data"].drop([f"mass loss in {UNITS['sample_mass']}"], axis=1).dropna(
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
                            1 / mlr.coef_[i] * Substance(gas).mass,
                            0,
                            np.nan,
                            np.nan,
                            np.nan,
                        ]
                    ],
                    index=[gas],
                    columns=cols,
                )

            # create new x_cali and y_cali
            gases = cali["linreg"].index
            cali[f"x in {UNITS['molar_amount']}"] = Y_cali.copy()
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
                    cali[f"x in {UNITS['molar_amount']}"].loc[(slice(None), step), gas] = (
                        X_cali[(slice(None), step)].values * step_gas / step_sums
                    ) / cali["linreg"]["molmass"][gas]
            cali[f"y in {UNITS['int_ir']}"] = Y_cali

            # from new x_cali and y_cali a linear regression is calculated
            # comment this lines out to get the result of mlr only
            for gas in gases:
                x = cali[f"x in {UNITS['molar_amount']}"][gas].dropna(axis=0).astype(float)
                y = cali[f"y in {UNITS['int_ir']}"][gas].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=[gas], columns=cols
                )

                if gas not in cali["linreg"].index:
                    cali["linreg"] = pd.concat([cali["linreg"], regression], verify_integrity=True)
                else:
                    cali["linreg"].loc[[gas]] = regression

        cali["stats"] = calibration_stats(cali[f"x in {UNITS['molar_amount']}"], cali[f"y in {UNITS['int_ir']}"], cali["linreg"])
        logger.info("Calibration finished.")
        # saving of
        try:
            path = PATHS["calibration"]/ "cali.xlsx"
            with pd.ExcelWriter(path) as writer:
                for name, data in cali.items():
                    pd.concat({profile:data}, names=["profile"]).to_excel(writer, sheet_name=name,merge_cells=MERGE_CELLS)

                logger.info(f"Calibration completed, data is stored under {path.as_posix()!r}")
        except PermissionError:
            logger.error(f"Could not write on {path.as_posix()!r}. Close file and rerun command.")

    # plotting
    if plot:
        
        plot_gases = [gas for gas in gases if gas in cali["linreg"].index]
        fig, axs = plt.subplots(len(plot_gases),2 ,gridspec_kw = {"hspace":.5, "wspace":.25}, **fig_args)
        for i,gas in enumerate(plot_gases):
            fig = plt.figure()
            x = cali[f"x in {UNITS['molar_amount']}"][gas]
            y = cali[f"y in {UNITS['int_ir']}"][gas]

            axs[i, 0].scatter(x, y, label=f"data (N = {len(x)})")

            x_bounds = np.array((min(x), max(x)))
            axs[i, 0].plot(
                x_bounds,
                x_bounds * cali["linreg"]["slope"][gas] + cali["linreg"]["intercept"][gas],
                label="regression",
                ls="dashed",
            )
            axs[i, 0].text(
                max(x),
                min(y),
                rf'y={cali["linreg"]["slope"][gas]:.1e}x{cali["linreg"]["intercept"][gas]:+.1e}, $R^2$={cali["linreg"]["r_value"][gas] ** 2:.3}',
                horizontalalignment="right",
            )
            label_mf = get_label(gas) if gas not in molecular_formulas else get_label(molecular_formulas[gas])
            label = f"{label_mf}" if gas not in molecular_formulas else f"{gas}\n({label_mf})"
            axs[i, 0].set_ylabel(label)
            axs[i, 0].set_xlim(0, max(x) + abs(min(x)))

            Y_cali = x.mul(cali["linreg"]["slope"][gas]).add(cali["linreg"]["intercept"][gas])
            axs[i, 1].scatter(Y_cali, y - Y_cali, label=f"data (N = {len(x)})")
            axs[i, 1].hlines(0, min(Y_cali), max(Y_cali))
            
            axs[i, 1].set_ylabel(f"$y_i-\\hat{{y}}_i$ {SEP} {UNITS['int_ir']}")
            if i == 0:
                axs[i, 0].set_title("Regression")
                axs[i, 1].set_title("Residuals")
            if i == len(plot_gases)-1:
                axs[i, 0].set_xlabel(UNITS["molar_amount"])
                axs[i, 1].set_xlabel(f"$\\hat{{y}}_i$ {SEP} {UNITS['int_ir']}")
            

            #axs[i, 1].set_title(f"Residual plot: {get_label(gas)}")


        for gas in gases:
            if gas not in cali["linreg"].index:
                continue
            x = cali[f"x in {UNITS['molar_amount']}"][gas]
            y = cali[f"y in {UNITS['int_ir']}"][gas]

            plt.scatter(x, y, label=f"data {get_label(gas)} (N = {len(x)})")

            x = np.array((min(x), max(x)))
            plt.plot(
                x,
                x * cali["linreg"]["slope"][gas] + cali["linreg"]["intercept"][gas],
                label=f"regression {get_label(gas)}",
                ls="dashed",
            )
            plt.xlabel(UNITS["molar_amount"])
            plt.ylabel(UNITS["int_ir"])
            plt.xlim(0, max(x) + abs(min(x)))
            plt.legend(loc=0)
        plt.show()

    os.chdir(PATHS["home"])
    return cali["linreg"], cali["stats"]


# try loading calibration data
linreg, stats = calibrate()
