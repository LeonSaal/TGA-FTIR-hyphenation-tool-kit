import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from chempy import Substance
from sklearn import linear_model

from .classes import Sample
from .config import COUPLING, MERGE_CELLS, PARAMS, PATHS, SAVGOL, SEP, UNITS
from .plotting import get_label

logger = logging.getLogger(__name__)


def integrate_peaks(
    FTIR_data, step_start, step_end, corr_baseline=None, plot=False, gases=None
):
    "integrating IR signal in between given bounds"
    gases = list(gases)
    # adsjusting step starts according to coupling delay
    step_start = step_start + int(60 * COUPLING.getfloat("coupling_delay"))
    step_end = step_end + int(60 * COUPLING.getfloat("coupling_delay"))
    integrals = pd.DataFrame(index=range(len(step_start)), columns=gases)

    # plotting
    if plot:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        x = FTIR_data["time"] / 60
        # setup figure and plot first gas
        graph = []
        fig, ax = plt.subplots()
        fig.subplots_adjust(right=0.8)
        graph.append(ax)

        graph[0].set_xlabel(f'{PARAMS["time"]} {SEP} { UNITS["time"]}')
        graph[0].set_ylabel(f'{get_label(gases[0])} {SEP} {UNITS["ir"]}')
        graph[0].yaxis.label.set_color(colors[0])
        graph[0].plot(x, FTIR_data[gases[0]])
        graph[0].set_ylim(0 - (max(FTIR_data[gases[0]]) / 20), max(FTIR_data[gases[0]]))

        # append secondary, third... y-axis on right side
        for i, gas in enumerate(gases[1:]):
            y = FTIR_data[gas]

            graph.append(graph[0].twinx())
            graph[i + 1].spines["right"].set_position(("axes", 1 + i * 0.1))
            graph[i + 1].plot(x, y, color=colors[i + 1])
            graph[i + 1].vlines(step_start / 60, 0, max(y), linestyle="dashed")
            graph[i + 1].vlines(step_end / 60, 0, max(y), linestyle="dashed")
            graph[i + 1].set_ylabel(f"{get_label(gas)} {SEP} {UNITS['ir']}")
            graph[i + 1].yaxis.label.set_color(colors[i + 1])
            graph[i + 1].set_ylim(0 - (max(y) / 20), max(y))

    # integration
    for gas in gases:
        for i in range(len(step_end)):
            subset = FTIR_data[gas][
                (FTIR_data["time"] >= step_start[i])
                & (FTIR_data["time"] <= step_end[i])
            ]
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

            integral = sp.integrate.simps(subset - baseline)
            integrals.loc[i, gas] = integral

            if plot:
                x = (
                    FTIR_data["time"][
                        (FTIR_data["time"] >= step_start[i])
                        & (FTIR_data["time"] <= step_end[i])
                    ]
                    / 60
                )
                graph[gases.index(gas)].plot(
                    x, baseline, color=colors[gases.index(gas)], linestyle="dashed"
                )

    if plot:
        plt.show()

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


def calibrate(plot=False, mode="load", method="max"):
    methods = ['max', "iter", "co_oxi", "co_oxi_iter", 'mlr']
    if method not in methods:
        logger.warn(f'{method=} not in {methods=}.')

    if mode == "load":
        # check if calibration folder is present
        if not PATHS["calibration"].exists():
            logger.warn(
                "No calibration data found. To obtain quantitative IR data supply an 'Calibration' folder in the home directory containing cali.xlsx or run TGA_FTIR_tools.calibrate(mode='recalibrate')!"
            )
            return None, None
        os.chdir(PATHS["calibration"])

        # try to load saved calibration
        try:
            cali = pd.read_excel("cali.xlsx", sheet_name=None, index_col=0)
            linreg = cali["linreg"]
            x_cali = cali[f"x in {UNITS['molar_amount']}"]
            y_cali = cali[f"y in {UNITS['int_ir']}"]
            stats = cali["stats"]
            data = pd.read_excel("cali.xlsx", sheet_name="data", index_col=[0, 1])
            gases = linreg.index
        except:
            logger.warn(
                "No calibration data found. To obtain quantitative IR data supply an 'Calibration' folder in the home directory containing cali.xlsx or run TGA_FTIR_tools.calibrate(mode='recalibrate')!"
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
        x_cali = pd.DataFrame()
        y_cali = pd.DataFrame()
        linreg = pd.DataFrame()
        stats = pd.DataFrame()
        data = pd.DataFrame()

        # imporitng sample list for calibration
        fname_sl = "Sample_list.txt"
        if os.path.exists(fname_sl):
            samples = pd.read_csv(fname_sl, delimiter="\t")
        else:
            with open("Sample_list.txt", "w") as file:
                file.write("Samples\tBaseline")
            logger.info(
                "'Sample_list.txt' was created in the 'Calibration' folder, please fill in calibration measurements and rerun this command."
            )
            os.chdir(PATHS["home"])
            return None, None

        # calculating mass steps and integrating FTIR_data signals for all samples
        logger.info("Calibrating...")
        for sample, baseline in zip(samples["Samples"], samples["Baseline"]):
            sample_data = Sample(sample)
            sample_data.corr(baseline)

            # calculating mass steps and integrating FTIR_data signals
            [steps, _, stepstart, stepend] = sample_data.mass_step(
                plot=plot,
                height=0.00025,
                width=50,
                prominence=0.00025,
                rel_height=0.98,
            )
            integrals = integrate_peaks(
                sample_data.ir,
                stepstart,
                stepend,
                plot=plot,
                corr_baseline=None,
                gases=sample_data.info.gases,
            )

            integrals.insert(
                loc=0, column=f"mass loss in {UNITS['sample_mass']}", value=steps,
            )
            data = pd.concat(
                [data, pd.concat({sample: integrals}, names=["samples", "step"])]
            )
        logger.info("Finished reading data.")

        # assigning gases to mass steps
        for sample in data.index.levels[0]:
            release_steps = []
            for i, step in enumerate(
                data.loc[sample, f"mass loss in {UNITS['sample_mass']}"]
            ):
                integrals = data.loc[sample].drop(
                    [f"mass loss in {UNITS['sample_mass']}"], axis=1
                )
                norm = integrals.divide(integrals.max(axis=0).values, axis=1).loc[i]
                gas = norm.loc[norm == 1].index.values
                release_steps.append(gas)
                if gas not in x_cali.columns.values:
                    x_cali[gas] = np.nan
                    y_cali[gas] = np.nan
                if sample not in x_cali.index:
                    x_cali = pd.concat([x_cali, pd.DataFrame(index=[sample])])
                    y_cali = pd.concat([y_cali, pd.DataFrame(index=[sample])])
                x_cali.loc[sample, gas] = step
                y_cali.loc[sample, gas] = integrals.loc[i, gas]

        cols = ["slope", "intercept", "r_value", "p_value", "std_error"]
        gases = x_cali.columns

        logger.info(f"Performing calibration based on {method=}")
        if method == "iter":
            for gas in gases:
                # slope, intercept, r_value, p_value, std_err=sp.stats.linregress(x_cali[gas].dropna(axis=0).astype(float),y_cali[gas].dropna(axis=0).astype(float))
                x = x_cali[gas].dropna(axis=0).astype(float)
                y = y_cali[gas].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=[gas], columns=cols
                )

                if gas not in linreg.index:
                    linreg = pd.concat([linreg, regression], verify_integrity=True)
                else:
                    linreg.loc[[gas]] = regression

            n_iter = 10
            X_cali = pd.DataFrame()  # attention: X_cali is not x_cali
            temp_linreg = pd.DataFrame(index=release_steps, columns=cols)
            for i in range(n_iter):
                for step in data.index.levels[1]:
                    gas = release_steps[step]
                    X_cali[gas] = data.loc[
                        (slice(None), step), f"mass loss in {UNITS['sample_mass']}",
                    ].droplevel(1)
                    for other in set(release_steps) - set([gas]):
                        corr = (
                            data.loc[(slice(None), step), other].droplevel(1)
                            - linreg["intercept"][other]
                        ) / linreg["slope"][other]
                        X_cali[gas] = np.subtract(X_cali[gas], corr * (corr > 0))

                    x = X_cali[gas].values.astype(float)
                    y = data.loc[(slice(None), step), gas].values.astype(float)
                    temp_linreg.update(
                        pd.DataFrame(
                            np.array([sp.stats.linregress(x, y)]),
                            columns=cols,
                            index=[gas],
                        )
                    )

                linreg.update(temp_linreg)

            for step in data.index.levels[1]:
                gas = release_steps[step]
                for other in set(release_steps) - set([gas]):
                    x_cali[gas] = x_cali[gas].subtract(
                        (
                            data.loc[(slice(None), step), other].droplevel(1)
                            - linreg["intercept"][other]
                        )
                        / linreg["slope"][other]
                    )

        # regression (method 'max')
        for gas in gases:
            x_cali.update(x_cali[gas] / (Substance.from_formula(gas).mass))
            x = x_cali[gas].dropna(axis=0).astype(float)
            y = y_cali[gas].dropna(axis=0).astype(float)
            regression = pd.DataFrame(
                [sp.stats.linregress(x, y)], index=[gas], columns=cols
            )

            if gas not in linreg.index:
                linreg = pd.concat([linreg, regression], verify_integrity=True)
            else:
                linreg.loc[[gas]] = regression

        if method == "co_oxi":
            co_corr = (
                (data.loc[(slice(None), 1), "CO2"] - linreg.loc["CO2", "intercept"])
                / linreg.loc["CO2", "slope"]
            ).values
            for i in range(
                len(co_corr)
            ):  # apply correction step wise to check each value
                j = x_cali.index[i]
                if (
                    co_corr[i] > 0.0
                ):  # if the correction would be negative, this is due to the intercept and should not be applied!
                    x_cali.loc[j, "CO"] = x_cali.loc[j, "CO"] - co_corr[i]

            x = x_cali["CO"].dropna(axis=0).astype(float)
            y = y_cali["CO"].dropna(axis=0).astype(float)
            regression = pd.DataFrame(
                [sp.stats.linregress(x, y)], index=["CO"], columns=cols
            )
            linreg.loc[["CO"]] = regression

        if method == "co_oxi_iter":
            n_iter = 10
            X_cali = x_cali.copy(deep=True)  # attention: X_cali is not x_cali
            for i in range(n_iter):
                co_corr = (
                    (data.loc[(slice(None), 1), "CO2"] - linreg.loc["CO2", "intercept"])
                    / linreg.loc["CO2", "slope"]
                ).values
                for j in range(
                    len(co_corr)
                ):  # apply correction step wise to check each value
                    k = X_cali.index[j]
                    if (
                        co_corr[j] > 0.0
                    ):  # if the correction would be negative, this is due to the intercept and should not be applied!
                        X_cali.loc[k, "CO"] = x_cali.loc[k, "CO"] - co_corr[j]
                x = X_cali["CO"].dropna(axis=0).astype(float)
                y = y_cali["CO"].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=["CO"], columns=cols
                )
                linreg.loc[["CO"]] = regression

                co2_corr = (
                    (data.loc[(slice(None), 2), "CO"] - linreg.loc["CO", "intercept"])
                    / linreg.loc["CO", "slope"]
                ).values
                for j in range(
                    len(co2_corr)
                ):  # apply correction step wise to check each value
                    k = X_cali.index[j]
                    if (
                        co2_corr[j] > 0.0
                    ):  # if the correction would be negative, this is due to the intercept and should not be applied!
                        X_cali.loc[k, "CO2"] = x_cali.loc[k, "CO2"] - co2_corr[j]
                x = X_cali["CO2"].dropna(axis=0).astype(float)
                y = y_cali["CO2"].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=["CO2"], columns=cols
                )
                linreg.loc[["CO2"]] = regression

            x_cali = X_cali.copy(deep=True)

        if method == "mlr":
            # Y_cali and X_cali is different to methods above, here it is more like the whole data dataframe
            # For better perfomance, a kind of x_cali and y_cali should be reconstructed, not to show confusing calibrations plots,
            # because the result ist quite good. But not needed for now...
            data = data.clip(lower=0)  # this is not necessary? (but reasonable)
            X_cali = data.loc[
                (slice(None), slice(None)), f"mass loss in {UNITS['sample_mass']}",
            ]  # attention: Y_cali is not y_cali
            Y_cali = data.drop([f"mass loss in {UNITS['sample_mass']}"], axis=1).dropna(
                axis=1
            )  # attention: X_cali is not x_cali, dropna() to exclude gases not available for every measurement
            mlr = (
                linear_model.LinearRegression()
            )  # RANSACRegressor()#fit_intercept=0.0)
            mlr.fit(Y_cali, X_cali)

            # fill linreg manually
            for i, gas in enumerate(Y_cali.columns):
                linreg.loc[[gas]] = pd.DataFrame(
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
            gases = linreg.index
            x_cali = Y_cali.copy()
            for step in data.index.levels[1]:
                step_sums = [0] * len(data.index.levels[0])
                for gas in gases:
                    step_sums += Y_cali.loc[(slice(None), step), gas].values / (
                        linreg.loc[gas, "slope"] * Substance(gas).mass
                    )
                for gas in gases:
                    step_gas = Y_cali.loc[(slice(None), step), gas].values / (
                        linreg.loc[gas, "slope"] * Substance(gas).mass
                    )
                    x_cali.loc[(slice(None), step), gas] = (
                        X_cali[(slice(None), step)].values * step_gas / step_sums
                    ) / Substance(gas).mass
            y_cali = Y_cali

            # from new x_cali and y_cali a linear regression is calculated
            # comment this lines out to get the result of mlr only
            for gas in gases:
                x = x_cali[gas].dropna(axis=0).astype(float)
                y = y_cali[gas].dropna(axis=0).astype(float)
                regression = pd.DataFrame(
                    [sp.stats.linregress(x, y)], index=[gas], columns=cols
                )

                if gas not in linreg.index:
                    linreg = pd.concat([linreg, regression], verify_integrity=True)
                else:
                    linreg.loc[[gas]] = regression

        stats = calibration_stats(x_cali, y_cali, linreg)
        logger.info("Calibration finished.")
        # saving of
        try:
            path = PATHS["calibration"]/ "cali.xlsx"
            with pd.ExcelWriter(path) as writer:
                linreg.to_excel(writer, sheet_name="linreg")
                x_cali.to_excel(writer, sheet_name=f"x in {UNITS['molar_amount']}")
                y_cali.to_excel(writer, sheet_name=f"y in {UNITS['int_ir']}")
                stats.to_excel(writer, sheet_name="stats")
                data.to_excel(writer, sheet_name="data",merge_cells=MERGE_CELLS)

                logger.info(f"Calibration completed, data is stored under {path=}")
        except PermissionError:
            logger.error(f"Could not write on {path=}. Close file and rerun command.")

    # plotting
    if plot:
        for gas in gases:
            fig = plt.figure()
            x = x_cali[gas]
            y = y_cali[gas]

            plt.scatter(x, y, label=f"data (N = {len(x)})")

            x_bounds = np.array((min(x), max(x)))
            plt.plot(
                x_bounds,
                x_bounds * linreg["slope"][gas] + linreg["intercept"][gas],
                label="regression",
                ls="dashed",
            )
            plt.text(
                max(x),
                min(y),
                f'y = {linreg["slope"][gas]:.3f} $\cdot$ x {linreg["intercept"][gas]:+.3f}, $R^2$ = {linreg["r_value"][gas] ** 2:.5}',
                horizontalalignment="right",
            )
            plt.xlabel(UNITS["molar_amount"])
            plt.ylabel(UNITS["int_ir"])
            plt.title(get_label(gas))
            plt.xlim(0, max(x) + abs(min(x)))
            plt.legend(loc=0)
            plt.show()

            Y_cali = x.mul(linreg["slope"][gas]).add(linreg["intercept"][gas])
            plt.scatter(Y_cali, y - Y_cali, label=f"data (N = {len(x)})")
            plt.hlines(0, min(Y_cali), max(Y_cali))
            plt.xlabel(f"$\\hat{{y}}_i$ {SEP} {UNITS['int_ir']}")
            plt.ylabel(f"$y_i-\\hat{{y}}_i$ {SEP} {UNITS['int_ir']}")
            plt.title(f"Residual plot: {get_label(gas)}")
            plt.legend(loc=0)
            plt.show()

        plt.figure()
        for gas in gases:
            x = x_cali[gas]
            y = y_cali[gas]

            plt.scatter(x, y, label=f"data {get_label(gas)} (N = {len(x)})")

            x = np.array((min(x), max(x)))
            plt.plot(
                x,
                x * linreg["slope"][gas] + linreg["intercept"][gas],
                label=f"regression {get_label(gas)}",
                ls="dashed",
            )
            plt.xlabel(UNITS["molar_amount"])
            plt.ylabel(UNITS["int_ir"])
            plt.xlim(0, max(x) + abs(min(x)))
            plt.legend(loc=0)
        plt.show()

    os.chdir(PATHS["home"])
    return linreg, stats


# try loading calibration data
linreg, stats = calibrate()
