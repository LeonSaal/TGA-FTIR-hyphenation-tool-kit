import copy
import logging
import os
import time as tm

import numpy as np
import pandas as pd

from ..config import BOUNDS, MERGE_CELLS, PATHS, DEFAULTS
from ..input_output.general import time
from .fitting import get_presets

logger = logging.getLogger(__name__)


def robustness(
    worklist,
    reference,
    T_max=None,
    save=True,
    var_T=10,
    var_rel=0.3,
    **kwargs,
):
    "perform robustness test on multiple Sample objects"

    # load default presets
    presets_rob = get_presets(reference)

    if presets_rob is None:
        return

    # setting up output DataFrames and lists of parameters to test
    params = ["center_0", "tolerance_center", "hwhm_max", "height_0", "hwhm_0"]
    results = dict()
    results["data"] = pd.DataFrame()
    variance = dict(zip(params, [var_T, var_T, var_T, var_rel, var_rel]))
    default = dict(
        zip(
            params,
            [
                0,
                BOUNDS.getfloat("tol_center"),
                BOUNDS.getfloat("hwhm_max"),
                BOUNDS.getfloat("height_0"),
                BOUNDS.getfloat("hwhm_0"),
            ],
        )
    )

    # Calculating initial results and timing initial fit for estimation of remaining time
    # silence loggers of sample and worklist to reduce amount of output
    logging.getLogger("TGA_FTIR_tools.classes.sample").setLevel(logging.WARNING)
    logging.getLogger("TGA_FTIR_tools.classes.worklist").setLevel(logging.WARNING)
    
    start = tm.time()
    res = worklist.fit(
        reference,
        save=False,
        T_max=T_max,
        presets=presets_rob,
        plot=False,
        mod_samples = False,
        **kwargs,
    )["n_rel"]
    t_fit = tm.time() - start
    logger.info(f"Initial results took {t_fit:.1f} s to compute.")
    results["init"] = res.rename("init")
    del res

    gases = [key for key in presets_rob]
    n_fits = 2 * (4 + sum([len(presets_rob[gas]) for gas in presets_rob]))
    t_total = (t_fit * n_fits / 60) * 1.1
    logger.info(f"Approximate remaining time: {t_total:.1f} min.")
    logger.info(f"Varying fitting parameters...")

    # cycling through parameters to vary
    fin_fits = 1
    for key in params:
        for i in [-1, 1]:
            logger.info(
                f"Fit {fin_fits} of {n_fits} (~{t_total*(n_fits-fin_fits)/n_fits:.1f} min remaining): {key} {variance[key] * i:+}"
            )
            temp_presets = copy.deepcopy(presets_rob)
            run = f"{key}{i*variance[key]:+}"
            if key == "center_0":
                results[run] = pd.Series(
                    index=results["init"].index, name=run, dtype=results["init"].dtype
                )
                for k, gas in enumerate(gases, start=1):
                    # vary center of each group
                    logger.info(f"  Gas {k} of {len(gases)}: {gas}")
                    for j, group in enumerate(
                        (groups := temp_presets[gas].index), start=1
                    ):
                        fin_fits += 1
                        cols = temp_presets[gas].drop("link", axis=1).columns
                        temp_presets = copy.deepcopy(presets_rob)
                        logger.info(f"    Group {j} of {len(groups)}: {group.capitalize()}")
                        temp_presets[gas].loc[group, cols] = presets_rob[gas].loc[
                            group, cols
                        ] + i * np.array(
                            [
                                variance[key],
                                0,
                                0,
                                variance[key],
                                0,
                                0,
                                variance[key],
                                0,
                                0,
                            ]
                        )
                        #logging.disable(logging.CRITICAL)

                        res = worklist.fit(
                            plot=False,
                            reference=reference,
                            save=False,
                            presets=temp_presets,
                            mod_samples=False,
                            **kwargs,
                        )
                        #logging.disable(logging.NOTSET)

                        results[run].update(
                            res.loc[
                                (
                                    slice(None),
                                    slice(None),
                                    slice(None),
                                    slice(None),
                                    slice(None),
                                    group,
                                    gas,
                                ),
                                "n_rel",
                            ]
                        )
                        del res
            else:
                for gas in gases:
                    cols = temp_presets[gas].drop("link", axis=1).columns
                    if key == "hwhm_max":
                        temp_presets[gas].loc[:, cols] += i * np.array(
                            [0, 0, 0, 0, 0, 0, 0, variance[key], 0]
                        )
                    elif key == "tolerance_center":
                        temp_presets[gas].loc[:, cols] += i * np.array(
                            [0, 0, 0, -variance[key], 0, 0, variance[key], 0, 0]
                        )
                    elif key == "height_0":
                        temp_presets[gas][key] = temp_presets[gas][
                            key[: key.rfind("_")] + "_max"
                        ] * (default[key] + i * variance[key])
                    elif key == "hwhm_0":
                        temp_presets[gas][key] += (
                            i * temp_presets[gas][key] * variance[key]
                        )

                        # temp_presets[gas][key[:key.rfind('_')]+'_max'] * (default[key]+i*variance[key])
                        # temp_presets[gas][key] = temp_presets[gas][key[:key.rfind('_')]+'_max'] * (default[key]+i*variance[key])
                fin_fits += 1
                res = worklist.fit(
                    reference=reference,
                    save=False,
                    T_max=T_max,
                    plot=False,
                    presets=temp_presets,
                    mod_samples=False,
                    **kwargs,
                )["n_rel"]
                results[run] = res.rename(run)

    logger.info(f"Fittings finished! Took {(tm.time() - start) / 60:.2f} min.")
    logger.info("Calculating statistic values.")

    data = pd.concat([df for df in results.values()], axis=1)
    data.index.rename(names=["reference","name","sample","alias", "run", "group", "gas"], inplace=True)

    # make further statistical data
    summary = data.pint.convert_object_dtype().astype(np.float64).groupby(["sample", "group", "gas"]).agg(["mean", "std", "min", "max"]).reset_index()
    # stat_names = ["err_dev", "rel_err_dev", "rel_stddev", "stddev", "mean"]
    # stat_cond = ~data.index.get_level_values(1).isin(stat_names)
    # columns = ["mean", "stddev", "meanstddev", "min", "max"]
    # stats = []
    # for (sample, group, gas), df in data[stat_cond].groupby(["alias", "group", "gas"]):
    #     print(sample, group, gas, df)
    #     yall = df.dropna(axis=1).values.flatten()
    #     ymean = data.loc[sample, "mean", group, gas].dropna().values.flatten()
    #     values = [
    #         [
    #             np.mean(ymean),
    #             np.std(yall),
    #             np.std(ymean),
    #             np.min(yall),
    #             np.max(yall),
    #         ]
    #     ]
    #     index = pd.MultiIndex.from_tuples(
    #         [(sample, group, gas)], names=["sample", "group", "gas"]
    #     )
    #     stat_df = pd.DataFrame(values, columns=columns, index=index)
    #     stats.append(stat_df)

    # summary = pd.concat(stats)

    # save results to excel file
    # if save:
    #     logger.info(f"Plots and results are saved.'{path=}'.")
    #     with pd.ExcelWriter("robustness_in_n_rel.xlsx") as writer:
    #         summary.to_excel(writer, sheet_name="summary", merge_cells=MERGE_CELLS)
    #         data.to_excel(writer, sheet_name="data", merge_cells=MERGE_CELLS)

    # os.chdir(PATHS["home"])
    logger.info("Robustness test finished!")

    # reestablish loglevel
    logging.getLogger("TGA_FTIR_tools.classes.sample").setLevel(DEFAULTS.get("logging_level", logging.INFO))
    logging.getLogger("TGA_FTIR_tools.classes.worklist").setLevel(DEFAULTS.get("logging_level", logging.INFO))
    
    data.sort_index(level="group", inplace=True)
    return data.drop("total", level="group").pint.convert_object_dtype(), summary.pint.convert_object_dtype()
