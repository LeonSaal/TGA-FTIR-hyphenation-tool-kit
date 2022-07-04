import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pint
from ..config import PATHS
from ..input_output import time
import os

import logging


logger = logging.getLogger(__name__)


def bar_plot_results(
    worklist,
    show_groups=[],
    group_by="sample",
    y_unit="µmol per g",
    res="robustness",
    save=False,
    title=True,
    w_group=0.9,
    w_bar=0.9,
    exclude_merged=True,
    exclude_total=True,
):
    "Returns a bar plot of the robustness tested results with units conversion and two group_by options."
    dev = "stddev"
    if res == "robustness":
        if res not in worklist.results:
            logger.warn(f"No results to plot. Run .robustness().")
        data = worklist.results["robustness"][1]

    elif res == "fit":
        if "fit" not in worklist.results:
            logger.warn(f"No results to plot. Run .fit().")

        data = worklist.results["fit"]
        use_cols = ["mean", "err_dev", "stddev"]
        use_index = data.index.get_level_values("run").isin(use_cols)
        index = ["sample", "group", "gas"]
        data = (
            data["mmol_per_mg"][use_index]
            .reset_index()
            .pivot(index=index, columns=["run"])["mmol_per_mg"]
        )
    else:
        options = ["fit", "robutsness"]
        logger.warn(f"{res=} not in {options=}.")
        return

    if show_groups:
        data = data[
            data.index.get_level_values("group").isin(show_groups)
        ]  # select and order df index/groups
    if exclude_merged:
        data = data[~data.index.get_level_values("gas").isin(["mean", "sum"])]

    if exclude_total:
        data = data[~(data.index.get_level_values("group") == "total")]

    factor = pint.Quantity("mmol per mg").to(y_unit).magnitude
    groups = data.groupby(["group", "gas"])
    samples = data.groupby("sample")
    if group_by == "group":
        labels = [
            f'{group.replace("_"," ").capitalize()}, {gas}'
            for group, gas in groups.groups.keys()
        ]
        grouped = samples
    elif group_by == "sample":
        labels = [group for group in samples.groups.keys()]
        grouped = groups
    else:
        options = ["sample", "group"]
        logger.warn(f"{group_by=} not in {options=}.")
        return

    N = len(labels)
    if N > 4:
        rot = 45
        ha = "center"
        va = "top"
    else:
        rot = 0
        ha = "center"
        va = "top"

    fig, ax = plt.subplots()
    x_ticks = np.arange(N)
    ax.set_xticks(x_ticks, labels=labels, rotation=rot, ha=ha, va=va)

    for i, (key, values) in enumerate(grouped):
        if group_by == "sample":
            group, gas = key
            subset = values.loc[(slice(None), group, gas)]
            y = subset["mean"] * factor
            yerr = subset[dev] * factor
            label = f'{group.replace("_"," ").capitalize()}, {gas}'

        if group_by == "group":
            subset = values.loc[(key, slice(None), slice(None))]
            y = subset["mean"] * factor
            yerr = subset[dev] * factor
            label = key

        x = np.arange(len(y)) + w_group * (1 / (len(grouped)) * (1 / 2 + i) - 1 / 2)
        ax.bar(
            x,
            y,
            yerr=yerr,
            width=(w_group / len(grouped)) * w_bar,
            label=label,
            capsize=5,
        )
    ax.legend(bbox_to_anchor=(1, 1), loc="upper left")
    ax.set_ylabel(f"SOG in {y_unit}")
    if title:
        ax.set_title(f"summary plot with errors from {res}")
    plt.show()
    if save:
        path_plots_eval = os.path.join(PATHS["plots"], "Evaluation")
        if os.path.exists(path_plots_eval) == False:
            os.mkdir(path_plots_eval)
        fig.savefig(
            os.path.join(path_plots_eval, f"{time()}_{worklist.name}_by_{group_by}.png"),
        )


# def bar_plot_results(
#     sample,
#     show_groups=[],
#     group_by="samples",
#     y_unit="mymol_per_g",
#     res = 'robustness',
#     x_gap=1,
#     save=False,
#     title=True,
# ):
#     "Returns a bar plot of the robustness tested results with units conversion and two group_by options."
#     if res == 'robustness':
#         data = sample.results["robustness"][1]
#     elif res== 'fit':
#         data = sample.results["fit"]

#     samples = data.index.get_level_values("sample").unique()

#     if len(show_groups) > 0:
#         data = data[
#             data.index.get_level_values("group").isin(show_groups)
#         ]  # select and order df index/groups

#     fig, ax = plt.subplots()
#     if group_by == "samples":
#         x = samples  # samples as x-ticks
#         groups = data.index.get_level_values("group")  # groups for each sample
#     elif group_by == "groups":
#         x = data.index.get_level_values("group")  # groups as x-ticks
#         groups = samples  # samples for each group

#     y = data[data.index.get_level_values('gas')=='mean']['mmol_per_mg']
#     yerr= data[data.index.get_level_values('gas')=='err_dev']['mmol_per_mg']

#     ax.bar(x, y, y_err= yerr)

#     x_num = np.arange(len(x))  # define x-ticks nummerically
#     width = 1 / (
#         len(groups) + x_gap
#     )  # calculate width of a bar; +1 is the gap to bars of the next sample

#     # find maximum y value in the data, for adjusting of labels on bars
#     y_max = 0.0
#     for group in groups:
#         if group_by == "samples":
#             max_group = data.loc[group, (slice(None), "mmol_per_mg")].max()
#         elif group_by == "groups":
#             max_group = data.loc[x, (group, "mmol_per_mg")].max()

#         if max_group > y_max:
#             y_max = max_group

#     # rotate x-axis labels and labels on bars dependent on length of data to plot
#     y_lab_rotation = "horizontal"
#     y_lab_space = y_max / 100.0
#     if (2 * len(x) + len(groups)) > 14:
#         y_lab_rotation = "vertical"
#         y_lab_space = y_max / 50.0

#     if y_unit == "mymol_per_g":
#         y_lab_space = y_lab_space * 1000000

#     if len(x) > 5:
#         if len(groups) > 4:
#             plt.setp(
#                 ax.get_xticklabels(), rotation=30, horizontalalignment="center"
#             )  # rotate x-axis labels, positioning 'center'
#         else:
#             plt.setp(
#                 ax.get_xticklabels(), rotation=30, horizontalalignment="right"
#             )  # rotate x-axis labels, positioning 'right'

#     # loop through the bars (defines as groups)
#     for i, group in enumerate(groups):
#         # define y values according to
#         if group_by == "samples":
#             y = df.loc[group, (slice(None), "mmol_per_mg")]
#             y_err = df.loc[group, (slice(None), "stddev")]
#             try:
#                 y_lab = df.loc[group, (slice(None), "label")]
#             except:
#                 pass
#         elif group_by == "groups":
#             y = df.loc[x, (group, "mmol_per_mg")]
#             y_err = df.loc[x, (group, "stddev")]
#             try:
#                 y_lab = df.loc[x, (group, "label")]
#             except:
#                 pass

#         x_i = x_num - width * (len(groups) - 1) / 2 + i * width

#         if y_unit == "mymol_per_g":
#             y = y * 1000000
#             y_err = y_err * 1000000
#         else:
#             pass

#         ax.bar(
#             x_i, y, yerr=y_err, capsize=5, width=width * 0.9, label=group, zorder=10
#         )  # color = c,

#         # add labels to bars
#         for j in range(len(x_i)):
#             x_j = x_i[j]
#             y_j = y.iloc[j] + y_lab_space
#             text = y_lab.iloc[j]
#             ax.text(
#                 x_j, y_j, text, horizontalalignment="center", rotation=y_lab_rotation
#             )

#     ax.legend()
#     ax.yaxis.grid(True)

#     if title == True:
#         ax.set(title="summary plot with errors from robustness testing")

#     if y_unit == "mymol_per_g":
#         ax.set_ylabel("surface oxygen groups in $\mu mol\,g^{-1}$")
#     else:
#         ax.set_ylabel("surface oxygen groups in $mmol\,mg^{-1}$")
#     ax.set_xticks(np.arange(len(x)))
#     ax.set_xticklabels(x)

#     fig.tight_layout()
#     plt.show()

#     if save:
#         path_plots_eval = os.path.join(PATHS["plots"], "Evaluation")
#         if os.path.exists(path_plots_eval) == False:
#             os.makedirs(path_plots_eval)
#         sample_names = "".join(
#             [x if (x.isalnum() or x in "._- ") else "" for x in str(samples)]
#         )  # to catch invalide sample names
#         # check path length and if necessary shorten file name by list of samples
#         path_save = os.path.join(
#             path_plots_eval, f"{time()}_{sample_names}_by_{group_by}.png"
#         )
#         if len(path_save) > 260:
#             sample_names = sample_names[: (len(sample_names) - (len(path_save) - 259))]
#         fig.savefig(
#             os.path.join(path_plots_eval, f"{time()}_{sample_names}_by_{group_by}.png"),
#         )

#     return
