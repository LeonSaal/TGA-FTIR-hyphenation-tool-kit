import matplotlib.pyplot as plt
from ..config import BOUNDS, UNITS
import pandas as pd

# labels and units for plots


def plot_robustness(data, save=False, ylim=[0, None], **kwargs):
    params = ["center_0", "tolerance_center", "hwhm_max", "height_0", "hwhm_0"]
    results = dict()
    results["data"] = pd.DataFrame()
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
    labels = [
        "$center$",
        "$tolerance\,center$",
        "$HWHM_{max}$",
        "$height_0$",
        "$HWHM_0$",
    ]
    units = ["°C", "°C", "°C", "$height_{max}$", "°C"]

    total_cond = data.index.get_level_values("group") != "total"
    stat_cond = ~data.index.get_level_values("gas").isin(["sum", "mean"])
    for sample, df in data[total_cond & stat_cond].groupby("sample"):
        for param, label, unit in zip(params, labels, units):
            fig = plt.figure()

            means = df[df.index.get_level_values("run") == "mean"]
            dev = df[df.index.get_level_values("run") == "err_dev"]
            x = [
                f"{group.replace('_',' ').capitalize()}, {gas}"
                for group, gas in zip(
                    means.index.get_level_values("group"),
                    means.index.get_level_values("gas"),
                )
            ]
            for run, i in zip(["+", "init", "-"], [-1, 0, 1]):
                if run == "init":
                    col = run
                    var = 0
                else:
                    prefix = f"{param}{run}"
                    (col,) = [col for col in df.columns if col.startswith(prefix)]
                    var = float(col.removeprefix(prefix))

                y = means[col]
                yerr = dev[col]
                if param == "hwhm_0":  # exception for correct hwhm_0 labeling
                    glabel = f"{default[param] + i * default[param] * var} {unit}"
                else:
                    glabel = f"{default[param] + i * var if param != 'center_0' else f'{default[param] + i * var:+}'} {unit}"
                plt.errorbar(
                    x, y, yerr, label=glabel, marker="x", capsize=10, ls="none"
                )
            plt.legend()
            plt.title(f"{sample}: {label}")
            plt.ylim(ylim)
            plt.ylabel(f"${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}$")
            plt.legend()
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            plt.show()

            if save:
                sample_name = "".join(
                    [x if (x.isalnum() or x in "._- ") else "_" for x in sample]
                )  # to catch invalide sample names
                fig.savefig(f'{sample_name}_{param}.png')
