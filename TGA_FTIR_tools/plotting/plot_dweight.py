import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from ..config import PARAMS, SEP, UNITS, PATHS
from .plotting import get_label
import os
import numpy as np


def plot_dweight(
    sample, save=False, xlim=[None, None], ylim=[None, None], title=True, how_dry="H2O"
):
    times, names = sample.info.step_time, sample.info.mass_steps
    weights = sample.tga[sample.tga.index.isin(times)]["sample_mass"].values
    mass_loss = abs(np.diff(weights))

    fig, ax = plt.subplots()
    x = sample.tga["sample_temp"]
    y = sample.tga["sample_mass"]
    ax.plot(x, y, label="TGA")

    if how_dry == "H2O":
        try:
            ref = sample.ir.filter(items=["sample_temp", "H2O"])
            ylabel = get_label(how_dry)
        except:
            how_dry = "sample_mass"

    if how_dry == "sample_mass":
        ref = sample.tga.filter(items=["sample_temp", "sample_mass"])
        ylabel = "DTG"

    for i in range(len(times) - 1):
        ax.annotate(
            text="",
            xy=(x[times[i + 1]], y[times[i]]),
            xytext=(x[times[i + 1]], y[times[i + 1]]),
            arrowprops=dict(arrowstyle="<->"),
        )
        ax.text(
            x[times[i + 1]] + 20,
            (y[times[i]] + y[times[i + 1]]) / 2,
            f"$ML$ {get_label(names[i])}: {mass_loss[i]:.2f} mg ({mass_loss[i] / sample.info[sample.info['reference_mass']] * 100:.1f} %)",
        )
    ax.hlines(weights[:-1], x[times[:-1]], x[times[1:]], color="black")
    ax.set_ylabel(f"{PARAMS['sample_mass']} {SEP} {UNITS['sample_mass']}")
    ax.set_xlabel(f"{PARAMS['sample_temp']} {SEP} {UNITS['sample_temp']}")
    ax.set_ylim(ylim)
    if type(how_dry) == str:
        ax2 = plt.twinx()
        ax2.plot(ref["sample_temp"], ref[how_dry], linestyle="dashed", label=ylabel)
        ax2.set_ylabel(ylabel)
        h2, l2 = ax2.get_legend_handles_labels()
    ax.set_xlim(xlim)

    ax.xaxis.set_minor_locator(
        ticker.AutoMinorLocator()
    )  # switch on minor ticks on each axis
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    if title == True:
        ax.set(title="Dry mass and mass steps determination")

    plt.legend()
    plt.show()
    if save:
        path_plots = PATHS["plots"]
        if not path_plots.exists():
            os.makedirs(path_plots)
        fig.savefig(path_plots/f"{sample.info.name}_mass_steps.png")
