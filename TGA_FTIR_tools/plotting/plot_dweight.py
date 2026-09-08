import os

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from ..config import PATHS, SEP, UNITS
from .plotting import get_label
import pint

ureg = pint.get_application_registry()


def plot_dweight(
    sample, save=False, xlim=[None, None], ylim=[None, None], title=True, how_dry="H2O"
):  
    step_data = sample.step_data()
    weights, names, times = step_data.sample_mass, step_data.step.to_list(), step_data.index
    mass_loss = abs(np.diff(weights))

    fig, ax = plt.subplots()
    x = sample.tga["sample_temp"]
    y = sample.tga["sample_mass"]
    ax.plot(x, y, label="TGA")

    if how_dry == "H2O":
        try:
            ref = sample.ega.filter(items=["sample_temp", "H2O"])
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
            x[times[i + 1]] + ureg.Quantity(20, "delta_degreeC"),
            (y[times[i]] + y[times[i + 1]]) / 2,
            f"$ML$ {get_label(names[i])}: {mass_loss[i]:.2f} ({(mass_loss[i] / sample.reference_mass * 100).magnitude:.1f}%)",
        )
    ax.hlines(weights[:-1], x[times[:-1]], x[times[1:]], color="black")
    ax.set_ylabel(f"{get_label('sample_mass')} {SEP} {UNITS.get('sample_mass', '?')}")
    ax.set_xlabel(f"{get_label('sample_temp')} {SEP} {UNITS.get('sample_temp', '?')}")
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
    if save:
        path_plots = PATHS["plots"]
        if not path_plots.exists():
            os.makedirs(path_plots)
        fig.savefig(path_plots/f"{sample.info.name}_mass_steps.png")
