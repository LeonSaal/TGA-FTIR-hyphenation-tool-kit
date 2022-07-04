import matplotlib.pyplot as plt
from ..config import PARAMS, SEP, UNITS, PATHS
import copy
import matplotlib.ticker as ticker
import os
from .utils import ylim_auto, get_label
from ..input_output.general import time
import logging


logger = logging.getLogger(__name__)


def plots(
    TG_IR_objs,
    plot,
    x_axis="sample_temp",
    y_axis="orig",
    ylim="auto",
    xlim=[None, None],
    gas=None,
    save=False,
    legend=True,
    reference_mass="reference_mass",
    linewidth=1,
):
    "overlay plots from different objects"

    # setting up axis-labels and catching possible input errors
    if plot == "TG":
        ylabel = "sample_mass"
    else:
        ylabel = plot.lower()
    if plot == "IR":
        if gas == None:
            logger.warn("Supply 'gas = '")
            return
        else:
            gas = gas.upper()
            if gas not in TG_IR_objs[0].ir.columns:
                logger.warn(f"{gas} was not found in IR data.")
                return

        # just to see if supplied gas is calibrated or not
        calibrated = set()
        for TG_IR in TG_IR_objs:
            try:
                calibrated.update(set(TG_IR.linreg.index))
            except AttributeError:
                if y_axis == "rel":
                    logger.warn(f"{gas} is not calibrated for {TG_IR.info['name']}")

        if y_axis == "rel":
            if calibrated == set():
                logger.warn(
                    f"{gas} is not calibrated. Change y_axis to 'orig' and rerun command."
                )
                return
    fig, ax = plt.subplots()
    ax.set_xlabel(f"{PARAMS[x_axis.lower()]} {SEP} ${UNITS[x_axis.lower()]}$")
    if plot != "IR":
        if y_axis == "orig":
            ax.set_ylabel(f"{PARAMS[ylabel]} {SEP} ${UNITS[ylabel]}$")
        elif y_axis == "rel":
            if plot == "DTG":
                ax.set_ylabel(f"{PARAMS[ylabel]} {SEP} $\\%\\,min^{{-1}}$")
            else:
                ax.set_ylabel(f"{PARAMS[ylabel]} {SEP} $\\%$")
    elif plot == "IR":
        if y_axis == "orig":
            ax.set_ylabel(f"{get_label(gas)} {SEP} ${UNITS[ylabel]}$")
        elif y_axis == "rel":
            ax.set_ylabel(
                f"{get_label(gas)} {SEP} ${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}$"
            )

    # actual plotting
    for obj in TG_IR_objs:
        if reference_mass == "reference_mass":
            ref_mass = obj.info[obj.info[reference_mass]]
        else:
            ref_mass = obj.info[reference_mass]

        if plot == "TG":
            x = copy.deepcopy(obj.tga[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = obj.tga["sample_mass"]
            elif y_axis == "rel":
                y = 100 * obj.tga["sample_mass"] / ref_mass
            if (
                ylim == "auto"
            ):  # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(
                x,
                y,
                linewidth=linewidth,
                label=f"{obj.info['alias']}, {obj.info['reference_mass']}: {obj.info[obj.info['reference_mass']]:.2f}$\\,{UNITS['sample_mass']}$",
            )
        if plot == "DTG":
            x = copy.deepcopy(obj.tga[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = obj.tga["dtg"] * 60
            elif y_axis == "rel":
                y = obj.tga["dtg"] * 60 / ref_mass * 100
            if (
                ylim == "auto"
            ):  # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(
                x,
                y,
                linewidth=linewidth,
                label=f"{obj.info['alias']}, {obj.info['reference_mass']}: {obj.info[obj.info['reference_mass']]:.2f}$\\,{UNITS['sample_mass']}$",
            )
        if plot == "heat_flow":
            x = copy.deepcopy(obj.tga[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = obj.tga["heat_flow"]
            elif y_axis == "rel":
                y = obj.tga["heat_flow"] / ref_mass
            if (
                ylim == "auto"
            ):  # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(
                x,
                y,
                linewidth=linewidth,
                label=f"{obj.info['alias']}, {obj.info['reference_mass']}: {obj.info[obj.info['reference_mass']]:.2f}$\\,{UNITS['sample_mass']}$",
            )
        if plot == "IR":
            x = copy.deepcopy(obj.ir[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = obj.ir[gas]
                if (
                    ylim == "auto"
                ):  # only select relevant range of x data, to auto-scale the y axis
                    x, y, ylim_temp = ylim_auto(x, y, xlim)
                ax.plot(
                    x,
                    y,
                    linewidth=linewidth,
                    label=f"{obj.info['alias']}, {obj.info['reference_mass']}: {obj.info[obj.info['reference_mass']]:.2f}$\\,{UNITS['sample_mass']}$",
                )
            elif y_axis == "rel":
                y = obj.ir[gas] / obj.linreg["slope"][gas] / ref_mass
                if (
                    ylim == "auto"
                ):  # only select relevant range of x data, to auto-scale the y axis
                    x, y, ylim_temp = ylim_auto(x, y, xlim)
                ax.plot(
                    x,
                    y,
                    linewidth=linewidth,
                    label=f"{obj.info['alias']}, {obj.info['reference_mass']}: {obj.info[obj.info['reference_mass']]:.2f}$\\,{UNITS['sample_mass']}$",
                )

    if ylim == "auto":  # reset ylim to [None,None]
        ylim = ylim_temp
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if legend:
        ax.legend()

    ax.xaxis.set_minor_locator(
        ticker.AutoMinorLocator()
    )  # switch on minor ticks on each axis
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    plt.show()

    if save:
        path_plots = PATHS["plots"]
        if os.path.exists(path_plots) == False:
            os.makedirs(path_plots)
        if gas == None:
            gas = ""
        fig.savefig(
            os.path.join(path_plots, "_".join([time(), plot, gas, y_axis])) + ".png",
        )

