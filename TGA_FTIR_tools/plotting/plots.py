import copy
import logging
import os
from typing import Literal

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from ..config import PATHS, SEP, UNITS
from ..input_output.general import time
from .utils import get_label, make_title, ylim_auto

logger = logging.getLogger(__name__)


def plots(
    samples,
    plot: Literal["TG", "EGA", "DTG", "heat_flow"],
    ax=None,
    x_axis="sample_temp",
    y_axis="orig",
    ylim="auto",
    xlim=[None, None],
    gas=None,
    legend=True,
    reference_mass_name=None,
    linewidth=1,
    **kwargs,
):
    "overlay plots from different samples"
    options = ["TG", "EGA", "DTG", "heat_flow"]
    if plot not in options:
        logger.warning(f"{plot=} not in {options=}")
        return

    # setting up axis-labels and catching possible input errors
    if plot == "TG":
        ylabel = "sample_mass"
    else:
        ylabel = plot.lower()
    if plot == "EGA":
        if gas == None:
            logger.warning("Supply 'gas = '")
            return
        else:
            gas = gas.upper()
            if gas not in samples[0].ega.columns:
                logger.warning(f"{gas} was not found in IR data.")
                return

        # just to see if supplied gas is calibrated or not
        calibrated = set()
        for TG_IR in samples:
            try:
                calibrated.update(set(TG_IR.linreg.index))
            except AttributeError:
                if y_axis == "rel":
                    logger.warning(f"{gas} is not calibrated for {TG_IR.name}")

        if y_axis == "rel":
            if calibrated == set():
                logger.warning(
                    f"{gas} is not calibrated. Change y_axis to 'orig' and rerun command."
                )
                return

    ax.set_xlabel(f"{get_label(x_axis.lower())} {SEP} ${UNITS[x_axis.lower()]}$")
    if plot != "EGA":
        if y_axis == "orig":
            ax.set_ylabel(f"{get_label(ylabel)} {SEP} ${UNITS[ylabel]}$")
        elif y_axis == "rel":
            if plot == "DTG":
                ax.set_ylabel(f"{get_label(ylabel)} {SEP} $\\%\\,min^{{-1}}$")
            else:
                ax.set_ylabel(f"{get_label(ylabel)} {SEP} $\\%$")
    elif plot == "EGA":
        if y_axis == "orig":
            ax.set_ylabel(f"{get_label(gas)} {SEP} ${UNITS[ylabel]}$")
        elif y_axis == "rel":
            ax.set_ylabel(
                f"{get_label(gas)} {SEP} ${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}$"
            )

    # actual plotting
    for sample in samples:
        if reference_mass_name:
            step_data = sample.step_data()
            if reference_mass_name in step_data.step:
                ref_mass = step_data[step_data.step==reference_mass_name].sample_mass
            else:
                logger.error(f"{reference_mass_name!r} is no valid option. Choose one of {step_data.step.to_list()!r}")
                continue
        else:
            ref_mass = sample.reference_mass

        label = make_title(sample)
        if plot == "TG":
            x = copy.deepcopy(sample.tga[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = sample.tga["sample_mass"]
            elif y_axis == "rel":
                y = 100 * sample.tga["sample_mass"] / ref_mass
            if (
                ylim == "auto"
            ):  # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(
                x,
                y,
                linewidth=linewidth,
                label=label,
            )
        if plot == "DTG":
            x = copy.deepcopy(sample.tga[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = sample.tga["dtg"] * 60
            elif y_axis == "rel":
                y = sample.tga["dtg"] * 60 / ref_mass * 100
            if (
                ylim == "auto"
            ):  # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(
                x,
                y,
                linewidth=linewidth,
                label=label)
        if plot == "heat_flow":
            x = copy.deepcopy(sample.tga[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = sample.tga["heat_flow"]
            elif y_axis == "rel":
                y = sample.tga["heat_flow"] / ref_mass
            if (
                ylim == "auto"
            ):  # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(
                x,
                y,
                linewidth=linewidth,
                label=label)
        if plot == "EGA":
            x = copy.deepcopy(sample.ega[x_axis])
            if x_axis == "time":
                x /= 60
            if y_axis == "orig":
                y = sample.ega[gas]
                if (
                    ylim == "auto"
                ):  # only select relevant range of x data, to auto-scale the y axis
                    x, y, ylim_temp = ylim_auto(x, y, xlim)
                ax.plot(
                    x,
                    y,
                    linewidth=linewidth,
                    label=label)
            elif y_axis == "rel":
                y = sample.ega[gas] / sample.linreg["slope"][gas] / ref_mass
                if (
                    ylim == "auto"
                ):  # only select relevant range of x data, to auto-scale the y axis
                    x, y, ylim_temp = ylim_auto(x, y, xlim)
                ax.plot(
                    x,
                    y,
                    linewidth=linewidth,
                    label=label)

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


