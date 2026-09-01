import copy
import logging
import os
from typing import Literal

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from ..config import PATHS, SEP, UNITS
from ..input_output.general import time
from .utils import get_label, make_title, ylim_auto, _validate_lim

import pint

logger = logging.getLogger(__name__)

ureg = pint.get_application_registry()



def plots(
    samples,
    plot: Literal["TG", "EGA", "DTG", "heat_flow"],
    ax=None,
    x_axis="sample_temp",
    y_axis=Literal["orig", "rel","rel_mol"],
    ylim="auto",
    xlim=(None, None),
    gas=None,
    legend=True,
    reference_mass_name=None,
    linewidth=1,
    **kwargs,
):
    "overlay plots from different samples"
    options = ["TG", "EGA", "DTG", "heat_flow"]
    if plot not in options:
        logger.warning(f"{plot=} not in {options}")
        return

    # setting up axis-labels and catching possible input errors
    if plot == "EGA":
        avail = {gas for sample in samples for gas in sample.ega.columns.to_list()}
        if gas == None:
            logger.warning(f"Supply 'gas = '. Available gases: {avail}.")
            return
        else:
            if gas not in avail:
                logger.warning(f"{gas} was not found in EGA data. Available gases: {avail}")
                return

        # just to see if supplied gas is calibrated or not
        calibrated = set()
        for sample in samples:
            try:
                calibrated.update(set(sample.linreg.index))
            except AttributeError:
                if y_axis == "rel_mol":
                    logger.warning(f"{gas} is not calibrated for {sample.name}")

        if y_axis == "rel_mol":
            if gas not in calibrated:
                logger.warning(
                    f"{gas} is not calibrated. Changing y_axis to 'orig'..."
                )
                y_axis = "orig"



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
        if x_axis not in sample.tga.columns and x_axis not in sample.ega.columns:
            logger.warning(f"{x_axis!r} not found in data.")
            continue

        # scale to default
        x = copy.deepcopy(sample.tga[x_axis])
        if x_axis == "time":
            x = x.pint.to(UNITS.get("time"))

        # select data to plot
        match plot:
            case "TG":
                y = sample.tga["sample_mass"]
            case "DTG":
                y = sample.tga["dtg"].pint.to("mg/min")
            case "heat_flow":
                y = sample.tga["heat_flow"]
            case "EGA":
                if gas not in sample.ega.columns:
                    logger.warning(f"{gas} was not found in IR data for {sample.name}.")
                    continue
                y = sample.ega[gas]

        # scale
        orig_units = y.pint.units
        if y_axis == "rel":
            y = y / ref_mass * 100
        elif y_axis == "rel_mol" and plot == "EGA":
            y = y / sample.linreg["slope"][gas] / ref_mass

        # get units or percent
        units = "%" if ((y.pint.units == "") and (orig_units != "")) else y.pint.units

        # only select relevant range of x data, to auto-scale the y axis
        x, y, ylim_temp = ylim_auto(x, y, xlim) if ylim == "auto" else (x, y, ylim)

        ax.plot(
            x,
            y,
            linewidth=linewidth,
            label=label)

        ax.set_xlim(_validate_lim(xlim, x))
        ax.set_ylim(_validate_lim(ylim_temp, y))

    # set up axes labels
    match plot:
        case "TG":
            ylabel = "sample_mass"      
        case "EGA": 
            ylabel = gas
        case _ :
            ylabel = plot.lower()
    ax.set_ylabel(f"{get_label(ylabel)} {SEP} {units}")

    if legend:
        ax.legend(
                loc="center left",
                bbox_to_anchor=(1, 0.5),
                frameon=False
            )

    # switch on minor ticks on each axis
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())  
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())


