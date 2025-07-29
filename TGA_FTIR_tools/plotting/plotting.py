import copy
import os
from turtle import color

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from chempy import Substance

from ..config import MERGE_CELLS, PATHS, SAVGOL, SEP, UNITS
from .utils import get_label, make_title, ylim_auto

WINDOW_LENGTH = int(SAVGOL.getfloat("window_length"))
POLYORDER =int( SAVGOL.getfloat("POLYORDER"))

import logging

logger = logging.getLogger(__name__)


def plot_TGA(
    sample,
    plot,
    ax,
    x_axis="sample_temp",
    y_axis="orig",
    ylim="auto",
    xlim=[None, None],
    legend=True,
    title=True,
):
    "plot TG data"

    DTG = ax.twinx()

    x = copy.deepcopy(sample.tga[x_axis])

    # if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis

    # adjusting y data and setting axis labels according to y_axis
    if y_axis == "rel":
        y = (sample.tga[plot] / sample.reference_mass) * 100
        yDTG = sample.tga["dtg"] * 60 / sample.reference_mass * 100
        ylabelDTG = rf'{get_label("dtg")} {SEP} $ \%\,min^{{-1}}$'
        if plot == "sample_mass":
            ylabel = f"{get_label('sample_mass')} {SEP} $\\%$"
        elif plot == "heat_flow":
            ylabel = f'{get_label("heat_flow")} {SEP} $ {UNITS["heat_flow"]}\,{UNITS["sample_mass"]}^{{-1}}$'

    elif y_axis == "orig":
        y = sample.tga[plot]
        yDTG = sample.tga["dtg"] * 60  # turning dtg from mg/s in mg/min
        ylabelDTG = f"{get_label('dtg')} {SEP} ${UNITS['sample_mass']}\\,min^{{-1}}$"
        ylabel = f"{get_label(plot)} {SEP} ${UNITS[plot]}$"

    ax.set_xlabel(f"{get_label(x_axis.lower())} {SEP} ${UNITS[x_axis.lower()]}$")

    # adjusting x data if x_axis == time and constructing y-axis for temperature
    if x_axis == "time":
        x = x / 60
        temp = ax.twinx()
        temp.plot(
            x,
            sample.tga["sample_temp"],
            label=f"{get_label('sample_temp')} {SEP} ${UNITS['sample_temp']}$",
            ls="dashed",
            color="black",
        )
        temp.spines["right"].set_position(("axes", 1.15))
        temp.set_ylabel(f"{get_label('sample_temp')} {SEP} ${UNITS['sample_temp']}$")
        if legend:
            temp.legend(loc=1)

    if ylim == "auto":  # only select relevant range of x data, to auto-scale the y axis
        x, y, ylim = ylim_auto(x, y, xlim)
        x, yDTG, ylim = ylim_auto(x, yDTG, xlim)

    # actual plotting
    (gTGA,) = ax.plot(x, y, ls="-", color='C0',label=ylabel)
    (gDTG,) = DTG.plot(x, yDTG, ls="--", color='C1',label="DTG")

    ax.set_ylabel(ylabel)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    DTG.set_ylabel(ylabelDTG)

    ax.yaxis.label.set_color(gTGA.get_color())
    DTG.yaxis.label.set_color(gDTG.get_color())

    DTG.set_yticks(np.linspace(DTG.get_yticks()[0], DTG.get_yticks()[-1], len(ax.get_yticks())))
    ax.set_yticks(np.linspace(ax.get_yticks()[0], ax.get_yticks()[-1], len(DTG.get_yticks())))
    if title == True:
        ax.set_title(make_title(sample))


def plot_FTIR(
    sample,
    ax,
    gases=[],
    x_axis="sample_temp",
    y_axis="orig",
    xlim=[None, None],
    legend=True,
    title=True,
):
    "plot IR data"
    #gases = set([gas.upper() for gas in gases])
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    x = copy.deepcopy(sample.ega[x_axis])

    # catching possible input errors
    if sample.linreg is not None:
        calibrated = set(sample.linreg.index)
    else:
        calibrated = set()
    on_axis = set(sample.info["gases"])
    if len(gases) == 0:
        if y_axis == "rel":
            intersection = calibrated & on_axis
            if len(on_axis - calibrated) != 0:
                logger.info(
                    f"{' and '.join([gas for gas in list(on_axis - calibrated)])} not calibrated. Proceeding with {' and '.join([gas for gas in intersection])}."
                )
            gases = intersection
        elif y_axis == "orig":
            gases = on_axis
    else:
        if y_axis == "rel":
            gases = set(gases)
            intersection = calibrated & on_axis & gases
            if len(gases - calibrated) != 0:
                logger.warning(
                    f"{' and '.join([gas for gas in gases - calibrated])} not calibrated."
                )
            if len(intersection) != 0:
                logger.info(
                    f"Proceeding with {' and '.join([gas for gas in intersection])}."
                )
                gases = intersection
            else:
                logger.error(
                    "None of supplied gases was found in IR data and calibrated."
                )
                return
        elif y_axis == "orig":
            gases = set(gases)
            intersection = on_axis & gases
            if len(gases - on_axis) != 0:
                logger.warning(
                    f"{' and '.join([gas for gas in gases - on_axis])} not found in IR data."
                )
            if len(intersection) != 0:
                logger.info(
                    f"Proceeding with {' and '.join([gas for gas in intersection])}."
                )
                gases = intersection
            else:
                logger.error("None of supplied gases was found in IR data.")
                return
    gases = list(gases)

    # setup figure
    graphs = []
    graphs.append(ax)
    graphs[0].set_xlabel(f"{get_label(x_axis.lower())} {SEP} ${UNITS[x_axis.lower()]}$")

    if x_axis == "time":
        x = x / 60
    if y_axis == "orig":
        graphs[0].set_ylabel(f"{get_label(gases[0])} {SEP} ${UNITS['ega']}$")
        graphs[0].yaxis.label.set_color(colors[0])
    elif y_axis == "rel":
        graphs[0].set_ylabel(
            f"${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}\\,{UNITS['time']}^{{-1}}$"
        )

    # plot data and append secondary, third... y-axis on right side if necessary
    for i, gas in enumerate(gases):
        if y_axis == "orig":
            y = sample.ega[gas]

            if i > 0:
                graphs.append(graphs[0].twinx())
                graphs[i].spines["right"].set_position(("axes", 1 + (i - 1) * 0.1))
            graphs[i].plot(x, y, color=colors[i])
            graphs[i].set_ylabel(f"{get_label(gas)} {SEP} {UNITS['ega']}")
            graphs[i].yaxis.label.set_color(colors[i])
            graphs[i].set_yticks(np.linspace(graphs[i].get_yticks()[0], graphs[i].get_yticks()[-1], len(graphs[0].get_yticks())))

        elif y_axis == "rel":
            tot_area = np.sum(sample.ega[gas])
            tot_mol = (
                (tot_area - sample.linreg["intercept"][gas])
                / sample.linreg["slope"][gas]
                / sample.reference_mass
            )
            y = sample.ega[gas] / tot_area * tot_mol
            graphs[0].plot(x, y, label=get_label(gas))

    if legend and y_axis == "rel":
        ax.legend()


    if title == True:
        graphs[0].set_title(make_title(sample))

    graphs[0].set_xlim(xlim)

def FTIR_to_DTG(
    sample,
    x_axis="sample_temp",
    save=False,
    gases=[],
    legend=True,
    y_axis=None,
    xlim=[None, None],
    title=True,
):
    "reconstructing DTG from calibrated IR data"
    gases_temp = set([gas.upper() for gas in gases])

    # catching possible input errors
    try:
        calibrated = set(sample.linreg.index)
    except AttributeError:
        calibrated = set()
    on_axis = set(sample.info["gases"])
    if len(gases) == 0:
        intersection = calibrated & on_axis
        if calibrated != on_axis:
            logger.warning(
                f"{' and '.join([gas for gas in list(on_axis - calibrated)])} not calibrated. "
            )
            logger.info(
                f'Proceeding with {" and ".join([gas for gas in intersection])}.'
            )
        gases_temp = list(intersection)

    else:
        intersection = calibrated & on_axis & gases_temp
        if len(gases_temp - calibrated) != 0:
            logger.warning(
                "{} not calibrated.".format(
                    " and ".join([gas for gas in (gases_temp - calibrated)])
                )
            )
        if len(intersection) != 0:
            logger.info(
                f"Proceeding with {' and '.join([gas for gas in intersection])}."
            )
            gases_temp = list(intersection)
            gases_temp.sort(key=lambda i: gases.index(i) if i in gases else len(gases))
        else:
            logger.error("None of supplied gases was found in IR data and calibrated.")
            return

    # setup IR data and calculating DTG
    gases = gases_temp
    data = pd.merge(
        sample.tga, sample.ega, how="left", on=["time", "sample_temp"]
    ).dropna()
    DTG = -sp.signal.savgol_filter(data["sample_mass"], 13, 3, deriv=1)

    x = data[x_axis]
    y = np.zeros((len(gases), len(sample.ega)))
    out = pd.DataFrame()
    out["time"] = data["time"]
    out["sample_temp"] = data["sample_temp"]
    cumul = np.zeros(len(sample.ega))

    # calibrating IR data
    for i, gas in enumerate(gases):
        tot_area = np.sum(sample.ega[gas])
        tot_mass = (
            (tot_area - sample.linreg["intercept"][gas])
            / sample.linreg["slope"][gas]
            * sample.linreg["molmass"][gas]
        )
        y[i][:] = sample.ega[gas] / tot_area * tot_mass
        out[gas] = y[i][:]
        cumul += y[i][:]

    # setup figure
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(8, 1)
    stack = fig.add_subplot(gs[:-1, 0])

    if title == True:
        stack.set_title(make_title(sample))

    error = fig.add_subplot(gs[-1, 0], sharex=stack)

    stack.set_xlabel(f"{get_label(x_axis.lower())} {SEP} ${UNITS[x_axis.lower()]}$")
    error.set_xlabel(f"{get_label(x_axis.lower())} {SEP} ${UNITS[x_axis.lower()]}$")

    # actual plotting
    if x_axis == "time":
        x = x / 60
        temp = stack.twinx()
        temp.plot(x, data["sample_temp"], ls="dashed", color="black", label="T")
        temp.set_ylabel(f"{get_label('sample_temp')} {SEP} ${UNITS['sample_temp']}$")
        if legend:
            temp.legend()

    stack.stackplot(x, y, labels=[get_label(gas) for gas in gases])
    stack.plot(x, DTG, label=get_label("dtg"))
    stack.set_ylabel(
        f"{get_label('dtg')}, {', '.join([get_label(gas) for gas in gases])} {SEP} ${UNITS['sample_mass']}\\,{UNITS['time']}^{{-1}}$"
    )



    if legend:
        stack.legend()

    # plot error of reconstruction
    error.plot(x, DTG - cumul)
    error.hlines(0, min(x), max(x), ls="dashed")
    error.set_ylabel(f"$\\Delta$ {get_label('dtg')}")
    stack.set_xlim(xlim)

    if save:
        path_plot_irdtg = PATHS["plots"]/ "IRDTG"
        if not path_plot_irdtg.exists():
            path_plot_irdtg.mkdir()
        fig.savefig(path_plot_irdtg/ f"{sample.name}_IRDTG.png")
        out["dtg"] = DTG
        out.to_excel(PATHS["output"]/ f"{sample.name}_IRDTG.xlsx", merge_cells=MERGE_CELLS)

