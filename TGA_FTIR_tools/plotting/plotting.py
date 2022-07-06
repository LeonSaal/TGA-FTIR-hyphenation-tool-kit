import os
from turtle import color
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import copy
from ..config import PATHS, SEP, UNITS, PARAMS, MOLAR_MASS, SAVGOL

from .utils import ylim_auto, get_label


WINDOW_LENGTH = int(SAVGOL.getfloat("window_length"))
POLYORDER =int( SAVGOL.getfloat("POLYORDER"))

import logging


logger = logging.getLogger(__name__)


def plot_TGA(
    sample,
    plot,
    save=False,
    x_axis="sample_temp",
    y_axis="orig",
    ylim="auto",
    xlim=[None, None],
    legend=True,
    title=True,
):
    "plot TG data"

    # setting up plot
    fig, TGA = plt.subplots()

    fig.subplots_adjust(right=0.8)

    DTG = TGA.twinx()

    x = copy.deepcopy(sample.tga[x_axis])

    # if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis

    # adjusting y data and setting axis labels according to y_axis
    if y_axis == "rel":
        y = (sample.tga[plot] / sample.info[sample.info["reference_mass"]]) * 100
        yDTG = sample.tga["dtg"] * 60 / sample.info[sample.info["reference_mass"]] * 100
        ylabelDTG = rf'{PARAMS["dtg"]} {SEP} $ \%\,min^{{-1}}$'
        if plot == "sample_mass":
            ylabel = f"{PARAMS['sample_mass']} {SEP} $\\%$"
        elif plot == "heat_flow":
            ylabel = f'{PARAMS["heat_flow"]} {SEP} $ {UNITS["heat_flow"]}\,{UNITS["sample_mass"]}^{{-1}}$'

    elif y_axis == "orig":
        y = sample.tga[plot]
        yDTG = sample.tga["dtg"] * 60  # turning dtg from mg/s in mg/min
        ylabelDTG = f"{PARAMS['dtg']} {SEP} ${UNITS['sample_mass']}\\,min^{{-1}}$"
        ylabel = f"{PARAMS[plot]} {SEP} ${UNITS[plot]}$"

    TGA.set_xlabel(f"{PARAMS[x_axis.lower()]} {SEP} ${UNITS[x_axis.lower()]}$")

    # adjusting x data if x_axis == time and constructing y-axis for temperature
    if x_axis == "time":
        x = x / 60
        temp = TGA.twinx()
        temp.plot(
            x,
            sample.tga["sample_temp"],
            label=f"{PARAMS['sample_temp']} {SEP} ${UNITS['sample_temp']}$",
            ls="dashed",
            color="black",
        )
        temp.spines["right"].set_position(("axes", 1.15))
        temp.set_ylabel(f"{PARAMS['sample_temp']} {SEP} ${UNITS['sample_temp']}$")
        if legend:
            temp.legend(loc=1)

    if ylim == "auto":  # only select relevant range of x data, to auto-scale the y axis
        x, y, ylim = ylim_auto(x, y, xlim)
        x, yDTG, ylim = ylim_auto(x, yDTG, xlim)

    # actual plotting
    (gTGA,) = TGA.plot(x, y, ls="-", color='C0',label=ylabel)
    (gDTG,) = DTG.plot(x, yDTG, ls="--", color='C1',label="DTG")

    TGA.set_ylabel(ylabel)
    TGA.set_ylim(ylim)
    TGA.set_xlim(xlim)
    DTG.set_ylabel(ylabelDTG)

    TGA.yaxis.label.set_color(gTGA.get_color())
    DTG.yaxis.label.set_color(gDTG.get_color())

    DTG.set_yticks(np.linspace(DTG.get_yticks()[0], DTG.get_yticks()[-1], len(TGA.get_yticks())))
    TGA.set_yticks(np.linspace(TGA.get_yticks()[0], TGA.get_yticks()[-1], len(DTG.get_yticks())))
    if title == True:
        TGA.set_title(
            f"{sample.info['alias']}, {sample.info['reference_mass']} = {sample.info[sample.info['reference_mass']]:.2f} ${UNITS['sample_mass']}$"
        )

    plt.show()

    if save:
        path_plots_tga = os.path.join(PATHS["plots"], "TGA")
        if os.path.exists(path_plots_tga) == False:
            os.makedirs(path_plots_tga)
        fig.savefig(
            os.path.join(path_plots_tga, f"{sample.info['name']}_TG_{y_axis}.png"),
        )


def plot_FTIR(
    sample,
    save=False,
    gases=[],
    x_axis="sample_temp",
    y_axis="orig",
    xlim=[None, None],
    legend=True,
    title=True,
):
    "plot IR data"
    gases = set([gas.upper() for gas in gases])
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    x = copy.deepcopy(sample.ir[x_axis])

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
                logger.warn(
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
                logger.warn(
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
    fig, g0 = plt.subplots()
    fig.subplots_adjust(right=0.8)
    graphs.append(g0)
    graphs[0].set_xlabel(f"{PARAMS[x_axis.lower()]} {SEP} ${UNITS[x_axis.lower()]}$")

    if x_axis == "time":
        x = x / 60
    if y_axis == "orig":
        graphs[0].set_ylabel(f"{get_label(gases[0])} {SEP} ${UNITS['ir']}$")
        graphs[0].yaxis.label.set_color(colors[0])
    elif y_axis == "rel":
        graphs[0].set_ylabel(
            f"${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}\\,{UNITS['time']}^{{-1}}$"
        )

    # plot data and append secondary, third... y-axis on right side if necessary
    for i, gas in enumerate(gases):
        if y_axis == "orig":
            y = sample.ir[gas]

            if i > 0:
                graphs.append(graphs[0].twinx())
                graphs[i].spines["right"].set_position(("axes", 1 + (i - 1) * 0.1))
            graphs[i].plot(x, y, color=colors[i])
            graphs[i].set_ylabel(f"{get_label(gas)} {SEP} {UNITS['ir']}")
            graphs[i].yaxis.label.set_color(colors[i])
            graphs[i].set_yticks(np.linspace(graphs[i].get_yticks()[0], graphs[i].get_yticks()[-1], len(graphs[0].get_yticks())))

        elif y_axis == "rel":
            tot_area = np.sum(sample.ir[gas])
            tot_mol = (
                (tot_area - sample.linreg["intercept"][gas])
                / sample.linreg["slope"][gas]
                / sample.info[sample.info["reference_mass"]]
            )
            y = sample.ir[gas] / tot_area * tot_mol
            graphs[0].plot(x, y, label=get_label(gas))

    if legend and y_axis == "rel":
        fig.legend()


    if title == True:
        graphs[0].set_title(
            f"{sample.info['alias']}, {sample.info['reference_mass']}: {sample.info[sample.info['reference_mass']]:.2f}$\\,{UNITS['sample_mass']}$"
        )

    graphs[0].set_xlim(xlim)
    plt.show()

    if save:
        path_plots_ir = os.path.join(PATHS["plots"], "IR")
        if os.path.exists(path_plots_ir) == False:
            os.makedirs(path_plots_ir)
        fig.savefig(
            os.path.join(path_plots_ir, f"{sample.info['name']}_IR_{y_axis}.png"),
        )


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
            logger.warn(
                f"{' and '.join([gas for gas in list(on_axis - calibrated)])} not calibrated. "
            )
            logger.info(
                f'Proceeding with {" and ".join([gas for gas in intersection])}.'
            )
        gases_temp = list(intersection)

    else:
        intersection = calibrated & on_axis & gases_temp
        if len(gases_temp - calibrated) != 0:
            logger.warn(
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
        sample.tga, sample.ir, how="left", on=["time", "sample_temp"]
    ).dropna()
    DTG = -sp.signal.savgol_filter(data["sample_mass"], 13, 3, deriv=1)

    x = data[x_axis]
    y = np.zeros((len(gases), len(sample.ir)))
    out = pd.DataFrame()
    out["time"] = data["time"]
    out["sample_temp"] = data["sample_temp"]
    cumul = np.zeros(len(sample.ir))

    # calibrating IR data
    for i, gas in enumerate(gases):
        tot_area = np.sum(sample.ir[gas])
        tot_mass = (
            (tot_area - sample.linreg["intercept"][gas])
            / sample.linreg["slope"][gas]
            * MOLAR_MASS.getfloat(gas)
        )
        y[i][:] = sample.ir[gas] / tot_area * tot_mass
        out[gas] = y[i][:]
        cumul += y[i][:]

    # setup figure
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(8, 1)
    stack = fig.add_subplot(gs[:-1, 0])

    if title == True:
        stack.set_title(
            f"{sample.info['alias']}, {sample.info['reference_mass']} = {sample.info[sample.info['reference_mass']]:.2f} ${UNITS['sample_mass']}$"
        )

    error = fig.add_subplot(gs[-1, 0], sharex=stack)

    stack.set_xlabel(f"{PARAMS[x_axis.lower()]} {SEP} ${UNITS[x_axis.lower()]}$")
    error.set_xlabel(f"{PARAMS[x_axis.lower()]} {SEP} ${UNITS[x_axis.lower()]}$")

    # actual plotting
    if x_axis == "time":
        x = x / 60
        temp = stack.twinx()
        temp.plot(x, data["sample_temp"], ls="dashed", color="black", label="T")
        temp.set_ylabel(f"{PARAMS['sample_temp']} {SEP} ${UNITS['sample_temp']}$")
        if legend:
            temp.legend()

    stack.stackplot(x, y, labels=[get_label(gas) for gas in gases])
    stack.plot(x, DTG, label=PARAMS["dtg"])
    stack.set_ylabel(
        f"{PARAMS['dtg']}, {', '.join([get_label(gas) for gas in gases])} {SEP} ${UNITS['sample_mass']}\\,{UNITS['time']}^{{-1}}$"
    )



    if legend:
        stack.legend()

    # plot error of reconstruction
    error.plot(x, DTG - cumul)
    error.hlines(0, min(x), max(x), ls="dashed")
    error.set_ylabel(f"$\\Delta$ {PARAMS['dtg']}")
    stack.set_xlim(xlim)
    plt.show()

    if save:
        path_plot_irdtg = os.path.join(PATHS["plots"], "IRDTG")
        if os.path.exists(path_plot_irdtg) == False:
            os.makedirs(path_plot_irdtg)
        fig.savefig(os.path.join(path_plot_irdtg, f"{sample.info['name']}_IRDTG.png"),)
        out["dtg"] = DTG
        out.to_excel(os.path.join(PATHS["output"], sample.info["name"] + "_IRDTG.xlsx"))

