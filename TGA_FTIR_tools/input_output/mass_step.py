from ..config import SAVGOL, UNITS, PARAMS, SEP
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import matplotlib.ticker as ticker


def mass_step(TGA_data, plot=False, **kwargs):  # rel_height=.963
    "deriving mass steps via peaks in DTG signal"
    # calculation and smoothing of DTG
    TG = TGA_data["sample_mass"] / TGA_data["sample_mass"][0]
    DTG = sp.signal.savgol_filter(
        TG, int(SAVGOL["window_length"]), int(SAVGOL["polyorder"]), deriv=1
    )

    # detect mass steps
    peaks, properties = sp.signal.find_peaks(-DTG, **kwargs)
    step_end = properties["right_ips"].astype(np.int, copy=False)
    step_start = properties["left_ips"].astype(np.int, copy=False)

    # calculate masssteps
    steps = np.zeros(len(peaks))
    samples = 20
    for i in range(len(peaks)):
        steps[i] = np.mean(TGA_data["sample_mass"][step_end[i] : step_end[i] + samples])

    # calculate step height
    step_height = np.zeros(len(steps))
    steps = np.insert(steps, 0, TGA_data["sample_mass"][0])

    step_height = np.zeros(len(step_start))
    for i in range(len(step_start)):
        step_height[i] = np.mean(
            TGA_data["sample_mass"][step_start[i] - samples : step_start[i]]
        ) - np.mean(TGA_data["sample_mass"][step_end[i] : step_end[i] + samples])

    rel_step_height = step_height / steps[0] * 100
    # plotting
    if plot == True:
        # plotting of rel. TG
        x = TGA_data["sample_temp"]
        fig, ax = plt.subplots()
        rel_steps = steps / steps[0] * 100
        ax.hlines(
            rel_steps[:-1],
            np.zeros(len(rel_steps) - 1),
            x[step_end],
            linestyle="dashed",
        )
        ax.vlines(x[step_end], rel_steps[1:], rel_steps[:-1], linestyle="dashed")
        for i in range(len(step_end)):
            ax.text(
                x[step_end[i]] + 5,
                rel_steps[i + 1] + rel_step_height[i] / 2,
                str(round(rel_step_height[i], 2)) + " %",
            )
        ax.plot(x, TGA_data["sample_mass"] / TGA_data["sample_mass"][0] * 100)
        ax.text(
            0.85 * max(TGA_data["sample_temp"]),
            100,
            f'sample mass: {TGA_data["sample_mass"][0]:.2f} {UNITS["sample_mass"]}',
            horizontalalignment="center",
        )
        ax.set_xlabel(f'{PARAMS["sample_temp"]} {SEP} {UNITS["sample_temp"]}')
        ax.set_ylabel(f'{PARAMS["sample_mass"]} {SEP} %')
        ax.xaxis.set_minor_locator(
            ticker.AutoMinorLocator()
        )  # switch on minor ticks on each axis
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.set(title="TG")
        plt.show()

        # plotting of DTG
        fig, ax = plt.subplots()
        y = -DTG
        ax.plot(x, y)
        ax.vlines(x[step_end], 0, max(y), linestyle="dashed")
        ax.vlines(x[step_start], 0, max(y), linestyle="dashed")
        ax.vlines(x[peaks], y[peaks] - properties["prominences"], y[peaks])
        ax.hlines(
            y[peaks] - kwargs["rel_height"] * properties["prominences"],
            x[step_end],
            x[step_start],
        )
        ax.set_xlabel(f'{PARAMS["sample_temp"]} {SEP} {UNITS["sample_temp"]}')
        ax.set_ylabel(
            f'{PARAMS["dtg"]} { SEP} {UNITS["sample_mass"]} ${UNITS["time"]}^{{-1}}$'
        )
        ax.xaxis.set_minor_locator(
            ticker.AutoMinorLocator()
        )  # switch on minor ticks on each axis
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.set(title="DTG")
        plt.show()

    return step_height, rel_step_height, step_start, step_end
