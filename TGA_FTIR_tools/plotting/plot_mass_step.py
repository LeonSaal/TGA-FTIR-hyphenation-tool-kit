from typing import Literal

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from ..config import SEP, UNITS
from .plotting import get_label
from .utils import make_title


def plot_mass_steps(sample, ax:plt.Axes, steps = [], y_axis:Literal['rel','orig']='rel', x_axis='sample_temp', title=False):
    """
    steps: (start, stop)
    """
    # get steps
    # if 0 not in steps:
    #     steps = [0, *steps]
    matches = sample.get_value(*steps, at=x_axis).T
    step_masses = matches.sample_mass.values
    steps = matches.index
    step_heights = np.diff(step_masses)
    
    # set up plot
    x,y = sample.tga[x_axis], sample.tga.sample_mass
    if y_axis =='rel':
        unit = '%'
        ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=max(y)))
    else:
        unit = UNITS["sample_mass"]

    # plot steps
    ax.hlines(
        step_masses[:-1],
        np.zeros(len(step_masses)-1),
        steps[1:],
        linestyle="dashed",
    )
    ax.vlines(steps[1:], step_masses[1:], step_masses[:-1], linestyle="dashed")
    

    #plot step heights
    for step_start, step_end, step_height in zip(steps[:-1], steps[1:], step_heights):
        y_mean = np.mean((matches.loc[step_start], matches.loc[step_end]))
        if y_axis =='orig':
            label = f'{step_height:.2f} {unit} ({step_height / step_masses[0]:.2%})'
        else:
            label = f'{step_height/y.max():.1%}'
        ax.text(step_end+ 5, y_mean ,label)

    #plot data and set axes labels
    ax.plot(x, y)
    ax.set_xlabel(f'{get_label("sample_temp")} {SEP} {UNITS["sample_temp"]}')
    ax.set_ylabel(f'{get_label("sample_mass")} {SEP} {unit}')

    if title:
        ax.set_title(make_title(sample))

# def plot_dtg_mass_step(x:pd.Series, y:pd.Series, steps:Tuple, peaks):
#     # plotting of DTG
#     step_starts, step_ends = steps
#     fig, ax = plt.subplots()
#     y = -DTG
#     ax.plot(x, y)
#     ax.vlines(x[step_ends], 0, max(y), linestyle="dashed")
#     ax.vlines(x[step_starts], 0, max(y), linestyle="dashed")
#     ax.vlines(x[peaks], y[peaks] - properties["prominences"], y[peaks])
#     ax.hlines(
#         y[peaks] - kwargs["rel_height"] * properties["prominences"],
#         x[step_ends],
#         x[step_starts],
#     )
#     ax.set_xlabel(f'{get_label("sample_temp")} {SEP} {UNITS["sample_temp"]}')
#     ax.set_ylabel(
#         f'{get_label("dtg")} { SEP} {UNITS["sample_mass"]} ${UNITS["time"]}^{{-1}}$'
#     )
#     ax.set(title="DTG")
#     plt.show()
