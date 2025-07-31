from typing import Literal
import matplotlib.pyplot as plt
import numpy as np
from .utils import make_title
import pint
ureg = pint.get_application_registry()

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
            percent = (step_height / step_masses[0]).to("percent")
            label = f'{step_height:.2f}({percent:.2f})'
        else:
            label = f'{step_height/y.max():.1%}'
        ax.text(step_end+ureg.Quantity(5, ureg.delta_degC), y_mean ,label)

    #plot data and set axes labels
    ax.plot(x, y)

    if title:
        ax.set_title(make_title(sample))