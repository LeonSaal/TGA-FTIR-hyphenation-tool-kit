from typing import Literal
import matplotlib.pyplot as plt
import numpy as np
from .utils import make_title
import pandas as pd
import pint
import pint_pandas
ureg = pint.get_application_registry()


def plot_mass_steps(sample, ax:plt.Axes, steps = [], y_axis:Literal['rel','orig']='rel', x_axis='sample_temp', title=False, thresh=.025):
    """
    steps: (start, stop)
    """
    # get steps
    # if 0 not in steps:
    #     steps = [0, *steps]
    matches = sample.get_value(*steps, at=x_axis)
    step_masses = matches.sample_mass
    steps = matches[x_axis]
    step_heights = step_masses.diff()
    
    # set up plot
    x,y = sample.tga[x_axis], sample.tga.sample_mass

    #plot data and set axes labels
    ax.plot(x, y)

    # plot steps
    ax.hlines(
        step_masses[:-1],
        np.zeros(len(step_masses)-1),
        steps[1:],
        linestyle="dashed",
    )
    ax.vlines(steps[1:], step_masses[1:], step_masses[:-1], linestyle="dashed")
    
    #plot step heights
    for step_end, step_startmass, step_endmass, step_height in zip(steps[1:],step_masses[:-1], step_masses[1:], step_heights[1:]):
        y_mean = (step_startmass + step_endmass)/2
        percent = (step_height / step_masses[0]).to("percent")
        if abs(percent) < ureg.Quantity(thresh*100, "percent"):
            continue
        if y_axis =='orig':
            label = f'{step_height:.2f~P}({percent:.2f~P})'
        else:
            label = f'{step_height/y.max():.1%}'
        ax.text(step_end+ureg.Quantity(5, ureg.delta_degC), y_mean ,label)

    

    if title:
        ax.set_title(make_title(sample))