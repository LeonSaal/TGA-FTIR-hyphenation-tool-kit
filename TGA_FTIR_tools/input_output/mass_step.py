import numpy as np
import scipy as sp

from ..config import SAVGOL


def mass_step(sample, samples= 20, **kwargs):  # rel_height=.963
    "deriving mass steps via peaks in DTG signal"
    # calculation and smoothing of DTG
    x = sample.tga['sample_temp'].to_numpy()
    y = sample.tga['sample_mass'].to_numpy()

    TG = y / y[0]
    DTG = sp.signal.savgol_filter(
        TG,
        int(SAVGOL.getfloat("window_length")),
        int(SAVGOL.getfloat("polyorder")),
        deriv=1,
    )
    
    # detect mass steps
    _, properties = sp.signal.find_peaks(-DTG, **kwargs)
    step_starts = np.insert(properties["left_ips"].astype(np.int, copy=False),0,0)
    step_ends = np.append(properties["right_ips"].astype(np.int, copy=False),len(x))

    # calculate mass steps
    start_masses = np.zeros(len(step_starts))
    for i, step_start in enumerate(step_starts):
        start_masses[i] = np.mean(y[step_start : step_start + samples])

    # calculate step height
    step_height = np.zeros(len(step_starts))
    for i, (start_mass, step_end) in enumerate(zip(start_masses, step_ends)):
        step_height[i] = start_mass - np.mean(y[step_end - samples : step_end])

    rel_step_height = step_height / sample.reference_mass

    return step_height, rel_step_height, step_starts, step_ends
