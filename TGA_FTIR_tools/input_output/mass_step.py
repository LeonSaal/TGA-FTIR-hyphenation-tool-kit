import numpy as np
import scipy as sp
import itertools as it
from ..config import SAVGOL
import pandas as pd

import pint
ureg = pint.get_application_registry()

def mass_step(sample, samples= 20, **kwargs):  # rel_height=.963
    "deriving mass steps via peaks in DTG signal"
    # calculation and smoothing of DTG
    y = sample.tga['sample_mass']
    dtg = sample.tga.dtg
    
    # detect mass steps
    peaks_idx, _ = sp.signal.find_peaks(dtg, **kwargs)
    peaks = [0]+ peaks_idx.tolist()+[dtg.size-1]

    # find bounds of peaks
    step_starts_idx = []
    step_ends_idx =  []
    
    # following data from peaks down in both directions and detect change in direction
    for i, (a,b) in enumerate(it.pairwise(peaks)):
        subset = dtg.diff()[a+1:b-1]
        subleft = subset[::-1] > 0
        left = b-np.argmin(subleft) if not subleft.all() else a
        right = np.argmax(subset>0)+a
        if i < len(peaks)-2:
            step_starts_idx.append(left)
        if i !=0:
            step_ends_idx.append(right)

    # calculate mass steps
    start_masses = pd.Series(np.zeros(len(step_starts_idx)), dtype=y.dtype)
    for i, step_start in enumerate(step_starts_idx):
        start_masses[i] = y[step_start : step_start + samples].mean()

    # calculate step height
    step_height = pd.Series(np.zeros(len(step_starts_idx)), dtype=y.dtype)
    for i, (start_mass, step_end) in enumerate(zip(start_masses, step_ends_idx)):
        step_height[i] = start_mass - y[step_end - samples : step_end].mean()

    rel_step_height = step_height / sample.reference_mass

    return step_height, rel_step_height, step_starts_idx, step_ends_idx, peaks_idx
