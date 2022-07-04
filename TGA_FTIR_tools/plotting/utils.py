import pandas as pd
from ..config import LABELS


def ylim_auto(x, y, xlim):
    "truncate x and y according to xlim"
    x_min = xlim[0]
    x_max = xlim[1]
    if pd.isnull(xlim[0]):
        x_min = x.min()
    if pd.isnull(xlim[1]):
        x_max = x.max()
    x = x[(x >= x_min) & (x <= x_max)]
    y = y[x.index]
    ylim = [None, None]  # reset ylim

    return x, y, ylim


def get_label(key):
    "get labels to put in plots"
    if key in LABELS:
        return LABELS[key]
    elif key.isdigit():
        if int(key) in LABELS:
            return LABELS[int(key)]
    return str(key)
