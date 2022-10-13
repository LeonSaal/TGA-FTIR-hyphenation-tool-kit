
import pandas as pd
from chempy import Substance
from pyparsing.exceptions import ParseException

from ..config import LABELS, UNITS


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
    try:
        substance = Substance.from_formula(key)
        return f'${substance.latex_name}$'
    except ParseException:
        return str(key)

def make_title(sample):
    alias = sample.alias
    type_ref_mass = sample.info['reference_mass']
    reference_mass = sample.info[type_ref_mass]
    label_ref_mass = get_label(type_ref_mass)
    unit = UNITS['sample_mass']
    
    return f"{alias}, {label_ref_mass} = {reference_mass:.2f} ${unit}$"
            