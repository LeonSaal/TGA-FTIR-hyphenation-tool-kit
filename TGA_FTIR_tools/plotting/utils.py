
import pandas as pd
import re
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

def format_mf_latex(mf: str) -> str:
    """_summary_

    Args:
        mf (str): _description_

    Returns:
        str: _description_
    """    
    pattern = "([A-Z][a-z]?)(\\d*)?"
    replacement = "\\1_{\\2}"
    return re.sub(pattern, replacement, mf).replace("_{}", "")


def get_label(key:str) -> str:
    """_summary_

    Args:
        key (str): _description_

    Returns:
        str: _description_
    """
    if key in LABELS:
        return LABELS[key]
    elif key.isdigit():
        if int(key) in LABELS:
            return LABELS[int(key)]

    mf_latex = format_mf_latex(key)
    return f'${mf_latex}$'


def make_title(sample):
    alias = sample.alias
    if sample.tga is not None:
        type_ref_mass = sample.info['reference_mass_name']
        reference_mass = sample.reference_mass
        label_ref_mass = get_label(type_ref_mass)
        
        return f"{alias}, {label_ref_mass} = {reference_mass:.2f~P}"
    else:
        return sample.alias
            