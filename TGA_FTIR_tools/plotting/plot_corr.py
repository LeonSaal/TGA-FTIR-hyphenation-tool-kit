import matplotlib.pyplot as plt
import pandas as pd
from types import NoneType

from ..config import SEP, UNITS
from .utils import get_label


def plot_corr(originalData: pd.DataFrame, BaselineData:pd.DataFrame, label: str) -> NoneType:
    try:
        x = originalData["sample_temp"]
    except:
        x = originalData["time"]
        x = x.pint.to(UNITS.get("time"))

    _, ax = plt.subplots()
    ax.plot(x, originalData[label], label="data")
    ax.plot(x, BaselineData[label], label="baseline")
    ax.plot(x, originalData[label].subtract(BaselineData[label]), label="corr. data")

    ax.legend()
    ax.set_xlabel(f"{get_label(x.name)} {SEP} {UNITS.get(x.name, '?')}")
    ax.set_ylabel(f"{get_label(label)} {SEP} {UNITS.get('EGA' if label not in UNITS else label, '?')}")
    ax.set(title=f"{get_label(label)} baseline correction")

