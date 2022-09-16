import matplotlib.pyplot as plt
import pandas as pd

from ..config import PARAMS, SEP, UNITS
from .utils import get_label


def plot_corr(originalData: pd.DataFrame, BaselineData:pd.DataFrame, label: str) -> None:
    try:
        x = originalData["sample_temp"]
    except:
        x = originalData["time"]
        x /= 60

    _, ax = plt.subplots()
    ax.plot(x, originalData[label], label="data")
    ax.plot(x, BaselineData[label], label="baseline")
    ax.plot(x, originalData[label].subtract(BaselineData[label]), label="corr. data")

    ax.legend()
    if x.name == "time":
        ax.set_xlabel(f"{PARAMS['time']} {SEP} {UNITS['time']}")
    elif x.name == "sample_temp":
        ax.set_xlabel(f"{PARAMS['sample_temp']} {SEP} {UNITS['sample_temp']}")
    ax.set_ylabel(f"{get_label(label)} {SEP} {UNITS['IR' if label not in UNITS else label]}")
    ax.set(title=f"{get_label(label)} baseline correction")
    plt.show()
