import matplotlib.pyplot as plt

from ..config import STYLE
from .plot_corr import plot_corr
from .plot_dweight import plot_dweight
from .plot_fit import plot_fit
from .plot_results import bar_plot_results
from .plot_rob import plot_robustness
from .plots import plots
from .plotting import FTIR_to_DTG, get_label, plot_FTIR, plot_TGA

plt.style.use(STYLE)
