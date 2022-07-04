from .plot_fit import plot_fit
from .plotting import get_label, plot_FTIR, plot_TGA, FTIR_to_DTG
from .plots import plots
from .plot_rob import plot_robustness
from .plot_results import bar_plot_results
from .plot_dweight import plot_dweight
from ..config import STYLE

import matplotlib.pyplot as plt
plt.style.use(STYLE)
