import matplotlib.pyplot as plt

from ..config import STYLE
from .plot_corr import plot_corr
from .plot_dweight import plot_dweight
from .plot_fit import plot_fit
from .plot_mass_step import plot_mass_steps
from .plot_results import bar_plot_results
from .plot_rob import plot_robustness
from .plots import plots
from .plot_cali import plot_integration, plot_calibration_single, plot_calibration_combined, plot_residuals_single
from .plotting import FTIR_to_DTG, plot_FTIR, plot_TGA
from .utils import get_label, make_title

plt.style.use(STYLE)
