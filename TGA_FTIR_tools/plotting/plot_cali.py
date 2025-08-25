import matplotlib as mpl
from .utils import get_label
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from ..config import SEP, UNITS
import pint
ureg = pint.get_application_registry()
FIGSIZE = np.array(plt.rcParams["figure.figsize"])

def plot_integration(ega_data, baselines, peaks_idx, step_starts_idx, step_ends_idx, gases, ax:mpl.axes.Axes):        
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    x = ega_data["sample_temp"]
    y = ega_data[gases].transform(lambda x: (x-x.min())/(x.max()-x.min()))
    step_starts = x[step_starts_idx]
    step_ends = x[step_ends_idx]
    peaks = x[peaks_idx]

    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.set_ylabel("relative intensity")
    ax.set_xlabel(x.pint.u)
    for step_start, peak,step_end in zip(step_starts, peaks,step_ends):
        ax.axvspan(step_start, step_end, alpha=.5)
        ax.axvline(peak, linestyle="dotted")

    # append secondary, third... y-axis on right side
    for i, gas in enumerate(gases):
        ax.plot(x, y[gas], color=colors[i], label=gas)

        # add baseline
        for j, (step_start_idx, step_end_idx) in enumerate(zip(step_starts_idx, step_ends_idx)):
            x_baseline = (
                ega_data["sample_temp"].iloc[step_start_idx:step_end_idx]
            )
            y_baseline = ((baselines[gas][j] - ega_data[gas].min()) / (ega_data[gas].max()- ega_data[gas].min()))
            ax.plot(
                x_baseline, y_baseline, color=colors[gases.index(gas)], linestyle="dashed"
            )

def plot_calibration_single(x,y, linreg, ax):
    x_unit = x.dtype.units
    y_unit = y.dtype.units
    ax.scatter(x, y)
    x_bounds = x.agg(["min", "max"]).astype(x.dtype)
    ax.plot(
        x_bounds,
        x_bounds * ureg.Quantity(linreg["slope"], y_unit / x_unit) + ureg.Quantity(linreg["intercept"], y_unit),
        label="regression",
        ls="dashed",
    )
    ax.text(
        max(x),
        min(y),
        f'y={linreg["slope"]:.1e}x{linreg["intercept"]:+.1e}\n$R^2$={linreg["r_value"] ** 2:.3}, N={len(x)}',
        horizontalalignment="right",
    )
    mf = linreg["molecular_formula"]
    label_mf = get_label(mf)
    label = f"{label_mf}" if linreg.name == mf  else f"{linreg.name}\n({label_mf})"
    ax.set_ylabel(label)
    ax.set_xlim(ureg.Quantity(0, x_unit), x.max() + x.abs().min())

def plot_calibration_combined(x,y, linreg, gases):
    y_units = set(dtype.units for dtype in y.dtypes)
    figdim = 1,len(y_units)
    fig, axs = plt.subplots(*figdim, squeeze=False, figsize = FIGSIZE * figdim[::-1])

    axdict =  {unit: ax for unit, ax in zip(y_units, axs[0])}

    for gas in gases:
        xgas = x[gas]
        ygas = y[gas]
        x_unit = xgas.dtype.units
        y_unit = ygas.dtype.units

        axdict[y_unit].scatter(xgas, ygas, label=f"data {get_label(gas)} (N = {len(x)})")
        xrange = xgas.agg(["min", "max"]).astype(xgas.dtype)
        axdict[y_unit].plot(
            xrange,
            xrange * ureg.Quantity(linreg["slope"][gas], y_unit / x_unit) + ureg.Quantity(linreg["intercept"][gas], y_unit),
            ls="dashed",
        )
        axdict[y_unit].set_xlim(0, max(xgas) + abs(min(xgas)))
        axdict[y_unit].legend(loc=0)
    return fig, axs

def plot_residuals_single(x,y, linreg, ax):
    x_unit = x.dtype.units
    y_unit = y.dtype.units
    Y_cali = x.mul(ureg.Quantity(linreg["slope"], y_unit/x_unit)).add(ureg.Quantity(linreg["intercept"], y_unit))
    ax.scatter(Y_cali, y - Y_cali, label=f"data (N = {len(x)})")
    ax.hlines(0, Y_cali.min(), Y_cali.max())
    ax.set_ylabel(f"$y_i-\\hat{{y}}_i$ {SEP} {UNITS['int_ega']}")