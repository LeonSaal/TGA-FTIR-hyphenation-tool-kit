import matplotlib as mpl
from .utils import get_label
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from ..config import SEP, UNITS
import pint
ureg = pint.get_application_registry()
import pandas as pd
FIGSIZE = np.array(plt.rcParams["figure.figsize"])

def plot_integration(ega_data, baselines, peaks_idx, step_starts_idx, step_ends_idx, gases, ax:mpl.axes.Axes):        
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
    for gas in gases:
        (l, ) = ax.plot(x, y[gas], label=gas)

        color = l.get_color()
        # add baseline
        for j, (step_start_idx, step_end_idx) in enumerate(zip(step_starts_idx, step_ends_idx)):
            x_baseline = (
                ega_data["sample_temp"].iloc[step_start_idx:step_end_idx]
            )
            y_baseline = ((baselines[gas][j] - ega_data[gas].min()) / (ega_data[gas].max()- ega_data[gas].min()))
            ax.plot(
                x_baseline, y_baseline, color=color, linestyle="dashed"
            )

def plot_calibration_single(x,y, linreg, ax):
    x_unit = x.dtype.units
    y_unit = y.dtype.units
    ax.scatter(x.to_numpy(), y.to_numpy())
    x_bounds = x.agg(["min", "max"]).astype(x.dtype)
    ax.plot(
        x_bounds,
        (x_bounds * linreg["slope"] + linreg["intercept"]) * y_unit,
        label="regression",
        ls="dashed",
    )
    ax.text(
        x.max(),
        y.min(),
        f'y={linreg["slope"]:.1e} x{linreg["intercept"]:+.1e}\n$R^2$={linreg["r_value"] ** 2:.3}, N={len(x)}',
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
        df = pd.DataFrame({"x":x[gas], "y":y[gas]}).dropna()
        xgas, ygas = df.x, df.y
        x_unit = xgas.dtype.units
        y_unit = ygas.dtype.units

        axdict[y_unit].scatter(xgas.to_numpy(), ygas.to_numpy(), label=f"data {get_label(gas)} (N={xgas.size})")
        xrange = xgas.agg(["min", "max"]).astype(xgas.dtype)
        axdict[y_unit].plot(
            xrange,
            (xrange * linreg["slope"][gas] + linreg["intercept"][gas]) * y_unit,
            ls="dashed",
        )
        axdict[y_unit].set_xlim(0, max(xgas) + abs(min(xgas)))
        axdict[y_unit].legend(loc=0)
    return fig, axs

def plot_residuals_single(x,y, linreg, ax):
    x_unit = x.dtype.units
    y_unit = y.dtype.units
    Y_cali = x.mul(linreg["slope"]).add(linreg["intercept"]) * y_unit
    ax.scatter(Y_cali.to_numpy(), (y - Y_cali).to_numpy(), label=f"data (N={len(x)})")
    ax.hlines(0, Y_cali.min(), Y_cali.max())
    ax.set_ylabel(f"$y_i-\\hat{{y}}_i$ {SEP} {UNITS.get('int_ega', '?')}")