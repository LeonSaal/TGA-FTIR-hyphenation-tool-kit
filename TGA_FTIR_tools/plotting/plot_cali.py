from .utils import get_label
import matplotlib.pyplot as plt
import numpy as np
from ..config import SEP, UNITS
import pint
ureg = pint.get_application_registry()

def plot_integration(ega_data, baselines, peaks_idx, step_starts_idx, step_ends_idx, gases, ax):        
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    graph = []
    #fig.subplots_adjust(right=0.8)
    x = ega_data["time"]
    y = ega_data[gases[0]]
    step_starts = x[step_starts_idx]
    step_ends = x[step_ends_idx]
    peaks = x[peaks_idx]
        # setup figure and plot first gas
    graph.append(ax)
    # graph[0].set_xlabel(f'{get_label("time")} {SEP} { UNITS["time"]}')
    graph[0].set_ylabel(f'{get_label(gases[0])} {y.pint.u}')
    graph[0].yaxis.label.set_color(colors[0])
    graph[0].plot(x, y)
    graph[0].set_ylim(0 - (max(y) / 20), max(y))
    for step_start, step_end in zip(step_starts, step_ends):
        graph[0].axvspan(step_start, step_end, alpha=.5)
    graph[0].vlines(peaks, 0, max(y), linestyle="dotted")

    # append secondary, third... y-axis on right side
    for i, gas in enumerate(gases[1:]):
        y = ega_data[gas]

        graph.append(graph[0].twinx())
        graph[i + 1].spines["right"].set_position(("axes", 1 + i * 0.1))
        graph[i + 1].plot(x, y, color=colors[i + 1])

        graph[i + 1].set_ylabel(f"{get_label(gas)} {y.pint.u}")
        graph[i + 1].yaxis.label.set_color(colors[i + 1])
        graph[i + 1].set_ylim(0 - (max(y) / 20), max(y))

        # add baseline
        for j, (step_start_idx, step_end_idx) in enumerate(zip(step_starts_idx, step_ends_idx)):
            x_baseline = (
                ega_data["time"].iloc[step_start_idx:step_end_idx]
            )
            graph[gases.index(gas)].plot(
                x_baseline, baselines[gas][j], color=colors[gases.index(gas)], linestyle="dashed"
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
    _, axs = plt.subplots(1,len(y_units), squeeze=False)

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

def plot_residuals_single(x,y, linreg, ax):
    x_unit = x.dtype.units
    y_unit = y.dtype.units
    Y_cali = x.mul(ureg.Quantity(linreg["slope"], y_unit/x_unit)).add(ureg.Quantity(linreg["intercept"], y_unit))
    ax.scatter(Y_cali, y - Y_cali, label=f"data (N = {len(x)})")
    ax.hlines(0, Y_cali.min(), Y_cali.max())
    ax.set_ylabel(f"$y_i-\\hat{{y}}_i$ {SEP} {UNITS['int_ega']}")