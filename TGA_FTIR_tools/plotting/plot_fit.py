import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from ..config import SEP, UNITS
from ..utils import gaussian, multi_gauss
from .plotting import get_label, make_title


def plot_fit(sample, reference, title=False, y_axis="orig", **kwargs):
    if y_axis == "rel":
        ir_values = "mmol_per_mg"
    else:
        ir_values = "area"

    fit_data = sample.results["fit"][reference][["center", "height", "hwhm"]].dropna()

    for gas, params in fit_data.groupby("gas"):
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(8, 1)
        fitting = fig.add_subplot(gs[:-1, 0])
        if title:
            fitting.set_title(make_title(sample))
        error = fig.add_subplot(gs[-1, 0], sharex=fitting)
        # fitting.xaxis.set_ticks(np.arange(0, 1000, 50))

        # plotting of fit
        # x, data = sample.ir.sample_temp, sample.ir[gas]

        num_curves = params.index.size
        fitting.plot(
            sample.ir.sample_temp,
            sample.ir[gas],
            label="data",
            lw=2,
            zorder=num_curves + 1,
        )  # ,ls='',marker='x',markevery=2,c='cyan')

        x, y_data = sample.ir.sample_temp, sample.ir[gas]

        yall = multi_gauss(
            x.values, *params.height.values, *params.center.values, *params.hwhm.values
        )
        fitting.plot(x, yall, label="fit", lw=2, zorder=num_curves + 2)

        for i, ((sample_name, alias, run, group, gas), row) in enumerate(params.iterrows()):
            y = gaussian(x.values, row.height, row.center, row.hwhm)
            fitting.text(
                row.center,
                row.height,
                group,
                zorder=num_curves + 3 + i,
            )
            fitting.plot(x, y, linestyle="dashed", zorder=i)  #

        fitting.legend()
        fitting.set_xlabel(f"{get_label('sample_temp')} {SEP} ${UNITS['sample_temp']}$")
        if y_axis == "orig":
            fitting.set_ylabel(f"{get_label(gas)} {SEP} ${UNITS['ir']}$")
        elif y_axis == "rel":
            fitting.set_ylabel(
                f"{get_label(gas)} {SEP} ${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}\\,{UNITS['time']}^{{-1}}$"
            )

        # mark center on x-axis
        fitting.scatter(
            params.center,
            np.zeros(num_curves),
            marker=7,
            color="k",
            s=100,
            zorder=num_curves + 3,
        )

        # plotting of absolute difference
        abs_max = 0.05 * max(y_data)
        sqerr = sample.results["fit"][reference].loc[(sample_name, alias, run, "total", gas), "sumsqerr"]
        total = sample.results["fit"][reference].loc[(sample_name, alias, run, "total", gas), ir_values]
        error.text(
            0,
            abs_max,
            f"SQERR: {sqerr:.2e} ({100 * sqerr / total:.2f} %)",
        )  # percentage SQERR

        diff = y_data - yall
        error.plot(x, diff)
        error.hlines(0, min(x), max(x), ls="dashed")
        error.set_xlabel(f"{get_label('sample_temp')} {SEP} ${UNITS['sample_temp']}$")
        error.set_ylabel("error")
        error.set_ylim(-abs_max, abs_max)

        fitting.xaxis.set_minor_locator(
            ticker.AutoMinorLocator()
        )  # switch on minor ticks on each axis
        fitting.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        error.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        fig.savefig(f'{sample.info["name"]}_{gas}.png')
    plt.show()
