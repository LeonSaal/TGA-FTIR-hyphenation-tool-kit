import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from ..config import SEP, UNITS
from ..utils import gaussian, multi_gauss
from .plotting import get_label, make_title


def plot_fit(sample, reference, title=False, y_axis="orig", **kwargs):
    if y_axis == "rel":
        ega_values = "mmol_per_mg"
    else:
        ega_values = "area"

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
        # x, data = sample.ega.sample_temp, sample.ega[gas]

        num_curves = params.index.size
        fitting.plot(
            sample.ega.sample_temp,
            sample.ega[gas],
            label="data",
            lw=2,
            zorder=num_curves + 1,
        )  # ,ls='',marker='x',markevery=2,c='cyan')

        x, y_data = sample.ega.sample_temp.to_numpy(dtype=np.float64), sample.ega[gas].to_numpy(dtype=np.float64)

        yall = multi_gauss(
            x, *params.height.values, *params.center.values, *params.hwhm.values
        )
        fitting.plot(x, yall, label="fit", lw=2, zorder=num_curves + 2)

        for i, ((ref, sample_name, alias, run, group, gas), row) in enumerate(
            params.iterrows()
        ):
            y = gaussian(x, row.height, row.center, row.hwhm)
            # fitting.text(
            #     row.center,
            #     row.height,
            #     group,
            #     zorder=num_curves + 3 + i,
            #     rotation = 45
            # )
            fitting.annotate(group, (row.center, row.height), (row.center, yall.max()*1.1),
            arrowprops=dict(facecolor='black', shrink=0.1, width=1, headwidth=3, headlength=3), rotation=45, zorder=num_curves+3)
            fitting.plot(x, y, linestyle="dashed", zorder=i)  #

        fitting.legend()
        fitting.set_xlabel(f"{get_label('sample_temp')} {SEP} ${UNITS['sample_temp']}$")
        if y_axis == "orig":
            fitting.set_ylabel(f"{get_label(gas)} {SEP} ${UNITS['int_ega']}$")
        elif y_axis == "rel":
            fitting.set_ylabel(
                f"{get_label(gas)} {SEP} ${UNITS['molar_amount']}\\,{UNITS['sample_mass']}^{{-1}}\\,{UNITS['time']}^{{-1}}$"
            )

        # mark center on x-axis
        # fitting.scatter(
        #     params.center,
        #     np.zeros(num_curves),
        #     marker=7,
        #     color="k",
        #     s=100,
        #     zorder=num_curves + 3,
        # )

        # fitting.vlines(
        #     params.center,
        #     np.zeros(num_curves),
        #     params.height,
        #     color="k",
        #     zorder=num_curves + 3,
        # )

        # plotting of absolute difference
        abs_max = 0.05 * max(y_data)
        sqerr = sample.results["fit"][reference].loc[
            (ref,sample_name, alias, run, "total", gas), "sumsqerr"
        ]
        total = sample.results["fit"][reference].loc[
            (ref,sample_name, alias, run, "total", gas), ega_values
        ]
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

