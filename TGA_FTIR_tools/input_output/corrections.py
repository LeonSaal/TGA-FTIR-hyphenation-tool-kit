import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from .general import find_files
from ..config import PATHS, PARAMS, UNITS, SEP, IR_NOISE
from .FTIR import read_FTIR
from ..plotting import get_label


def corr_TGA(TGA, file_baseline, plot=False):
    "corrects TG data by subtracting buoyancy blank value of crucible"
    corr_data = TGA.copy()
    path_baseline = find_files(file_baseline, ".txt", PATHS["dir_data"])[0]

    # opens the buoyancy blank value 'baseline' and substracts them from the original data
    try:
        reference_mass = pd.read_csv(
            path_baseline,
            delim_whitespace=True,
            decimal=",",
            names=["Index", "time", "sample_temp", "reference_temp", "sample_mass"],
            skiprows=13,
            skipfooter=11,
            converters={"sample_mass": lambda x: float(x.replace(",", "."))},
            engine="python",
            encoding="latin-1",
        ).drop(columns="Index")
        corr_data["sample_mass"] = corr_data["sample_mass"].subtract(
            reference_mass["sample_mass"]
        )
    except:
        print(">", path_baseline, " was not found.")
        return None
    try:
        path_mW = find_files(file_baseline, "_mW.txt", PATHS["dir_data"])[0]
        reference_heat_flow = pd.read_csv(
            path_mW,
            delim_whitespace=True,
            decimal=",",
            names=["Index", "time", "sample_temp", "reference_temp", "heat_flow"],
            skiprows=13,
            skipfooter=11,
            converters={"sample_mass": lambda x: float(x.replace(",", "."))},
            usecols=["heat_flow"],
            engine="python",
            encoding="latin-1",
        )
        corr_data["heat_flow"] = corr_data["heat_flow"].subtract(
            reference_heat_flow["heat_flow"]
        )
    except:
        pass

    # plotting of data, baseline and corrected value
    if plot == True:
        fig, ax = plt.subplots()
        x = TGA["sample_temp"]
        y = TGA["sample_mass"]
        ax.plot(x, y, label="data")
        ax.plot(x, reference_mass["sample_mass"][: len(TGA)], label="baseline")
        ax.plot(x, corr_data["sample_mass"], label="corrected")
        ax.set_xlabel(
            "{} {} {}".format(PARAMS["sample_temp"], SEP, UNITS["sample_temp"])
        )
        ax.set_ylabel(
            "{} {} {}".format(PARAMS["sample_mass"], SEP, UNITS["sample_mass"])
        )
        ax.legend()
        ax.xaxis.set_minor_locator(
            ticker.AutoMinorLocator()
        )  # switch on minor ticks on each axis
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.set(title="TGA baseline correction")
        plt.show()

    return corr_data


def corr_FTIR(Sample, file_baseline, plot=False):
    "corrects IR data by setting minimal adsorption to 0 "
    FTIR = Sample.ir
    # setting up output DataFrame
    corr_data = pd.DataFrame(
        index=FTIR.index,
        columns=FTIR.columns.drop(
            ["time", "sample_temp", "reference_temp"], errors="ignore"
        ),
    )
    # opens FTIR data of the baseline
    try:
        baseline = read_FTIR(file_baseline)
        gases = set(baseline.columns.drop(["time"]).values).intersection(
            Sample.info["gases"]
        )
        print("Baseline found for {}".format(", ".join(gases)))
    except:
        print("No baseline data found.")
        baseline = pd.DataFrame(
            index=FTIR.index,
            columns=FTIR.columns.drop(
                ["time", "sample_temp", "reference_temp"], errors="ignore"
            ),
        )

    # cycling through gases
    for gas in gases:

        # load threshold of noise from settings.ini or determine it from baseline
        if gas.lower() in IR_NOISE:
            thresh = IR_NOISE.getfloat(gas.lower())
        else:
            try:
                thresh = np.median(baseline[gas] - min(baseline[gas]))
            except:
                thresh = 0

        # special CO2 correction due to periodic fluctuations in signal
        if gas == "CO2":
            try:
                co2_baseline = np.array(baseline["CO2"])

                # determining the peaks and valleys in the baseline as well as its amplitude
                peaks_baseline, properties_baseline = sp.signal.find_peaks(
                    co2_baseline, height=[None, None]
                )
                valleys_baseline, valley_properties_baseline = sp.signal.find_peaks(
                    -co2_baseline, height=[None, None]
                )
                amplitude_baseline = np.mean(
                    properties_baseline["peak_heights"]
                ) + np.mean(valley_properties_baseline["peak_heights"])

                # in the original data the peaks and valleys that have similar height as the baseline are determined
                tol = 1.5
                peaks, properties = sp.signal.find_peaks(
                    FTIR["CO2"],
                    height=[-tol * amplitude_baseline, tol * amplitude_baseline],
                )
                valleys, valley_properties = sp.signal.find_peaks(
                    -FTIR["CO2"],
                    height=[None, None],
                    prominence=amplitude_baseline * 0.05,
                )

                # the median distance between between baseline-peaks, aka the period is determined
                dist_peaks = np.diff(peaks_baseline)
                len_period = int(np.median(dist_peaks))

                # determination of the phase shift in x direction by checking if there is also a valley in the baseline in proximity.
                # the necessary x shift is calculated as the mode of the differences
                dists = []
                j = 0
                for valley in valleys:
                    while j < len(valleys_baseline) - 1 and (
                        valleys_baseline[j] - valley <= 0
                    ):
                        j = j + 1
                    if valleys_baseline[j] - valley >= 0:
                        dists.append(valleys_baseline[j] - valley)

                x_shift = int(sp.stats.mode(dists)[0])

                # elongating the baseline by one period
                period = co2_baseline[:len_period]
                co2_baseline = np.concatenate((period, co2_baseline), axis=None)

                # shifting the baseline in x direction
                c = []
                for x_offs in range(-1, len_period % x_shift + 1):
                    peaks, props = sp.signal.find_peaks(
                        FTIR["CO2"]
                        - co2_baseline[x_shift + x_offs : len(FTIR) + x_shift + x_offs],
                        height=[None, None],
                        prominence=amplitude_baseline * 0.02,
                    )
                    c.append(len(peaks))
                x_offs = np.where(c == np.min(c))[0][0] - 1

                co2_baseline = co2_baseline[
                    x_shift + x_offs : len(FTIR) + x_shift + x_offs
                ]
                corr_data[gas] = co2_baseline
                # in this case the corrected data has to be provided to const_baseline() to average the noise threshold
                corr_data[gas] += const_baseline(
                    FTIR[gas].subtract(co2_baseline) - min(FTIR[gas]), thresh
                ) + min(FTIR[gas].subtract(co2_baseline))

            except:
                print("WARNING: Unable to align CO2 baseline with measurement.")
                corr_data[gas] = np.zeros(len(FTIR))
                corr_data[gas] += const_baseline(
                    FTIR[gas] - min(FTIR[gas]), thresh
                ) + min(FTIR[gas])
        else:
            corr_data[gas] = np.zeros(len(FTIR))
            corr_data[gas] += const_baseline(FTIR[gas] - min(FTIR[gas]), thresh) + min(
                FTIR[gas]
            )

        # plotting of baseline, data and the corrected data
        if plot:
            try:
                x = FTIR["sample_temp"]
            except:
                x = FTIR["time"]
                x /= 60

            fig, ax = plt.subplots()
            ax.plot(x, FTIR[gas], label="data")
            ax.plot(x, corr_data[gas], label="baseline")
            ax.plot(x, FTIR[gas].subtract(corr_data[gas]), label="corr. data")

            ax.legend()
            if x.name == "time":
                ax.set_xlabel("{} {} {}".format(PARAMS["time"], SEP, UNITS["time"]))
            elif x.name == "sample_temp":
                ax.set_xlabel(
                    "{} {} {}".format(PARAMS["sample_temp"], SEP, UNITS["sample_temp"])
                )
            ax.set_ylabel("{} {} {}".format(get_label(gas), SEP, UNITS["IR"]))
            ax.xaxis.set_minor_locator(
                ticker.AutoMinorLocator()
            )  # switch on minor ticks on each axis
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            ax.set(title="{} baseline correction".format(get_label(gas)))
            plt.show()

    return FTIR[list(gases)].subtract(corr_data)


def const_baseline(data, thres):
    "calculate y-offset so the integral of noise in data is equal to zero"
    baseline = data[data < thres]

    if len(baseline) == 0:
        return 0
    else:
        return np.sum(baseline) / len(baseline)

