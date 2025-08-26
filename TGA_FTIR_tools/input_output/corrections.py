import logging
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
import scipy as sp

from ..config import EGA_NOISE
from ..plotting import plot_corr

logger = logging.getLogger(__name__)


def corr_TGA_Baseline(TGA: pd.DataFrame, baseline, plot: bool = False) -> pd.DataFrame:
    corr_data = TGA.copy()
    baseline = baseline.tga
    
    for value in ['sample_mass','heat_flow']:
        if (value in TGA.columns) and (value in baseline.columns):
            corr_data.loc[:,value] = corr_data.loc[:,value].subtract(baseline.loc[:,value])
            logger.info(f'Corrected "{value}".')

            if plot:
                plot_corr(TGA, baseline, label = value)
    
    return corr_data



def corr_FTIR(Sample, baselineData, plot: bool | Iterable | Mapping=False, co2_offs = 0):
    "corrects IR data by setting minimal adsorption to 0 "
    FTIR = Sample.ega
    # setting up output DataFrame
    corr_data = pd.DataFrame(
        index=FTIR.index,
        columns=FTIR.columns.drop(
            ["time", "sample_temp", "reference_temp"], errors="ignore"
        ),
    )
    # opens FTIR data of the baseline
    baseline = baselineData.ega
    gases = baselineData.info.gases

    # cycling through gases
    if plot == True:
        plot = gases
    elif plot == False:
        plot = []

    for gas in gases:

        # load threshold of noise from settings.ini or determine it from baseline
        if gas.lower() in EGA_NOISE:
            thresh = EGA_NOISE.getfloat(gas.lower())
        else:
            try:
                thresh = np.median(baseline[gas] - min(baseline[gas]))
            except:
                thresh = 0

        # special CO2 correction due to periodic fluctuations in signal
        if gas == "CO2":
            try:
                co2_baseline = np.array(baseline["CO2"])

                # determining the peaks and valleys in the baseline as well as their amplitude
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
                peaks, _ = sp.signal.find_peaks(
                    FTIR["CO2"],
                    height=[-tol * amplitude_baseline, tol * amplitude_baseline],
                )
                valleys, _ = sp.signal.find_peaks(
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
                    peaks, _ = sp.signal.find_peaks(
                        FTIR["CO2"]
                        - co2_baseline[x_shift + x_offs : len(FTIR) + x_shift + x_offs],
                        height=[None, None],
                        prominence=amplitude_baseline * 0.02,
                    )
                    c.append(len(peaks))
                x_offs = np.where(c == np.min(c))[0][0] - 1 + co2_offs

                co2_baseline = co2_baseline[
                    x_shift + x_offs : len(FTIR) + x_shift + x_offs
                ]
                corr_data[gas] = co2_baseline
                # in this case the corrected data has to be provided to const_baseline() to average the noise threshold
                corr_data[gas] += const_baseline(
                    FTIR[gas].subtract(co2_baseline) - min(FTIR[gas]), thresh
                ) + min(FTIR[gas].subtract(co2_baseline))

            except:
                logger.warning("Unable to align CO2 baseline with measurement.")
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
        if gas in plot:
            plot_corr(originalData=FTIR, BaselineData=corr_data, label=gas)

    return FTIR[list(gases)].subtract(corr_data)


def const_baseline(data, thres):
    "calculate y-offset so the integral of noise in data is equal to zero"
    baseline = data[data < thres]

    if len(baseline) == 0:
        return 0
    else:
        return np.sum(baseline) / len(baseline)

