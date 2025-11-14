import numpy as np
import scipy as sp

def corr_CO2(x, y, co2_offs = 0):
    "corrects IR data by setting minimal adsorption to 0 "
    # special CO2 correction due to periodic fluctuations in signal
    try:
        co2_baseline = y

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
            x,
            height=[-tol * amplitude_baseline, tol * amplitude_baseline],
        )
        valleys, _ = sp.signal.find_peaks(
            -x,
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
                x
                - co2_baseline[x_shift + x_offs : len(x) + x_shift + x_offs],
                height=[None, None],
                prominence=amplitude_baseline * 0.02,
            )
            c.append(len(peaks))
        x_offs = np.where(c == np.min(c))[0][0] - 1 + co2_offs

        co2_baseline = co2_baseline[
            x_shift + x_offs : len(x) + x_shift + x_offs
        ]
        # in this case the corrected data has to be provided to const_baseline() to average the noise threshold
        z = x-co2_baseline
        z = z-z.min()

    except:
        print("Unable to align CO2 baseline with measurement.")
        z = np.ones_like(x)*x.min()

    return z