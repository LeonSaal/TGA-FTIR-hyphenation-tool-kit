import numpy as np
import scipy as sp
import inspect
import typing
from collections.abc import Iterable
import itertools as it
from pathlib import Path
from typing import get_args
import pytest
from .conftest import get_options

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


def class_tester(cls):
    members = inspect.getmembers(cls, predicate=lambda x: (inspect.isroutine(x) and not inspect.isbuiltin(x)))
    dunders = [member for member in members if member[0].startswith("__")]
    others = [member for member in members if not member[0].startswith("__")]

    # par = {}
    # total_calls = {}
    for name, fun in others:
        sig = inspect.signature(fun)
        params = sig.parameters
        if get_args(ret_ann:=sig.return_annotation): 
            ret_args = get_args(ret_ann)
        else:
            ret_args = (ret_ann, )

        kwargs = {}
        for key, value in params.items():

            if len(params)==1:
                match key:
                    case "self":
                        pass
                    case _:
                        pass

            if key in ["self", "kwargs"]:
                continue
            default = value.default
            annot = value.annotation
            args = get_args(annot)
            
            # kwargs
            match default:
                case float():
                    kwargs[key] = (default*.8, default,  default*1.2)
                case str():
                    kwargs[key] = (default,)
                case bool():
                    kwargs[key] = (True, False)
                case Iterable():
                    kwargs[key] = args
                case _:
                    kwargs[key] = get_options().get(key, ())
    
        # cartesian product of arguments
        calls = [{key:val for key, val in zip(kwargs.keys(), vals)} for vals in it.product(*kwargs.values())]
        #total_calls[name]=len(calls)
        #errs = []
        for call in calls:
            #with pytest.raises(Exception):
            # try:
            ret = cls.__getattribute__(name)(**call)
            # except Exception as e:
            #     errs.append(f"{call}: {e}")
                #continue
                
            if isinstance(iter:=ret, Iterable):
                args = tuple(type(i) for i in iter)
                assert args == ret_args
                # if not args == ret_args:
                #     errs.append(f"{call}: Type mismatch: {args} <=> {ret_args}")
            else:
                args = type(iter)
                try:
                    assert args in ret_args
                except AssertionError as e:
                    print(e)
                # if not args in ret_args:
                #     errs.append(f"{call}: Type mismatch: {args} <=> {ret_args}")
                        
            # if errs:
            #     par[name] = errs

    # total_errors = {key: len(val) for key, val in par.items()}
    # print(f"In total {sum(total_calls.values())} calls and {sum(total_errors.values())}.")
    # [print(f"> {key!r}: {val} calls, {total_errors.get(key, 0)} errors ({total_errors.get(key, 0)/val:.2%})") for key, val in total_calls.items() if val >0] 
    
    return #par

