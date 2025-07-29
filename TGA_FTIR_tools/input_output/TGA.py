import logging
import re

import numpy as np
import pandas as pd
import scipy as sp

from ..config import DEFAULTS, PATHS, SAVGOL
from .general import find_files_re, read_profile_json

logger = logging.getLogger(__name__)


WINDOW_LENGTH = int(SAVGOL.getfloat("window_length"))
POLYORDER = int(SAVGOL.getfloat("POLYORDER"))


def TGA_info(file, TGA, profile=DEFAULTS["profile"]):
    from ..classes import SampleInfo

    profile = read_profile_json(profile)["data"]['tga']

    "extract TG info e.g. measurement time, initial mass... from TG file"
    # open file from TGA in given directory and make a DataFrame from it
    path = find_files_re(file, profile["ext"], PATHS["data"])[0] 
    with open(path, encoding=profile.get("kwargs").get("encoding", "utf-8")) as f:
        text = f.readlines()

    text = "".join(text)
    info = SampleInfo(file)

    for key, pat in profile['info_pattern'].items():
        if m := re.findall(pat, text):
            for (k, val) in m:
                try:
                    val = pd.to_numeric(val)
                except ValueError:
                    val = val

                if len(m)> 1:
                    info[k] = val
                else:
                    info[key]=val
            
    if not info["initial_mass"] and "sample_mass" in TGA.columns:
        info["initial_mass"]=TGA["sample_mass"].iloc[0]

    return info


def default_info(name, tga):
    from ..classes import SampleInfo
    info = SampleInfo(
        name=name,
        initial_mass=tga["sample_mass"].iloc[0] if "sample_mass" in tga.columns else None,
        final_mass = tga["sample_mass"].iloc[-1] if "sample_mass" in tga.columns else None,
        steps_idx = {"initial":0, "final":tga.index.size-1}
    )
    return info


def dry_weight(sample, how_dry="H2O"):
    "determine dry point and mass from TG data"

    # if how_dry is None, no dry point is determined
    if how_dry == None:
        dry_point_idx = 0

    if how_dry=="H2O" and sample.ega is None:
        logger.error(f"If {how_dry=!r}, EGA data is required.")

    if col:="sample_temp" not in sample.ega.columns:
        logger.error(f"{col!r} must be present in EGA data.")
        return

    min_temp, max_temp = sample.tga["sample_temp"].min() , sample.tga["sample_temp"].max()
    
    match how_dry:
        # if how_dry is a number, the dry point is set to that temperature
        case float() | int():
            # check if value is within the temperature range
            if (min_temp < how_dry < max_temp):
                dry_point_idx = sample.tga["time"][sample.tga["sample_temp"] >= how_dry].argmax()
            else:
                logger.error(f"Supplied value is out of range. Must be within {min_temp:.2f} to {max_temp:.2f}")
                return 
                # check if how_dry value could not be found (very scarce, but possible..)

        # if how_dry is 'H2O' or 'sample_mass', the dry point is determined from the respective data
        case str():
            if how_dry == "H2O":
                ref = sample.ega.copy()#.filter(items=["sample_temp", "H2O"])
            elif how_dry == "sample_mass":
                ref = sample.tga.copy()#.filter(items=["sample_temp", "sample_mass"])
                ref["sample_mass"] = -sp.signal.savgol_filter(
                    sample.tga["sample_mass"] / sample.info.initial_mass,
                    WINDOW_LENGTH,
                    POLYORDER,
                    deriv=1,
                )
            else:
                logger.info(f"'H2O' not in EGA-data. Defaulting to '{how_dry}")
            
            # look for signal peak between 50 and 200 °C
            peak_signal_idx = ref[how_dry][(50 < ref["sample_temp"]) & (ref["sample_temp"] < 200)].argmax()
            min_T = ref["sample_temp"].iloc[peak_signal_idx]
            max_T = min_T + 50
            range_T = (min_T < ref["sample_temp"]) & (ref["sample_temp"]< max_T)

            # fit line to slope of peak and find intersection with temperature signal
            x = ref["sample_temp"][range_T]
            y = ref[how_dry][range_T]
            slope, intercept, _, _, _ = sp.stats.linregress(x, y)

            intersection = (ref["sample_temp"] >= -intercept / slope)
            if intersection.any():
                dry_point_idx = ref["sample_temp"][intersection].argmax()
            else:
                logger.error(f"Unable to determine dry point from {how_dry}-trace.")
                return
        
        case _:
            logger.error(f"{how_dry=!r} is of the wrong type. Must be either a number between {min_temp:.2f} to {max_temp:.2f} or one of ['H2O' and 'sample_mass'].")
            pass

    # getting the dry_mass at the dry_point as well as the final weight and calculating the relative
    # mass-loss and the water content from it
    dry_temp = sample.tga["sample_temp"][dry_point_idx]
    info = {}

    info["steps_idx"] = {"dry":dry_point_idx}
    info['reference_mass_name'] = "dry_mass"
    info['dry_mass'] = sample.tga["sample_mass"][dry_point_idx]
    info['dry_temp'] = dry_temp
    sample.info.update(info)

