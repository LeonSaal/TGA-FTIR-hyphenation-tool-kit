import logging
import re

import numpy as np
import pandas as pd
import scipy as sp

from ..config import COUPLING, PATHS, SAVGOL
from .general import find_files_re, read_profile_json

logger = logging.getLogger(__name__)


WINDOW_LENGTH = int(SAVGOL.getfloat("window_length"))
POLYORDER = int(SAVGOL.getfloat("POLYORDER"))


def TGA_info(file, TGA, profile=COUPLING["profile"]):
    from ..classes import SampleInfo

    profile = read_profile_json(profile)['tga']

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
        initial_mass=tga.loc[0, "sample_mass"] if "sample_mass" in tga.columns else None,
        step_temp=[max(tga["reference_temp"])],
    )
    return info


def dry_weight(
    sample, how_dry="H2O", ref_mass="dry_mass",
):

    "determine dry point and mass from TG data"
    # make h2o to H2O
    if how_dry == "h2o":
        how_dry = how_dry.upper()

    # if how_dry is None, no dry point is determined
    if how_dry == None:
        dry_point = 1

    # if how_dry is a number, the dry point is set to that temperature
    elif type(how_dry) != str:
        # check if value is within the temperature range
        if (how_dry > max(sample.tga["sample_temp"])) or (
            how_dry <= min(sample.tga["sample_temp"])
        ):
            how_dry = None
            dry_point = 1
            logger.warning(
                "Supplied value is out of range. 'how_dry' is set to None. Be aware that 'reference_mass' is also set to 'initial_mass'."
            )
        else:
            dry_point = sample.tga["time"][sample.tga["sample_temp"] >= how_dry].values[
                0
            ]
            # check if how_dry value could not be found (very scarce, but possible..)
            if dry_point == 0:
                how_dry = None
                dry_point = 1
                logger.warning(
                    "Supplied value is out of range. 'how_dry' is set to None. Be aware that 'reference_mass' is also set to 'initial_mass'."
                )

    # if how_dry is 'H2O' or 'sample_mass', the dry point is determined from the respective data
    else:
        if how_dry == "H2O":
            try:
                ref = sample.ir.filter(items=["sample_temp", "H2O"])
            except:
                how_dry = "sample_mass"

        if how_dry == "sample_mass":
            ref = sample.tga.filter(items=["sample_temp", "sample_mass"])
            ref["sample_mass"] = -sp.signal.savgol_filter(
                sample.tga["sample_mass"] / sample.info.initial_mass,
                WINDOW_LENGTH,
                POLYORDER,
                deriv=1,
            )

        min_T = ref["sample_temp"][
            ref[how_dry]
            >= max(ref[how_dry][(ref["sample_temp"] > 50) & (ref["sample_temp"] < 200)])
        ].values[0]
        max_T = min_T + 50

        x = ref["sample_temp"][
            (ref["sample_temp"] > min_T) & (ref["sample_temp"] < max_T)
        ]
        y = ref[how_dry][(ref["sample_temp"] > min_T) & (ref["sample_temp"] < max_T)]
        slope, intercept, _, _, _ = sp.stats.linregress(x, y)

        dry_point = ref["sample_temp"][ref["sample_temp"] >= -intercept / slope].index[
            0
        ]

    # getting the dry_mass at the dry_point as well as the final weight and calculating the relative
    # mass-loss and the water content from it
    dry_temp = sample.tga["sample_temp"][dry_point]
    info = {}

    if not how_dry:
        info['reference_mass'] = "initial_mass"
        info['initial_mass'] = sample.info['initial_mass']
        times = [0] + list(
            sample.tga.index[sample.tga["sample_temp"].isin(sample.info.step_temp)]
        )
        names = sample.info.mass_steps

    elif (
        (how_dry == "H2O") or (how_dry == "sample_mass") or (type(how_dry) != str)
    ):  # turns to be true everytime
        times = [0, dry_point] + list(
            sample.tga.index[sample.tga["sample_temp"].isin(sample.info.step_temp)]
        )
        names = ["dry"] + sample.info.mass_steps
        info['reference_mass'] = "dry_mass"
        info['dry_mass'] = sample.tga["sample_mass"][dry_point]
        info['dry_temp'] = dry_temp
        info['dry_time'] = dry_point

    if ref_mass != "dry_mass":
        info['reference_mass'] = ref_mass
    weights = sample.tga["sample_mass"][sample.tga.index.isin(times)].values
    mass_loss = abs(np.diff(weights))

    info['reference'] = sample.reference

    info['final_mass'] = sample.tga["sample_mass"][len(sample.tga) - 1]
    for name, ml in zip(names, mass_loss):
        info["ML_" + name] = ml
        info["rel_ML_" + name] = ml / info[ref_mass]

    info['step_temp'] = sample.tga["sample_temp"][sample.tga.index.isin(times)].unique()
    info['step_time'] = times
    info['mass_steps'] = names
    sample.info.update(info)

