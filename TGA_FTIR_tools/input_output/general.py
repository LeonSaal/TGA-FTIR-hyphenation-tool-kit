import json
import logging
import os
import re
from datetime import datetime
from pathlib import Path
from typing import List, Mapping

import pandas as pd
import pint_pandas

from ..config import DEFAULTS, PATH_SET, PATHS

logger = logging.getLogger(__name__)


def time():
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def find_files_re(file: str, suffix: str, parent_dir: str) -> List[str]:
    files = []
    for dirpath, _, filenames in os.walk(parent_dir):
        for filename in filenames:
            if re.match(f"{re.escape(file)}{suffix}", filename, flags=re.IGNORECASE):
                filepath = os.path.join(dirpath, filename)
                files.append(filepath)
    return files


def read_profile_json(profile: str) -> Mapping:
    path = PATH_SET/ "import_profiles"/ "profiles"
    filename = path / f"{profile}.json"

    # if not path.exists():
    #     download_supplementary(directory='import_profiles', filename=filename, dst=path)

    if not filename.exists():
        logger.error(f"Cannot find '{profile}.json' in {path!r}")
    else:
        logger.debug(f"Reading {profile}.json from {path!r}")

        with open(filename, encoding="UTF-8") as json_file:
            profile = json.load(json_file)
        for device, file in profile["data"].items():
            with open(file) as json_file:
                profile["data"][device] = json.load(json_file)
        return profile
            
def read_data(sample_name: str, profile=DEFAULTS["profile"]) -> pd.DataFrame:
    out = {}
    profile_specs = read_profile_json(profile)
    if not profile_specs:
        logger.error(f"Profile {profile!r} not found or empty.")
        return out

    # iterate over devices
    for key, values in profile_specs["data"].items():
        paths = find_files_re(sample_name, values["ext"], PATHS["data"])
        
        if not paths:
            continue

        # convert strings to functions
        kwargs = values["kwargs"]
        if (arg := "converters") in kwargs:
            for col, converter in kwargs[arg].items():
                kwargs[arg][col] = eval(converter)

        frames = []
        for path in paths:
            filename = Path(path).name

            # load data from path
            try:
                data = pd.read_csv(path, **kwargs)

            except PermissionError:
                logger.error(f"Failed to read {key}-data from {path}")

            # rename columns
            if "(?P<suffix>" in values["ext"]:
                # extract suffix from filename
                pat = f'{re.escape(sample_name)}{values["ext"]}'
                if m := re.match(pat, filename, flags=re.I):
                    suffix = m.group("suffix")
                else:
                    break

                if "map_suffix" in values:
                    data.rename(
                        {"suffix": values["map_suffix"][suffix]}, axis=1, inplace=True
                    )
                else:
                    data.rename({"suffix": suffix.upper()}, axis=1, inplace=True)

            frames.append(data)

        # concatenate data and remove duplicate columns
        concat = pd.concat(frames, axis=1)
        concat = concat.loc[:, ~concat.columns.duplicated()]

        # select specific columns if specified
        if "usecols" in values and isinstance(values["usecols"], list):
            usecols = concat.columns[values["usecols"]]
            concat = concat[usecols]
        
        # get units
        if values.get("units"):
            match values["units"]:
                # determine units from columns via regex
                case str():
                    pat = values["units"]
                    units = {col: re.match(pat, col).group("unit") for col in concat.columns if re.match(pat, col)}
                
                # pass units as list alongside columns
                case list():
                    if len(units:=values["units"]) != concat.columns.size:
                        logger.warning(f"Length of supplied units ({len(units)}) doesn't match length of columns ({concat.columns.size})!")
                        logger.warning(f"You must supply units for {concat.columns!r}.")
                        units= {}
                    else:
                        units = {col:unit for col, unit in zip(concat.columns, units)}
                    
                # map units to columns
                case dict():
                    units = {name: value for name, value in values["units"].items() if name in concat.columns}
                case _:
                    units={}

        # rename columns if specified
        if values.get("rename"):
            # handle different formats for rename
            match values["rename"]:
                case dict():
                    rename = values["rename"]
                case str():
                    rename = eval(values["rename"])
                case list():
                    if len(values["rename"]) != concat.columns.size:
                        logger.error(f"Rename list for {key} in profile {profile!r} does not match number of columns.")
                    else:
                        rename = {col: new_col if new_col else col for col, new_col in zip(concat.columns, values["rename"])}
                case _:
                    logger.error(f"Invalid rename format for {key} in profile {profile!r}.")
                    rename = {}
            oldnames = concat.columns
            concat.rename(columns=rename, inplace=True)
            
            # rename unit indices
            newnames = concat.columns
            mapper = {old:new for old, new in zip(oldnames, newnames)}
            units = {mapper[oldname]:unit for oldname, unit in units.items()}

        # assign units to columns
        concat = concat.transform({col: (lambda x, unt=unit: x.astype(f"pint[{unt}]")) if unit else (lambda y: y) for col, unit in units.items()})

        out[key] = concat

    return out
