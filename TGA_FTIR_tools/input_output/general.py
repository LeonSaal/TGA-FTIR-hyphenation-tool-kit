import datetime as dt
import json
import logging
import os
import re
from pathlib import Path
from typing import List, Mapping

import pandas as pd

from ..config import COUPLING, PATH_SET, PATHS

logger = logging.getLogger(__name__)


def time():
    return (
        str(dt.datetime.now().date())
        + "_"
        + str(dt.datetime.now().hour)
        + "-"
        + str(dt.datetime.now().minute).zfill(2)
        + "-"
        + str(dt.datetime.now().second).zfill(2)
        + "_"
    )


def find_path(file, parent_dir):
    "return directory of file in parent_dir and its subdirectories"
    for dirpath, _, filenames in os.walk(parent_dir):
        for filename in filenames:
            if filename[: filename.rfind(".")].find(file) != -1:
                return dirpath


def find_files(file, suffix, parent_dir):
    "return absolute path of file with suffix in parent_dir and its subdirectories"
    files = []
    for dirpath, _, filenames in os.walk(parent_dir):
        for filename in filenames:
            if (re.search(re.escape(file), filename) != None) and (
                re.search(suffix + "$", filename, flags=re.IGNORECASE)
            ):
                files.append(os.path.join(dirpath, filename))
    return files


def find_files_re(file: str, suffix: str, parent_dir: str) -> List[str]:
    files = []
    for dirpath, _, filenames in os.walk(parent_dir):
        for filename in filenames:
            if re.match(f"{re.escape(file)}{suffix}", filename, flags=re.IGNORECASE):
                filepath = os.path.join(dirpath, filename)
                files.append(filepath)
    return files


def read_profiles(profile):
    profiles = pd.read_excel(PATHS["import_profiles"], index_col=0, sheet_name=None)

    return {key: df.loc[profile, :].to_dict() for key, df in profiles.items()}


def read_profile_json(profile: str) -> Mapping:
    path = os.path.join(PATH_SET, f"{profile}.json")
    if not os.path.exists(path):
        logger.error(f"Cannot find '{profile}.json' in {PATH_SET!r}")
    else:
        return json.load(open(path, encoding="UTF-8"))


def read_data(sample_name: str, profile=COUPLING["profile"]) -> pd.DataFrame:
    out = {}
    profile = read_profile_json(profile)
    if not profile:
        return out
    
    # iterate over devices
    for key, values in profile.items():
        paths = find_files_re(sample_name, values["ext"], PATHS["data"])

        # convert strings to functions
        kwargs = values["kwargs"]
        if (arg := "converters") in kwargs:
            for col, converter in kwargs[arg].items():
                kwargs[arg][col] = eval(converter)

        frames = []
        for path in paths:
            filename = Path(path).name

            # extract suffix from filename
            pat = f'{re.escape(sample_name)}{values["ext"]}'
            if m := re.match(pat, filename, flags=re.I):
                suffix = m.group("suffix")

            # load data from path
            try:
                data = pd.read_csv(path, **kwargs)

                # rename columns
                if "map_suffix" in values:
                    data.rename(
                        {"suffix": values["map_suffix"][suffix]}, axis=1, inplace=True
                    )
                else:
                    data.rename({"suffix": suffix.upper()}, axis=1, inplace=True)

            except PermissionError:
                logger.error(f"Failed to read {key}-data from {path}")
            frames.append(data)

        # concatenate data and remove duplicate columns
        concat = pd.concat(frames, axis=1)
        concat = concat.loc[:, ~concat.columns.duplicated()]
        out[key] = concat

    return out
