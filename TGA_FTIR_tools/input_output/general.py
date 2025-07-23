import json
import logging
import os
import re
from datetime import datetime
from pathlib import Path
from typing import List, Mapping

import pandas as pd

from ..config import COUPLING, PATH_SET, PATHS
from ..utils import download_supplementary

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
            return json.load(json_file)


def read_data(sample_name: str, profile=COUPLING["profile"]) -> pd.DataFrame:
    out = {}
    profile_specs = read_profile_json(profile)
    if not profile_specs:
        logger.error(f"Profile {profile!r} not found or empty.")
        return out

    # iterate over devices
    for key, values in profile_specs.items():
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
        if values["rename"]:
            rename = eval(values["rename"][3:]) if values.rename.startswith("fn:") else values["rename"]
            concat.rename(columns=rename, inplace=True)
        out[key] = concat

    return out
