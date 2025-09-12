import json
import logging
import os
import re
from datetime import datetime
from pathlib import Path
from typing import List, Mapping

import pandas as pd

import pint_pandas
import pint
ureg = pint.get_application_registry()
Q_ = ureg.Quantity

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
    
def read_info(path:str, profile:dict):
    kwargs = profile.get("kwargs")
    info_pattern = profile.get('info_pattern')
    if kwargs and info_pattern:
        skiprows = kwargs.get("skiprows")
        skipfooter = kwargs.get("skipfooter")
    else: 
        return {}

    if not skiprows and not skipfooter:
        return {}
    
    with open(path, encoding=kwargs.get("encoding", "cp1252")) as f:
        if skipfooter:
            text = f.readlines()
        elif skiprows:
            text = [f.readline() for _ in range(skiprows)]

    text = "".join(text)
    info = {}

    for key, pat in info_pattern.items():
        req_groups = ["name", "value", "unit"]
        pat_comp = re.compile(pat)
        mapper = pat_comp.groupindex
        if len(req_groups & mapper.keys()) != 3:
            logger.warning(f"Regex {key!r} doesn't have 3 required capture groups ({req_groups!r}).")
            continue
        if matches := re.findall(pat, text):
            for match in matches:
                val = match[mapper["value"]-1]
                unit = match[mapper["unit"]-1]
                try:
                    val = pd.to_numeric(val)
                    if not pd.isna(val):
                        val = Q_(val, unit)
                except ValueError:
                    val = val

                name = match[mapper["name"]-1] if len(matches) > 1 else key
                info[name] = val
    return info
            
def read_data(sample_name: str, profile=DEFAULTS["profile"]) -> pd.DataFrame:
    out = {}
    info = {}
    init_paths = []
    profile_specs = read_profile_json(profile)
    if not profile_specs:
        logger.error(f"Profile {profile!r} not found or empty.")
        return out

    info["corrections"] = profile_specs.get("corrections", {})

    # iterate over devices
    for key, values in profile_specs["data"].items():
        paths = find_files_re(sample_name, values["ext"], PATHS["data"])
        
        if not paths:
            continue

        # convert strings to functions
        kwargs = values.get("kwargs", {})
        if (arg := "converters") in kwargs:
            for col, converter in kwargs[arg].items():
                kwargs[arg][col] = eval(converter)

        frames = []
        for path in paths:
            filename = Path(path).name

            # load data from path
            try:
                data = pd.read_csv(path, **kwargs)
                init_paths.append(path)
                info.update(read_info(path, values))

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

        # fill na
        concat = concat.interpolate(method="cubicspline")

        rename = {}

        # select specific columns if specified
        if values.get("usecols") and isinstance(values["usecols"], dict):
            indices = set()
            for val, repl in values.get("usecols").items():
                match val:
                    case int():
                        indices.add(val)
                        if repl:
                            rename[concat.columns[val]] = repl
                    case str():
                        # match slice range from string
                        if m:=re.match("^(?P<start>\\d+):(?P<stop>\\d+)$", val):
                            start, stop = int(m.group("start")), int(m.group("stop"))
                            stop = stop if stop > 0 else stop+concat.columns.size
                            indices_new = list(range(start,stop ))
                            indices.update(indices_new)
                        else:
                            # match name
                            if val in concat.columns:
                                indices.add(concat.columns.get_loc(val))
                                if repl:
                                    rename[val] = repl
                            # match regular expression to names
                            elif (bidxs:=concat.columns.str.match(val)).any():
                                indices.update({concat.columns.get_loc(col) for col, bidx in zip(concat.columns, bidxs) if bidx})
                                rename.update({col:re.sub(val, repl, col) for col, bidx in zip(concat.columns, bidxs) if bidx and repl})
                            else:
                                logger.debug(f"{val!r} is neither a valid slice, a matching regular expression nor a valid column name.")
                    case _:
                        logger.warning(f"Invalid type of element in 'usecols': {val}, {type(val)}")
            concat = concat.iloc[:,list(indices)]
        
        # get units
        if values.get("units"):
            match values["units"]:
                # determine units from columns via regex
                case str():
                    pat = values["units"]
                    units = {col:re.match(pat, col).group("unit") for col in concat.columns if re.match(pat, col)}
                
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
        else:
            units= {}

        if units:
            pass

        # rename columns if specified
        oldnames = concat.columns
        concat.rename(columns=rename, inplace=True)

        # rename unit indices
        newnames = concat.columns
        mapper = {old:new for old, new in zip(oldnames, newnames)}
        units = {mapper[oldname]: units.get(oldname) if units.get(oldname) else "dimensionless" for oldname in oldnames}

        # assign units to columns
        concat = concat.transform({col: (lambda x, unt=unit: x.astype(f"pint[{unt}]")) for col, unit in units.items()})
        out[key] = concat.pint.convert_object_dtype()
        
    info["paths"] = set(init_paths)
    out["_info"] = info

    return out
