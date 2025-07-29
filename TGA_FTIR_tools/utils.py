import os
from typing import Literal
from .config import PATH_SET
import numpy as np
from numba import jit
import json
import logging
from os import PathLike
from pathlib import Path
from functools import lru_cache
from tabulate import tabulate
import pandas as pd
from .config import list_gh_profiles, download_supplementary
from molmass import Formula, ELEMENTS

logger = logging.getLogger(__name__)
ln2 = np.log(2)

@jit
def gaussian(x, height, center, hwhm):
    "evaluate gaussian function with height, center and HWHM at x"
    return height * np.exp(-ln2 * np.power((x - center) / hwhm, 2))


def multi_gauss(x, *args):
    "evaluate sum of multiple gaussian functions with height, center and HWHM at x"
    n = int(len(args) / 3)
    heights = args[:n]
    centers = args[n : 2 * n]
    hwhms = args[2 * n : len(args)]
    s = 0
    for i in range(n):
        s = s + gaussian(x, heights[i], centers[i], hwhms[i])
    return s


# Read a JSON file and cache the result for performance
@lru_cache
def read_json(file: str|PathLike) -> dict:
    with open(file) as f:
        return json.load(f)

def check_profile_exists(profile: str, directory: Path=PATH_SET / "import_profiles") -> bool:
    """
    Check if a profile exists in the specified directory.
    
    Args:
        profile (str): The name of the profile to check.
        directory (Path): The directory where profiles are stored.
        
    Returns:
        bool: True if the profile exists, False otherwise.
    """
    path = directory / "profiles" / f"{profile}.json"

    if not path.exists():
        logger.error(f"Profile '{profile}' does not exist in {path.parent.as_posix()!r}")
        return False
    return path.exists() and path.is_file()

# Select an import profile from the local or remote directory
def select_import_profile(directory: Path=PATH_SET / "import_profiles", loc: Literal["local", "remote"]="local"):
    profiledir = directory / "profiles"

    if not profiledir.exists():
        profiledir.mkdir()
    
    if (files:=os.listdir(profiledir)):
        filelist = [["", "profile", "file"]]+[[i,file.removesuffix(".json"), file] for i, file in enumerate(files) if file.endswith(".json")]
        print(tabulate(filelist, headers="firstrow", tablefmt="pretty"))
        while True:
            selection = input("Select profile by number or 'n' to create new profile:")
            if selection in [str(i) for i in range(len(files))]:
                profile = files[int(selection)]
                break
            elif selection == "n":
                profile = create_import_profile(directory, loc=loc)
                break

    else:
        logger.info("No profiles found.")
        profile = create_import_profile(directory, loc=loc)
    
    if not profile:
        profile = select_import_profile(directory, loc=loc)

    return profile.removesuffix(".json")

# Create a new import profile by selecting devices and their respective files
def create_import_profile(directory: Path=PATH_SET / "import_profiles", loc: Literal["local", "remote"]="local"):
    data = {}
    logger.info(f"Creating new import profile in {directory.as_posix()!r}.")

    # iterate over devices
    for device in ["tga", "ega"]:
        folder = directory / device
        if not folder.exists():
            continue
        options = []
        if loc == "local":
            for file in folder.iterdir():
                if not file.suffix ==".json":
                    continue
                try:
                    contents = read_json(file)["spec"]
                except KeyError as e:
                    logger.error(f"File {file.name!r} does not contain a valid profile. Skipping.")
                    continue
                contents["file"] = file.name
                options+= [contents]
        elif loc == "remote":
            profiles = list_gh_profiles("import_profiles")
            options = [{"spec" :spec[device]["spec"], "file": fname} for fname, spec in profiles.items()]

        df = pd.DataFrame.from_dict(options)
        print(f"Choose {device.upper()}:")
        print(tabulate(df, headers="keys", tablefmt="pretty", showindex=True))
        while True:
            inp = input("choose number or 's' to skip or 'n' to make new profile")
            if inp in [str(i) for i in range(df.index.size)]:
                file = df.file[int(inp)]
                if not (folder / file).exists():
                    download_supplementary(f"import_profiles/{device}", filename=file, dst=f"import_profiles/{file}")
                data[device] = (folder / file).as_posix()
                break
            elif inp == "n":
                raise NotImplementedError
            elif inp == "s":
                break

    if not data:
        logger.error("No profile selected. Aborting.")
        return  
    
    # make profile from parts
    profile = {"data": data}
    fields = ["name_pattern", "mass_resolution_ug"]
    profile["supplementary"] = {field: "" for field in fields}
    # save profile
    while True:
        name = input("Enter profile name!")
        if name:
            fname = name + ".json"
            path = directory / "profiles" 
            if not path.exists():
                path.mkdir() 
            with open(path/ fname, "w") as f:
                json.dump(profile, f, indent = 4)
            logger.info(f"Saved profile {name!r} to {path.as_posix()!r}")
            break
    return fname

def validate_mf(mf:str) -> bool:
    symbols = [e.symbol for e in ELEMENTS]
    pattern = f"(({'|'.join(symbols)})\\d*)+"
    return re.match(pattern, mf) is not None