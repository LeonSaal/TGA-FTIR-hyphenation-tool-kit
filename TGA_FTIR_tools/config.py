# loading settings from settings.ini

import configparser
import logging
import os
import shutil as sh
from pathlib import Path
from typing import Literal, Union, Dict, Any
from os import PathLike
import requests

# set up logging
fmt = "[{levelname:^7s}] {module:}.{funcName}: {message}"
logging.basicConfig(level=logging.INFO, format=fmt, style="{")
logger = logging.getLogger(__name__)

# Define links to GitHub repository and wiki
class LINKS():
    API = 'https://api.github.com/repos/LeonSaal/TGA-FTIR-hyphenation-tool-kit/contents'
    WIKI = 'https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit/wiki'
    REPO = 'https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit#readme'

# Download supplementary files from the GitHub repository
def download_supplementary(directory: str, filename: str, dst: str|PathLike):
    base_url = f'{LINKS.API}/{directory}'
    resp = requests.get(base_url)
    if not resp.ok:
        return
    
    contents = resp.json()
    for file in contents:
        if file['name']==filename:
            download_url = file['download_url']
            if (resp:=requests.get(download_url)).ok:
                with open(dst, "wb") as file:
                    file.write(resp.content)
                    logger.info(f"Downloaded {filename!r} from repository and stored in {dst!r}.")
            break
    
    else:
        logger.error(f"Unable to download {filename!r} from repository.")

# List all import profiles available in the GitHub repository
def list_gh_profiles(device: str):
    base_url = f'{LINKS.API}/{device}'
    response = requests.get(base_url)
    profiles = {}
    if not response.ok:
        return
    
    contents = response.json()
    for file in contents:
        download_url = file.get('download_url', None)
        if not download_url:
            continue
        if (resp:=requests.get(download_url)).ok:
            profiles.update({file["name"]: resp.json()})

    return profiles

# define package path and other paths
PATH_DIR = Path(os.path.dirname(__file__))
PATH_SET = PATH_DIR / "settings"

config ={"ini":"settings.ini", "fitting_params":"Fitting_parameter.xlsx"}
PATHS = {name: PATH_SET/ file for name, file in config.items()}

# copy fresh settings into working directory         
if not os.path.exists(config["ini"]):
    if PATHS["ini"].exists():
        sh.copy(PATHS["ini"], config["ini"])
    else:
        logger.error("Unable to find default settings.")

# read settings
cfg = configparser.ConfigParser()
cfg.read(config["ini"], encoding='ANSI')

PLOTTING = cfg["plotting"]
STYLE = cfg["plotting"]["mpl-style"]
LABELS = cfg["labels"]
DEFAULTS = cfg["defaults"]
SAVGOL = cfg["savgol"]
BOUNDS = cfg["fitting"]
EGA_NOISE = cfg["correction"]
MERGE_CELLS = False
UNITS = cfg["units"]
SEP = UNITS["sep"]

# update settings and paths for data and working directory
logging.basicConfig(level=DEFAULTS["logging_level"], format=fmt, style="{")
PATHS.update({key: Path(value) for key, value in cfg["paths"].items()})
if PATHS["home"] == Path() or not PATHS["home"].exists():
    cwd = os.getcwd().replace(os.sep, os.altsep)
    logger.info(f"No valid home path was supplied in {config['ini']!r}. home was set to {cwd!r}")
    cfg["paths"]["home"] = cwd
if PATHS["data"] == Path() or not PATHS["data"].exists():
    logger.warning(f"No valid data path was supplied in {config['ini']!r}.")
    cfg["paths"]["data"] = input("Supply directory of Data:").replace(
        os.sep, os.altsep
    )
    if not os.path.exists(cfg["paths"]["data"]):
        logger.error(
            f"Supplied directory does not exist. Revise path in {config['ini']!r} before continuing"
        )
    else:
        logger.info(f"Data path set to {cfg['paths']['data']!r}.")

PATHS.update({key: Path(value) for key, value in cfg["paths"].items()})

with open(config["ini"], "w", encoding='latin-1') as configfile:
    cfg.write(configfile)

# define default units
keys = ["sample_mass", "time", "sample_temp", "molar_amount", "heat_flow", "dtg"]
units = ["mg", "min", "°C", "mmol", "mW", "mg\\,min^{{-1}}"]
for key, val in zip(keys, units):
    UNITS[key] = val

# hint location of fit parameters
def fit_references(open=False):
    path = PATHS['fitting_params']
    if open:
        os.startfile(path.as_posix())
    return path

# read, write, update settings
def read_config(cfg_file: Union[str, PathLike] = config["ini"]) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser()
    cfg.read(cfg_file, encoding='ANSI')
    return cfg

def write_config(cfg) -> configparser.ConfigParser:
    with open(config["ini"], "w", encoding='latin-1') as configfile:
        cfg.write(configfile)
    logger.info(f"Configuration written to {config['ini']!r}.")
    return cfg

def update_config(section: str, key: str, value: Union[str, int, float, bool]):
    cfg = read_config()
    if section not in cfg:
        cfg[section] = {}
    cfg[section][key] = str(value)
    write_config(cfg)
    logger.info(f"Updated {section!r} section: {key} = {value!r}.")
    return cfg[section][key]