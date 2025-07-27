# loading settings from settings.ini

import configparser
import logging
import os
import shutil as sh
from pathlib import Path
from typing import Literal, Union, Dict, Any
from os import PathLike
import requests

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


names = ["ini", "fitting_params"]
config_files = ["settings.ini", "Fitting_parameter.xlsx"]
config = dict(zip(names, config_files))

PATH_DIR = Path(os.path.dirname(__file__))
PATH_SET = PATH_DIR / "settings"

PATHS = {name: PATH_SET/ file for name, file in config.items()}

for name, path in PATHS.items():
    if not path.exists():
        dst = PATH_SET/ config[name]
        filename = config[name]
        download_supplementary(directory='TGA_FTIR_tools/settings',filename=filename, dst=dst)
            
if not os.path.exists(config["ini"]):
    if PATHS["ini"].exists():
        sh.copy(PATHS["ini"], config["ini"])
    else:
        logger.error("Unable to find default settings.")

cfg = configparser.ConfigParser()
cfg.read(config["ini"], encoding='ANSI')

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

UNITS = cfg["units"]
keys = ["sample_mass", "time", "sample_temp", "molar_amount", "heat_flow", "dtg"]
units = ["mg", "min", "°C", "mmol", "mW", "mg\\,min^{{-1}}"]
for key, val in zip(keys, units):
    UNITS[key] = val

SEP = UNITS["sep"]

PLOTTING = cfg["plotting"]

STYLE = cfg["plotting"]["mpl-style"]

LABELS = cfg["labels"]

COUPLING = cfg["coupling"]

SAVGOL = cfg["savgol"]

BOUNDS = cfg["fitting"]

IR_NOISE = cfg["correction"]

MERGE_CELLS = False

def fit_references(open=False):
    path = PATHS['fitting_params']
    if open:
        os.startfile(path.as_posix())
    return path

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