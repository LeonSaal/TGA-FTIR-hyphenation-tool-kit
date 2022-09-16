# loading settings from settings.ini

import configparser
import logging
import os
import shutil as sh

import requests

from .links import LINKS

fmt = "[{levelname:^7s}] {module:}.{funcName}: {message}"
logging.basicConfig(level=logging.INFO, format=fmt, style="{")

logger = logging.getLogger(__name__)


names = ["ini", "import_profiles", "fitting_params"]
config_files = ["settings.ini", "TGA_import_profiles.xlsx", "Fitting_parameter.xlsx"]
config = dict(zip(names, config_files))

PATH_DIR = os.path.realpath(os.path.dirname(__file__))
PATH_SET = os.path.join(PATH_DIR, "settings")

PATHS = {name: os.path.join(PATH_SET, file) for name, file in config.items()}

for name, path in PATHS.items():
    if not os.path.exists(path):
        url = f"{LINKS.SETTINGS}/{config[name]}"
        resp = requests.get(url)
        if resp.ok:
            dst = os.path.join(PATH_SET, config[name])
            with open(dst, "wb") as file:
                file.write(resp.content)
        else:
            logger.err(f"Unable to download '{url}'.")

if not os.path.exists(config["ini"]):
    if os.path.exists(PATHS["ini"]):
        sh.copy(PATHS["ini"], config["ini"])
    else:
        logger.error("Unable to find default settings.")

cfg = configparser.ConfigParser()
cfg.read(config["ini"])

PATHS.update(cfg["paths"])

if PATHS["home"] == "":
    cwd = os.getcwd().replace(os.sep, os.altsep)
    logger.info(f"No home path was supplied in {config['ini']}. home was set to '{cwd}'")
    cfg["paths"]["home"] = cwd
if PATHS["data"] == "" or os.path.exists(PATHS["data"]) == False:
    logger.warn(f"No valid data path was supplied in '{config['ini']}'.")
    cfg["paths"]["data"] = input("Supply directory of Data:").replace(
        os.sep, os.altsep
    )
    if os.path.exists(cfg["paths"]["data"]) == False:
        logger.error(
            "Supplied directory does not exist. Revise path in 'settings.ini' prior to continue."
        )

PATHS.update(cfg["paths"])

with open(config["ini"], "w") as configfile:
    cfg.write(configfile)

UNITS = cfg["units"]
keys = ["sample_mass", "time", "sample_temp", "molar_amount", "heat_flow", "dtg"]
units = ["mg", "min", "°C", "mmol", "mW", "mg\,min^{{-1}}"]
for key, val in zip(keys, units):
    UNITS[key] = val

SEP = UNITS["sep"]

PARAMS = cfg["parameters"]

PLOTTING = cfg["plotting"]

STYLE = cfg["plotting"]["mpl-style"]

LABELS = cfg["labels"]

COUPLING = cfg["coupling"]

SAVGOL = cfg["savgol"]

BOUNDS = cfg["fitting"]

IR_NOISE = cfg["correction"]

MERGE_CELLS = False

