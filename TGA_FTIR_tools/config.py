# loading settings from settings.ini

import configparser
import logging
import os
import shutil as sh
from pathlib import Path

from .links import download_supplementary

fmt = "[{levelname:^7s}] {module:}.{funcName}: {message}"
logging.basicConfig(level=logging.INFO, format=fmt, style="{")

logger = logging.getLogger(__name__)



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
    logger.warn(f"No valid data path was supplied in {config['ini']!r}.")
    cfg["paths"]["data"] = input("Supply directory of Data:").replace(
        os.sep, os.altsep
    )
    if not os.path.exists(cfg["paths"]["data"]):
        logger.error(
            f"Supplied directory does not exist. Revise path in {config['ini']!r} prior to continue."
        )

PATHS.update({key: Path(value) for key, value in cfg["paths"].items()})

with open(config["ini"], "w", encoding='latin-1') as configfile:
    cfg.write(configfile)

UNITS = cfg["units"]
keys = ["sample_mass", "time", "sample_temp", "molar_amount", "heat_flow", "dtg"]
units = ["mg", "min", "°C", "mmol", "mW", "mg\,min^{{-1}}"]
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
        os.startfile(path)
    return path
