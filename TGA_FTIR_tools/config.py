# loading settings from settings.ini

import configparser
import os
import shutil as sh

names=['ini','import_profiles','fitting_params']
cfg_files=["settings.ini", "TGA_import_profiles.xlsx", "Fitting_parameter.xlsx"]
cfg = dict(zip(names, cfg_files))

PATH_DIR = os.path.realpath(os.path.dirname(__file__))
PATH_SET = os.path.join(PATH_DIR, "settings")

PATHS ={name: os.path.join(PATH_SET, file) for name, file in cfg.items()}

if not os.path.exists(cfg['ini']):
    if os.path.exists(PATHS['ini']):
        sh.copy(PATHS['ini'], cfg['ini'])
    else:
        print('Unable to find default settings.')

config = configparser.ConfigParser()
config.read(cfg['ini'])

PATHS.update(config["paths"])

if PATHS["dir_home"] == "":
    cwd = os.getcwd().replace(os.sep, os.altsep)
    print(
        "No dir_home path was supplied in {}. dir_home was set to '{}'".format(
            cfg['ini'], cwd
        )
    )
    config["paths"]["dir_home"] = cwd
if PATHS["dir_data"] == "" or os.path.exists(PATHS["dir_data"]) == False:
    print("\nNo valid dir_data path was supplied in '{}'.".format(cfg['ini']))
    config["paths"]["dir_data"] = input("Supply directory of Data:").replace(
        os.sep, os.altsep
    )
    if os.path.exists(config["paths"]["dir_data"]) == False:
        print(
            "\n!!! Supplied directory does not exist. Revise path in 'settings.ini' prior to continue. !!!\n"
        )

with open(cfg['ini'], "w") as configfile:
    config.write(configfile)

UNITS = config["units"]
keys = ["sample_mass", "time", "sample_temp", "molar_amount", "heat_flow", "dtg"]
units = ["mg", "min", "°C", "mmol", "mW", "mg\,min^{{-1}}"]
for key, val in zip(keys, units):
    UNITS[key] = val

SEP = UNITS["sep"]

PARAMS = config["parameters"]

MOLAR_MASS = config["molar_mass"]

PLOTTING = config["plotting"]

DPI = PLOTTING.getint("dpi")

LABELS = config["labels"]

COUPLING = config["coupling"]

SAVGOL = config["savgol"]

BOUNDS = config["fitting"]

IR_NOISE = config["ir_noise"]

