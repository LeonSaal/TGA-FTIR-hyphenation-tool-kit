import re

import numpy as np
import pandas as pd

from ..config import PATHS
from .general import find_files


def read_FTIR(file_name):
    "read IR data from files"

    files = find_files(file_name, ".csv", PATHS["data"])

    # find the gases by looking at the suffix of the files
    gases = []
    paths = []
    for file in files:
        gases.append(file[file.rfind("_") + 1 : file.rfind(".")].upper())
        paths.append(file)

    if gases == []:
        return
    else:
        # make DataFrame with the first gas, keeping the time column
        data = pd.read_csv(
            files[0], delimiter=";", decimal=",", names=["time", gases[0]]
        )

        # append the IR data from the other gases as new columns
        for i in range(1, len(gases)):
            data[gases[i]] = pd.read_csv(
                files[i],
                delimiter=";",
                decimal=",",
                names=["time", gases[i]],
                usecols=[gases[i]],
            )

        # convert time column for minutes to seconds
        data["time"] = (data["time"] * 60).astype(int)
        return data.dropna()


def FTIR_info(TG_IR):
    "determine total area of detected gases as well as the molar amount (if calibrated) of gases and elements"
    info = {}

    # calculate total area of each gas
    for gas in TG_IR.info["gases"]:
        info[f"area_{gas}"] = np.sum(TG_IR.ir[gas])

    if hasattr(TG_IR, "linreg"):
        # calculate molar amount of calibrated gases
        for gas in [gas for gas in TG_IR.info["gases"] if gas in TG_IR.linreg.index]:
            info[f"mmol_{gas}"] = (
                info[f"area_{gas}"] - TG_IR.linreg["intercept"][gas]
            ) / TG_IR.linreg["slope"][gas]

        # calculate molar amount of elements in gases, assuming the elemental formaula of gases does not exceed 5 characters
        try:
            elems = list(
                set(
                    re.sub(
                        "\d",
                        "",
                        "".join([gas for gas in TG_IR.info["gases"] if len(gas) < 5]),
                    )
                )
            )
            for elem in elems:
                temp = 0
                for gas in TG_IR.linreg.index:
                    if elem in gas:
                        if (n_elem:=re.search("(?<=" + elem + ")\d", gas)) != None:
                            n = int(n_elem.group())
                        else:
                            n = 1
                        temp += n * info[f"mmol_{gas}"]
                if temp != 0:
                    info[f"mmol_{elem}"] = temp
        except:
            pass

    return info
