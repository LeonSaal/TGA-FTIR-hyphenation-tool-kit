import os
import pandas as pd
import numpy as np

from ..config import PATHS


def samplelog(info=None, create=True, overwrite=False):
    "load and write samplelog file with obj.info"
    path = os.path.join(PATHS["dir_home"], "Samplelog.xlsx")

    # try to load samplelog file
    if os.path.exists(path) == False:
        samplelog = pd.DataFrame(columns=["alias", "reference"])
        samplelog.index.name = "name"
        if create == True:  # create new Samplelox.xlsx file
            samplelog.to_excel(path)
            print("Empty 'Samplelog.xlsx' created in", path)
    else:
        samplelog = pd.read_excel(path, index_col=0)

    # update existing samplelog file
    if info != None:
        name = info["name"]
        data = pd.DataFrame.from_dict(info, orient="index", columns=[name]).T.drop(
            ["name"], 1
        )
        data.index.name = "name"

        for key in data.columns:
            if key not in samplelog.columns:
                samplelog[key] = np.nan
        if name in samplelog.index:
            if overwrite == False:
                samplelog = samplelog.fillna(data)
            else:
                samplelog.loc[[name]] = data
        else:
            samplelog = samplelog.append(data)

        try:
            samplelog.to_excel(path)
            print("> Successfully updated 'Samplelog.xlsx'.")
        except:
            print(
                "> Unable to write on 'Samplelog.xlsx'. Please close file and try again!"
            )

    return samplelog
