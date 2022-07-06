import os
import pandas as pd
import numpy as np

from ..config import PATHS
import logging


logger = logging.getLogger(__name__)


def samplelog(info=None, create=True, overwrite=False,**kwargs) -> pd.DataFrame:
    "load and write samplelog file with obj.info"
    path = os.path.join(PATHS["home"], "Samplelog.xlsx")

    # try to load samplelog file
    if not os.path.exists(path):
        samplelog = pd.DataFrame(columns=["alias", "reference"])
        samplelog.index.name = "name"
        if create:  # create new Samplelox.xlsx file
            samplelog.to_excel(path)
            logger.info(f"Empty 'Samplelog.xlsx' created in {path}")
    else:
        samplelog = pd.read_excel(path, index_col=0)

    # update existing samplelog file
    if info != None:
        name = info["name"]
        data = pd.DataFrame.from_dict(info, orient="index", columns=[name]).T.drop(
            ["name"], axis=1
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
            samplelog = pd.concat([samplelog, data])

        try:
            samplelog.to_excel(path)
            logger.info("Successfully updated 'Samplelog.xlsx'.")
        except:
            logger.error(
                "Unable to write on 'Samplelog.xlsx'. Please close file and try again!"
            )

    return samplelog
