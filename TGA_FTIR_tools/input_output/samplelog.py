import logging
import os

import numpy as np
import pandas as pd

from ..config import MERGE_CELLS, PATHS

logger = logging.getLogger(__name__)


def samplelog(data=None, create=True, overwrite=False, sheet_name=0,**kwargs) -> pd.DataFrame:
    "load and write samplelog file with obj.info"
    path = PATHS["home"]/ "Samplelog.xlsx"

    # try to load samplelog file
    if not path.exists():
        samplelog = pd.DataFrame(columns=["alias", "sample", "run", "reference", "profile"])
        samplelog.index.name = "name"
        if create:  # create new Samplelox.xlsx file
            samplelog.to_excel(path, merge_cells=MERGE_CELLS)
            logger.info(f"Empty 'Samplelog.xlsx' created in {path}")
    else:
        samplelog = pd.read_excel(path, index_col=0, sheet_name=sheet_name)

    # update existing samplelog file
    if isinstance(data, pd.DataFrame):
        # convert data to string
        samplelog = pd.concat([samplelog, data])
        samplelog = samplelog[~samplelog.index.duplicated(keep="last" if overwrite else "first")]

        try:
            with pd.ExcelWriter(path, if_sheet_exists="replace", mode="a") as writer:
                samplelog.to_excel(writer, merge_cells=MERGE_CELLS, sheet_name=sheet_name)
            logger.info("Successfully updated 'Samplelog.xlsx'.")
        except PermissionError:
            logger.error(
                "Unable to write on 'Samplelog.xlsx'. Please close file and try again!"
            )

    return samplelog
