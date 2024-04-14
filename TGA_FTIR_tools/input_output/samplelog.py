import logging
import os

import numpy as np
import pandas as pd

from ..config import MERGE_CELLS, PATHS

logger = logging.getLogger(__name__)


def samplelog(info=None, create=True, overwrite=False,**kwargs) -> pd.DataFrame:
    "load and write samplelog file with obj.info"
    path = PATHS["home"]/ "Samplelog.xlsx"

    # try to load samplelog file
    if not path.exists():
        samplelog = pd.DataFrame(columns=["alias", "reference"])
        samplelog.index.name = "name"
        if create:  # create new Samplelox.xlsx file
            samplelog.to_excel(path, merge_cells=MERGE_CELLS)
            logger.info(f"Empty 'Samplelog.xlsx' created in {path}")
    else:
        samplelog = pd.read_excel(path, index_col=0)

    # update existing samplelog file
    if info != None:
        data = pd.DataFrame.from_dict(info)
        data.set_index("name", inplace=True)

        samplelog = pd.concat([samplelog, data])
        samplelog = samplelog[~samplelog.index.duplicated("last" if overwrite else "first")]

        try:
            samplelog.to_excel(path, merge_cells=MERGE_CELLS)
            logger.info("Successfully updated 'Samplelog.xlsx'.")
        except PermissionError:
            logger.error(
                "Unable to write on 'Samplelog.xlsx'. Please close file and try again!"
            )

    return samplelog
