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
        samplelog = pd.DataFrame(columns=["alias", "sample", "run", "reference", "profile"]).pint.quantify()
        samplelog.index.name = "name"
        if create:  # create new Samplelox.xlsx file
            samplelog.pint.dequantify().rename({"":"No Unit"}, level=1, axis=1).to_excel(path, merge_cells=MERGE_CELLS)
            logger.info(f"Empty 'Samplelog.xlsx' created in {path}")
    else:
        with pd.ExcelFile(path) as ef:
            sheet_names = ef.sheet_names
            if isinstance(sheet_name, int):
                if 0  <= sheet_name < len(sheet_names):
                    sheet_name = sheet_names[sheet_name]
                else:
                    logger.warning(f"{sheet_name} is an invalid sheet index. {path.name!r} has {len(sheet_names)} sheet(s) (zero-indexed).")
                    return
            if sheet_name in sheet_names:
                samplelog = pd.read_excel(ef, index_col=0, header=[0,1], sheet_name=sheet_name).pint.quantify()
            else:
                samplelog = pd.DataFrame(columns=["alias", "sample", "run", "reference", "profile"]).pint.quantify()

    # update existing samplelog file
    if isinstance(data, pd.DataFrame):
        samplelog = data
        try:
            with pd.ExcelWriter(path, if_sheet_exists="replace", mode="a", engine="openpyxl") as writer:
                sheet_name = f"Sheet{sheet_name + 1}" if isinstance(sheet_name, int) else sheet_name
                samplelog.pint.dequantify().rename({"":"No Unit"}, level=1, axis=1).to_excel(writer, merge_cells=MERGE_CELLS, sheet_name=sheet_name)
            logger.info("Successfully updated 'Samplelog.xlsx'.")
        except PermissionError:
            logger.error(
                "Unable to write on 'Samplelog.xlsx'. Please close file and try again!"
            )

    return samplelog
