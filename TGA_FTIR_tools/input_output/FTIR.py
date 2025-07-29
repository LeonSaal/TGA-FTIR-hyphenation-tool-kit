import numpy as np
import pandas as pd
from chempy import Substance
import logging
logger = logging.getLogger(__name__)


def FTIR_info(sample):
    "determine total area of detected gases as well as the molar amount (if calibrated) of gases and elements"
    info = {}
    gases = sample.info["gases"]
    # calculate total area of each gas
    for gas in sample.info["gases"]:
        area = sample.ega[gas].sum()
        info[f"area_{gas}"] = area if area > 0 else pd.NA

        # add molar amount for calibrated traces
        if gas in sample.linreg.index:
            if pd.isna(area):
                continue
            molar_amount = (area - sample.linreg["intercept"][gas]) / sample.linreg["slope"][gas]
            info[f"mmol_{gas}"] = molar_amount if molar_amount >= 0 else pd.NA
            

    # calculate molar amount of elements in gases, assuming the elemental formaula of gases does not exceed 5 characters
    try:
        cali_substances = sample.linreg.filter(gases, axis=0).molecular_formula
        substances = [Substance.from_formula(mf) for _, mf in cali_substances.items()]
        elem_keys = Substance.composition_keys(substances)
        for elem_key in elem_keys:
            temp = 0
            for name, mf in cali_substances.items():
                gas = Substance.from_formula(mf)
                if elem_key in Substance.composition_keys([gas]):
                    n = gas.composition[elem_key]
                    temp += n * info[f"mmol_{name}"]
            if temp != 0:
                info[f"mmol_{elem_key}"] = temp
    except Exception as e:
        logger.warning("Unable to calculate released molar masses.")
        raise e


    return info
