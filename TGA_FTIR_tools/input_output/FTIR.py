import numpy as np
import pandas as pd
from molmass import Formula
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
            molar_amount = (area.magnitude - sample.linreg["intercept"][gas]) / sample.linreg["slope"][gas]
            info[f"mmol_{gas}"] = molar_amount if molar_amount >= 0 else 0
            

    # calculate molar amount of elements in gases, assuming the elemental formaula of gases does not exceed 5 characters
    try:
        cali_substances = sample.linreg.filter(gases, axis=0).molecular_formula
        substances = [Formula(mf) for _, mf in cali_substances.items()]
        elems = {f.composition().keys() for f in substances}
        for elem in elems:
            temp = 0
            for name, mf in cali_substances.items():
                gas = Formula(mf)
                if elem in gas.composition().keys():
                    n = gas.composition[elem].count
                    temp += n * info[f"mmol_{name}"]
            if temp != 0:
                info[f"mmol_{elem}"] = temp
    except Exception as e:
        logger.warning("Unable to calculate released molar masses.")
        raise e


    return info
