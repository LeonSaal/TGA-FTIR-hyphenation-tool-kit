import numpy as np
from chempy import Substance


def FTIR_info(TG_IR):
    "determine total area of detected gases as well as the molar amount (if calibrated) of gases and elements"
    info = {}

    # calculate total area of each gas
    for gas in TG_IR.info["gases"]:
        info[f"area_{gas}"] = np.sum(TG_IR.ega[gas])

    if TG_IR.linreg is not None:
        # calculate molar amount of calibrated gases
        gases = [gas for gas in TG_IR.info["gases"] if gas in TG_IR.linreg.index]
        for gas in gases:
            info[f"mmol_{gas}"] = (
                info[f"area_{gas}"] - TG_IR.linreg["intercept"][gas]
            ) / TG_IR.linreg["slope"][gas]

        # calculate molar amount of elements in gases, assuming the elemental formaula of gases does not exceed 5 characters
        try:
            substances = [Substance.from_formula(gas) for gas in gases]
            elem_keys = Substance.composition_keys(substances)
            for elem_key in elem_keys:
                temp = 0
                for gas in TG_IR.linreg.index:
                    Gas = Substance.from_formula(gas)
                    if elem_key in Gas.composition_keys:
                        n = Gas.composition[elem_key]
                        temp += n * info[f"mmol_{gas}"]
                if temp != 0:
                    info[f"mmol_{elem_key}"] = temp
        except:
            pass


    return info
