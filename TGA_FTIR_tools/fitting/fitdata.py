import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..config import MERGE_CELLS

logger = logging.getLogger(__name__)


@dataclass
class FitData:
    from ..classes import Sample

    sample: Sample
    presets: dict
    profiles: dict = field(default_factory=dict)
    peaks: pd.DataFrame = field(default=None)

    def __post_init__(self):
        self.ir = self.sample.ir
        self.info = self.sample.info
        self.linreg = self.sample.linreg
        self.gases, self.gas_links, self.links = link_groups(self.presets)
        index = pd.MultiIndex.from_tuples(
            [
                (group, gas)
                for gas in self.presets.keys()
                for group in self.presets[gas].index
            ]
            + [("total", gas) for gas in self.presets.keys()],
            names=["group", "gas"],
        )
        columns = [
            "sumsqerr",
            "center",
            "height",
            "hwhm",
            "area",
            "mmol",
            "mmol_per_mg",
        ]
        self.peaks = pd.DataFrame(columns=columns, index=index).sort_index()
        del self.sample

    def xy(self, gas):
        return self.ir.sample_temp, self.ir[gas]

    def get_bounds(self, gas, sample, predef_tol):
        param_keys = {'c': 'center', 'w':'hwhm', 'h':'height'}

        for key in ["0", "min", "max"]:
            self.presets[gas].loc[:, f"height_{key}"] = (
                self.presets[gas].loc[:, f"height_{key}"].multiply(max(sample.ir[gas]))
            )
        if gas in self.gas_links:
            df = self.links.dropna(thresh=2).replace("0", np.nan).dropna(thresh=1)
            for group in df.index:
                other = (
                    self.links.loc[group, :]
                    .index[self.links.loc[group] == "0"]
                    .values[0]
                )
                for letter in df.loc[group, gas]:
                    param = param_keys[letter]

                    if param == "height":
                        if self.linreg is not None:
                            preset = (
                                self.peaks.loc[(group, other), param]
                                / self.linreg.slope[other]
                                * self.linreg.slope[gas]
                            )
                        else:
                            continue
                    else:
                        preset = self.peaks.loc[(group, other), param]
                    self.presets[gas].loc[group, f"{param}_0"] = preset
                    self.presets[gas].loc[group, f"{param}_min"] = preset * (
                        1 - predef_tol
                    )
                    self.presets[gas].loc[group, f"{param}_max"] = preset * (
                        1 + predef_tol
                    )
        # guesses
        params_0 = np.concatenate(
            (
                [
                    self.presets[gas].loc[:, f"{key}_0"]
                    for key in ["height", "center", "hwhm"]
                ]
            )
        )

        # ...and bounds
        params_min = np.concatenate(
            (
                [
                    self.presets[gas].loc[:, f"{key}_min"]
                    for key in ["height", "center", "hwhm"]
                ]
            )
        )
        params_max = np.concatenate(
            (
                [
                    self.presets[gas].loc[:, f"{key}_max"]
                    for key in ["height", "center", "hwhm"]
                ]
            )
        )
        return params_0, params_min, params_max

    def update_peaks(self, popt, gas, xy):
        from .fitting import gaussian, multi_gauss

        x, y = xy
        num_curves = self.presets[gas].index.size
        for i in range(num_curves):
            group = self.presets[gas].index[i]
            self.peaks.loc[(group, gas), "height"] = popt[i]
            self.peaks.loc[(group, gas), "center"] = popt[i + num_curves]
            self.peaks.loc[(group, gas), "hwhm"] = popt[i + 2 * num_curves]
            self.peaks.loc[(group, gas), "area"] = np.sum(
                gaussian(x, popt[i], popt[i + num_curves], popt[i + 2 * num_curves])
            )

        profiles = pd.DataFrame()
        fit = multi_gauss(x, *popt)
        diff = y - fit
        self.peaks.loc[("total", gas), "area"] = np.sum(y)
        self.peaks.loc[("total", gas), "sumsqerr"] = np.sum(np.power(diff, 2))

        profiles["sample_temp"] = x
        profiles["data "] = y
        profiles["fit"] = fit
        profiles["diff"] = diff
        for i in range(0, num_curves):
            profiles[self.presets[gas].index[i]] = gaussian(
                x, popt[i], popt[i + num_curves], popt[i + 2 * num_curves]
            )
        self.profiles[gas] = profiles
        self.popt = popt

    def calculate_molar_mass(self):
        tot_area = (
            pd.Series(self.peaks.index.get_level_values(1))
            .apply(lambda x: self.info[f"area_{x}"])
            .values
        )
        tot_mmol = (
            pd.Series(self.peaks.index.get_level_values(1))
            .apply(lambda x: self.info[f"mmol_{x}"])
            .values
        )
        self.peaks.mmol = self.peaks.area / tot_area * tot_mmol
        self.peaks.mmol_per_mg = (
            self.peaks.mmol / self.info[self.info["reference_mass"]]
        )

    def summarize(self):
        groups = self.peaks.index.get_level_values("group").unique().to_list()
        all_groups = groups + list(
            set(
                [
                    basegroup[0]
                    for group in groups
                    if len((basegroup := group.rsplit("_", maxsplit=1))) == 2
                ]
            )
        )

        sums, means = [], []
        for group in all_groups:
            sel = self.peaks.index.get_level_values("group").str.startswith(group)
            subset = self.peaks[sel][["mmol", "mmol_per_mg"]]
            summa = subset.sum(min_count=2).dropna()
            mean = subset.mean()
            sums.append(
                pd.concat(
                    [summa.to_frame().T], keys=[(group, "sum")], names=["group", "gas"]
                )
            )
            means.append(
                pd.concat(
                    [mean.to_frame().T], keys=[(group, "mean")], names=["group", "gas"]
                )
            )

        sums = pd.concat(sums)
        means = pd.concat(means)
        self.peaks = pd.concat([self.peaks, sums, means])
        self.peaks.sort_index(level="group", inplace=True)

    def save(self, y_axis):
        f_name = f"{self.info.name}_{y_axis}.xlsx"
        f_presets = "presets.xlsx"

        with pd.ExcelWriter(f_name, engine="openpyxl") as writer, pd.ExcelWriter(
            f_presets, engine="openpyxl"
        ) as writer_presets:
            for gas, profiles in self.profiles.items():
                profiles.to_excel(writer, sheet_name=gas, merge_cells=MERGE_CELLS)
                self.presets[gas].to_excel(writer_presets, sheet_name=gas, merge_cells=MERGE_CELLS)
            self.peaks.dropna(axis=1, how="all").to_excel(writer, sheet_name="summary", merge_cells=MERGE_CELLS)


def link_groups(presets):
    data = []
    for key in presets:
        data += [presets[key].loc[:, "link"].rename(key)]
    links = pd.concat(data, axis=1)
    gas_links = links.replace("0", np.nan).dropna(thresh=1).dropna(thresh=1, axis=1)
    if (
        gas_links.dropna(axis=1).empty
        and not links.replace("0", np.nan).dropna(thresh=1).empty
    ):
        logger.warn("You cannnot predefine fitting parameters for all supplied gases!")
        gas_links = pd.DataFrame()
        links = pd.DataFrame()
    gases = [key for key in presets]
    for gas in gas_links.columns:
        gases.remove(gas)
        gases.append(gas)

    return gases, gas_links, links
