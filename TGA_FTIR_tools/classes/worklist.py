import logging
import pickle
import re
from dataclasses import dataclass, field
from typing import List, Optional

import pandas as pd

from ..classes import Sample
from ..config import PATHS
from ..fitting import fits, robustness
from ..plotting import bar_plot_results, plot_robustness, plots

logger = logging.getLogger(__name__)
import os


@dataclass
class Worklist:
    samples: List[Sample] = field(default_factory=list)
    name: Optional[str] = 'worklist'
    results: dict = field(default_factory=dict)

    def __add__(self, other):
        if isinstance(other, Sample):
            return Worklist([other] + self.samples)
        else:
            logger.warn("Can only add samples to Worklist")

    def __repr__(self) -> str:
        return (
            f"Worklist(name='{self.name}', samples= \n"
            + "\n".join(
                [f"[{i}] {sample.__repr__()}" for i, sample in enumerate(self.samples)]
            )
            + "\n)"
        )

    def __getitem__(self, i):
        if type(i) == int:
            return self.samples[i]
        elif type(i) == slice:
            return Worklist(samples = self.samples[i], name= f"{self.name} [{i.start}:{i.stop}]")
        elif type(i) == str:
            if i in (d := {sample.name: sample for sample in self.samples}):
                return d[i]
        elif type(i) == list and i is not []:
            samples = []
            for elem in i:
                samples.append(self.__getitem__(elem))
            if len(samples) == 1:
                return samples[0]
            return Worklist(samples = samples, name= f"{self.name} {repr(i)}")

    def __iter__(self):
        yield from self.samples

    def get(self, key: str, attr: str = "name"):
        return Worklist(
            samples = 
            [
                sample
                for sample in self.samples
                if re.search(key, sample.__dict__[attr])
            ],
            name = key,
        )

    def append(self, other: Sample) -> None:
        self.samples.append(other)

    def fit(self, reference: str, **kwargs) -> pd.DataFrame:
        fits(self, reference, **kwargs)
        return self.results["fit"]

    def robustness(self, reference: str, plot=True, **kwargs)-> pd.DataFrame:
        self.results["robustness"] = robustness(self, reference, **kwargs)
        if plot:
            self.plot('robustness')
        return self.results["robustness"]

    def plot(self, plot, **kwargs)-> None:
        if plot == "robustness" and "robustness" in self.results:
            plot_robustness(self.results["robustness"][0])
        elif plot == "results":
            bar_plot_results(self, **kwargs)
        else:
            plots(self.samples, plot, **kwargs)

    def __len__(self) -> int:
        return len(self.samples)

    def corr(self,references=None, plot=False, **kwargs) -> None:
        if type(references) == list:
            if (ls:=len(self)) != (lr:=len(references)):
                logger.error(f'Lengths of samples ({ls}) and references ({lr}) do not match.')
                return
        elif type(references)==str:
            references = [references for _ in self.samples]

        elif not references:
            references = [None for _ in self.samples]

        for sample, baseline in zip(self.samples, references):
            sample.corr(baseline, plot=plot, **kwargs)

    def pop(self, i: int) -> None:
        self.samples.pop(i)

    def save(self, fname:str=None, **kwargs) -> None:
        path_output = PATHS["output"]
        if not path_output.exists(path_output) :
            os.makedirs(path_output)
        if not fname:
            fname= self.name
        path = path_output/ f'{fname}.wkl'
        with open(path, "wb") as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        for sample in self.samples:
            sample.save(**kwargs)

    def load(self, fname:str) -> None:
        with open(PATHS["output"]/ f'{fname}.wkl', "rb") as inp:
            obj = pickle.load(inp)
        for key in obj.__dict__:
            self.__dict__[key] = obj.__dict__[key]

