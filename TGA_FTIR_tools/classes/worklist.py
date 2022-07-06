from dataclasses import dataclass, field
from ..classes import Sample
from typing import List
from ..fitting import fits, robustness
from ..plotting import plots, plot_robustness, bar_plot_results
import logging
import re
import pickle
from ..config import PATHS
logger = logging.getLogger(__name__)
import os


@dataclass
class Worklist:
    name: str
    samples: List[Sample] = field(default_factory=list)
    results: dict = field(default_factory=dict)

    def __add__(self, other):
        if isinstance(other, Sample):
            return Worklist([other] + self.samples)
        else:
            logger.warn("Can only add samples to Worklist")

    def __repr__(self):
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
            return Worklist(f"{self.name} [{i.start}:{i.stop}]", self.samples[i])
        elif type(i) == str:
            if i in (d := {sample.name: sample for sample in self.samples}):
                return d[i]
        elif type(i) == list:
            samples = []
            for elem in i:
                samples.append(self.__getitem__(elem))
            if len(samples) == 1:
                return samples[0]
            return Worklist(f"{self.name} {repr(i)}", samples)

    def __iter__(self):
        yield from self.samples

    def get(self, key: str, attr: str = "name"):
        return Worklist(
            key,
            [
                sample
                for sample in self.samples
                if re.search(key, sample.__dict__[attr])
            ],
        )

    def append(self, other: Sample):
        self.samples.append(other)

    def fit(self, reference: str, **kwargs):
        fits(self, reference, **kwargs)
        return self.results["fit"]

    def robustness(self, reference: str, plot=True, **kwargs):
        self.results["robustness"] = robustness(self, reference, **kwargs)
        if plot:
            self.plot('robustness')
        return self.results["robustness"]

    def plot(self, plot, **kwargs):
        if plot == "robustness" and "robustness" in self.results:
            plot_robustness(self.results["robustness"][0])
        elif plot == "results":
            bar_plot_results(self, **kwargs)
        else:
            plots(self.samples, plot, **kwargs)

    def __len__(self):
        return len(self.samples)

    def corr(self,references=None, plot=False, **kwargs):
        if type(references) == list:
            if (ls:=len(self)) != (lr:=len(references)):
                logger.error(f'Lengths of samples ({ls}) and references ({lr}) do not match.')
                return
        elif type(references)==str:
            references = [references for _ in self.samples]

        for sample, baseline in zip(self.samples, references):
            sample.corr(baseline, plot=plot, **kwargs)

    def pop(self, i):
        self.samples.pop(i)

    def save(self, fname:str=None, **kwargs):
        path_output = PATHS["output"]
        if os.path.exists(path_output) == False:
            os.makedirs(path_output)
        if not fname:
            fname= self.name
        path = os.path.join(path_output, f'{fname}.wkl')
        with open(path, "wb") as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        for sample in self.samples:
            sample.save()

    def load(self, fname:str):
        with open(os.path.join(PATHS["output"], f'{fname}.wkl'), "rb") as inp:
            obj = pickle.load(inp)
        for key in obj.__dict__:
            self.__dict__[key] = obj.__dict__[key]

