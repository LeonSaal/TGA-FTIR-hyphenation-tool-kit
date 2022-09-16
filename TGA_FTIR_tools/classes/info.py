from dataclasses import dataclass, field
from typing import Mapping, Optional

from ..config import COUPLING


@dataclass
class SampleInfo:
    name: Optional[str]
    initial_mass: Optional[float]
    reference_mass: Optional[str] = "initial_mass"
    background_delay: float = int(COUPLING["background_delay"])
    step_temp: list = field(default_factory=list)
    mass_steps: list = field(default_factory=list)

    def __post_init__(self):
        self.alias = self.name

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __iter__(self):
        yield from self.__dict__.keys()

    def update(self, other: Mapping):
        for key, value in other.items():
            self.__dict__[key] = value

    def __repr__(self) -> str:
        return (
            "{"
            + ",\n".join(
                [f"{repr(key)}:{repr(val)}" for key, val in self.__dict__.items()]
            )
            + "}"
        )
