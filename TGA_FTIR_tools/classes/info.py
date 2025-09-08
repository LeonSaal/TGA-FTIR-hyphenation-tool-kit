from dataclasses import dataclass, field
from typing import Mapping, Optional
import pandas as pd
from ..config import DEFAULTS


@dataclass
class SampleInfo:
    name: Optional[str]
    reference: str = None
    alias: str = None
    profile: str = None
    corrected: bool = False
    initial_mass: float = None
    final_mass: float=None
    reference_mass_name: Optional[str] = "initial_mass"
    steps_idx: dict = field(default_factory=lambda: {"initial_mass":0, "final_mass" : -1})
    info: dict = field(default_factory=dict)

    def __post_init__(self):
        self.alias = self.name
        for key, value in self.info.items():
            self[key] = value
        del self.info

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

    def get(self, key):
        if key in self.__dict__:
            return self.__dict__[key]

    
    def to_dict(self):
        return {
            key: str(value)
            for key, value in self.__dict__.items()
            if value is not None
        }

    def to_row(self):
        return (pd.DataFrame
                .from_dict(self.to_dict(), orient="index")
                .T
                .set_index("name"))