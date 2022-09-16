from .sample import Sample


class Baseline(Sample):
    def __init_subclass__(cls) -> None:
        return super().__init_subclass__()
