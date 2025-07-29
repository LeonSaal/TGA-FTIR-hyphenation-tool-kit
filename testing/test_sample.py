from TGA_FTIR_tools import Sample, Worklist
import pandas as pd

class TestInit:
    def test_init_sample():
        path = ""
        s = Sample(path)
        
        # check all class attributes
        assert isinstance(s.tga, pd.DataFrame)
        assert isinstance(s.ega, pd.DataFrame)

    def test_init_worklist():
        path = ""
        s = Worklist(path)

        assert isinstance(s.samples, list)
        assert isinstance(s.name, str)