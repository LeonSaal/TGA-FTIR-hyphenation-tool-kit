from TGA_FTIR_tools import Worklist
import pytest
import pandas as pd

@pytest.mark.usefixtures("get_cali_worklist")
#@pytest.mark.parametrize("get_sample, get_cali_worklist", [((sample_names[0], "Netzsch"), cali_sample_names)], indirect=True)
class TestSampleCali:
        def test_cali(self, get_cali_worklist):    
            wl = get_cali_worklist
            linreg, stats, xcali, ycali= wl.calibrate(molecular_formulas = {'m:18': 'H2O', 'm:44': 'CO2'},plot=True)
            assert isinstance(stats , pd.DataFrame)
            assert isinstance(linreg , pd.DataFrame)
            assert isinstance(xcali , pd.DataFrame)
            assert isinstance(ycali , pd.DataFrame)

@pytest.mark.usefixtures("get_worklist")
class TestWorklistInit:
    def test_init_worklist(self, get_worklist):
        for worklist in get_worklist():
            print(worklist)
            assert isinstance(worklist.samples, list)
            assert isinstance(worklist.name, str)
            assert isinstance(worklist.info, pd.DataFrame)

class TestWorklistClassMethods():
    def test_saving(get_worklist):
        get_worklist.to_pickle()

        assert Worklist.from_pickle(get_worklist.name) == get_worklist