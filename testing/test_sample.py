from TGA_FTIR_tools import Sample, SampleInfo
import pandas as pd
import inspect as ins
from typing import get_args

import pytest
import matplotlib as mpl
from .utils import corr_CO2

@pytest.mark.usefixtures("get_sample")
class TestSampleInit:
    def test_init_sample(self, get_sample):        
        # check all class attributes
        assert isinstance(get_sample.tga, pd.DataFrame)
        assert isinstance(get_sample.ega, pd.DataFrame)
        assert isinstance(get_sample.name, str)
        assert isinstance(get_sample.alias, str)
        #assert isinstance(get_sample.reference, str) or None

        # check loading of linreg before calibration
        # assert get_sample.linreg is None
        # assert get_sample.stats is None
        # assert get_sample.xcali is None
        # assert get_sample.ycali is None
        assert isinstance(get_sample.info,SampleInfo)

@pytest.mark.usefixtures("get_sample")
class TestSamplePlot:
    def test_plot_sample(self, get_sample):
        plot_params = ins.signature(get_sample.plot).parameters["plot"].default
        plot_opts = get_args(plot_params)
        for plot in plot_opts:
            if plot not in ["EGA_to_DTG", "calibration", "fit"]:
                assert isinstance(get_sample.plot(plot), mpl.axes.Axes) 
            elif plot=="fit":
                assert get_sample.plot(plot) == None
            else:
                assert get_sample.plot(plot) == True

@pytest.mark.usefixtures("get_sample", "get_baseline")
class TestSampleClassMethods():
    def test_saving(get_sample):
        get_sample.to_pickle()

        assert Sample.from_pickle(get_sample.name) == get_sample
    
    def test_corr(get_sample, get_baseline):
        if get_sample.profile != get_baseline.profile and get_sample.profile != "Otto":
            pytest.skip("Different profiles")

        get_sample.corr(get_baseline, {"ega": {"CO2":corr_CO2}})
        assert get_sample.ega != get_sample.raw.ega

    def test_fit(get_sample):
        assert isinstance(get_sample.fit("test"), pd.DataFrame)

    def test_calibrate(get_sample):
        assert get_sample.calibrate == None    

#  'dry_weight',
#  'get_value',
#  'mass_step',
#  'robustness',
#  'save',
#  'step_data'
    pass

