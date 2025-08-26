from TGA_FTIR_tools import Sample, Worklist, SampleInfo
import pandas as pd
import inspect as ins
from typing import get_args
from pytest import fixture
import pytest
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import repeat
import shutil as sh
import os

cali_sample_names = [['ExpDat_250605_dd_Calciumoxalat_02_spline',
 'ExpDat_250605_dd_Calciumoxalat_03_spline',
 'ExpDat_250605_dd_Calciumoxalat_04_spline',
 'ExpDat_250605_dd_Calciumoxalat_05_spline',
 'ExpDat_250605_dd_Calciumoxalat_06_spline',
 'ExpDat_250605_dd_Calciumoxalat_07_spline'], ["dd_220722_Ca-oxalat_02", "dd_220722_Ca-oxalat_01"]]
sample_names = [cali_sample_names[0][0], "dd_220721_AS5000_01"]
profiles = ["Netzsch", "Otto"]

@fixture(params = [(samples, profile) for samples, profile in zip(cali_sample_names, profiles)], ids=profiles)
def get_cali_worklist(request):
    names, profile = request.param
    yield Worklist(names, profile=profile)

@fixture(params = [(name, profile) for name, profile in zip(sample_names, profiles)], ids=profiles)
def get_sample(request):
    name, profile = request.param
    yield Sample(name, profile=profile)

@fixture(params = ((samplename, profile) for samplename, profile in zip(sample_names*3, profiles*3)), ids=profiles*3)
def get_worklist(request):
    def init_options():
        yield Worklist(request.param[0], profile = request.param[1])
        for fun in [lambda x, y: Sample(x, profile=y), lambda x, y: [Sample(x, profile=y)]]:
            init = fun(*request.param)
            yield Worklist(init)
    return init_options

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
        assert get_sample.linreg is None
        assert get_sample.stats is None
        assert get_sample.xcali is None
        assert get_sample.ycali is None
        assert isinstance(get_sample.info,SampleInfo)

@pytest.mark.usefixtures("get_sample")
class TestSamplePlot:
    def test_plot_sample(self, get_sample):
        plot_params = ins.signature(get_sample.plot).parameters["plot"].default
        plot_opts = get_args(plot_params)
        for plot in plot_opts:
            if plot not in ["IR_to_DTG", "calibration", "fit"]:
                assert isinstance(get_sample.plot(plot), mpl.axes.Axes) 
            elif plot=="fit":
                assert get_sample.plot(plot) == None
            else:
                assert get_sample.plot(plot) == True

@pytest.mark.usefixtures("get_cali_worklist", "get_sample")
#@pytest.mark.parametrize("get_sample, get_cali_worklist", [((sample_names[0], "Netzsch"), cali_sample_names)], indirect=True)
class TestSampleCali:    
    def test_calibration(self, get_sample, get_cali_worklist):
        if get_sample.profile != get_cali_worklist.profile:
            pytest.skip("Different profiles")
        
        sample = get_sample
        wl = get_cali_worklist
        sample.calibrate(worklist=wl, mode="recalibrate", molecular_formulas = {'QMID(s:1|m:18)/A': 'H2O', 'QMID(s:2|m:44)/A': 'CO2'},plot=True)
        sample = get_sample
        assert isinstance(get_sample.stats , pd.DataFrame)
        assert isinstance(get_sample.linreg , pd.DataFrame)
        assert isinstance(get_sample.xcali , pd.DataFrame)
        assert isinstance(get_sample.ycali , pd.DataFrame)

@pytest.mark.usefixtures("get_worklist")
class TestWorklistInit:
    def test_init_worklist(self, get_worklist):
        for worklist in get_worklist():
            print(worklist)
            assert isinstance(worklist.samples, list)
            assert isinstance(worklist.name, str)

class TestSampleClassMethods():
#     def 'calibrate',
#  'check_profile',
#  'corr',
#  'dry_weight',
#  'fit',
#  'get_value',
#  'mass_step',
#  'missing_tga_columns',
#  'plot',
#  'robustness',
#  'save',
#  'step_data'
    pass

class TestWorklistClassMethods():
    pass