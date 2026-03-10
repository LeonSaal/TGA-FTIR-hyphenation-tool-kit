from TGA_FTIR_tools import Sample, SampleInfo
import pandas as pd
import inspect as ins
from typing import get_args

import pytest
import matplotlib as mpl

# Sample
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

# Worklist