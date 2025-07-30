from TGA_FTIR_tools import Sample, Worklist, SampleInfo
import pandas as pd
import inspect as ins

class TestInit:
    def test_init_sample():
        path = ""
        s = Sample(path)
        
        # check all class attributes
        assert isinstance(s.tga, pd.DataFrame)
        assert isinstance(s.ega, pd.DataFrame)
        assert isinstance(s.name, str)
        assert isinstance(s.alias, str)
        assert isinstance(s.reference, str)

        # check loading of linreg
        assert isinstance(s.linreg, pd.DataFrame)
        assert isinstance(s.stats, pd.DataFrame)
        assert isinstance(s.xcali, pd.DataFrame)
        assert isinstance(s.ycali, pd.DataFrame)
        assert isinstance(s.info, SampleInfo)


        

    def test_other_import_profile():
        pass

    def test_init_worklist():
        path = ""
        s = Worklist(path)

        assert isinstance(s.samples, list)
        assert isinstance(s.name, str)

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