import pandas as pd
import pytest

from TGA_FTIR_tools import Sample, Worklist, SampleInfo
from TGA_FTIR_tools.config import PATHS
from ..utils import corr_CO2, class_tester


@pytest.mark.usefixtures("get_sample")
class TestSampleInit:
    def test_init_sample(self, get_sample):
        assert isinstance(get_sample.name, str)
        assert isinstance(get_sample.alias, str)
        assert isinstance(get_sample.tga, pd.DataFrame)
        assert isinstance(get_sample.ega, pd.DataFrame)
        assert isinstance(get_sample.info, SampleInfo)
        assert "sample_mass" in get_sample.tga.columns

    def test_repr_contains_properties(self, get_sample):
        repr_text = repr(get_sample)
        assert repr_text.startswith("Sample(")
        assert f"name={get_sample.name!r}" in repr_text
        assert f"profile={get_sample.profile!r}" in repr_text

    def test_get_returns_attribute(self, get_sample):
        assert get_sample.get("name") == get_sample.name
        assert get_sample.get("alias") == get_sample.alias


@pytest.mark.usefixtures("get_sample", "get_baseline")
class TestSampleClassMethods:
    def test_to_pickle_roundtrip(self, get_sample):
        get_sample.to_pickle()
        loaded = Sample.from_pickle(get_sample.name)
        assert isinstance(loaded, Sample)
        assert loaded.name == get_sample.name
        output_file = PATHS["output"] / f"{get_sample.name}.pkl"
        if output_file.exists():
            output_file.unlink()

    def test_corr_modifies_ega(self, get_sample, get_baseline):
        if get_sample.profile != get_baseline.profile and get_sample.profile != "Otto":
            pytest.skip("Different profiles")

        original_ega = get_sample.ega.copy(deep=True)
        get_sample.corr(get_baseline, {"ega": {"CO2": corr_CO2}})
        assert not get_sample.ega.equals(original_ega)

    def test_fit_returns_dataframe(self, get_sample):
        result = get_sample.fit("test", plot=False, save=False)
        assert isinstance(result, pd.DataFrame)

    def test_calibrate_callable(self, get_sample):
        assert callable(get_sample.calibrate)

    def test_plot_invalid_option(self, get_sample):
        assert get_sample.plot("unsupported_option") is None

    def test_step_data_and_reference_mass(self, get_sample):
        step_data = get_sample.step_data()
        assert isinstance(step_data, pd.DataFrame)
        assert "sample_mass" in step_data.columns
        assert "step" in step_data.columns

        if get_sample.info.reference_mass_name in step_data.step.values:
            assert get_sample.reference_mass == step_data.loc[
                step_data.step == get_sample.info.reference_mass_name, "sample_mass"
            ].iloc[0]

    def test_get_value_returns_dataframe(self, get_sample):
        temperature = get_sample.tga.sample_temp.iloc[0]
        output = get_sample.get_value(temperature, which="sample_mass", at="sample_temp")
        assert isinstance(output, pd.DataFrame)
        assert "sample_temp" in output.columns
        assert "sample_mass" in output.columns


@pytest.mark.usefixtures("get_sample")
class TestSampleMiscellaneous:
    def test_len_and_iter(self, get_sample):
        assert len(get_sample) == 1
        assert list(get_sample) == [get_sample]

    def test_add_returns_worklist(self, get_sample):
        combined = get_sample + get_sample
        assert isinstance(combined, Worklist)
        assert len(combined) == 2

    def test_sample_contracts(self, get_sample):
        class_tester(get_sample)
