from copy import deepcopy

import pandas as pd
import pytest

from TGA_FTIR_tools import Worklist, Sample
from TGA_FTIR_tools.config import PATHS


@pytest.mark.usefixtures("get_cali_worklist")
class TestWorklistCalibration:
    def test_calibrate_returns_dataframes(self, get_cali_worklist):
        wl = get_cali_worklist
        linreg, stats, xcali, ycali = wl.calibrate(
            molecular_formulas={"m:18": "H2O", "m:44": "CO2"},
            plot=False,
        )
        assert isinstance(linreg, pd.DataFrame)
        assert isinstance(stats, pd.DataFrame)
        assert isinstance(xcali, pd.DataFrame)
        assert isinstance(ycali, pd.DataFrame)


@pytest.mark.usefixtures("get_worklist")
class TestWorklistInitialization:
    def test_init_worklist(self, get_worklist):
        for worklist in get_worklist():
            assert isinstance(worklist.name, str)
            assert isinstance(worklist.names, list)
            assert isinstance(worklist.info, pd.DataFrame)
            assert all(isinstance(sample, Sample) for sample in worklist.names)
            assert len(worklist) == len(worklist.names)
            assert worklist[0] == worklist.names[0]

            slice_result = worklist[0:1]
            assert isinstance(slice_result, Worklist)
            assert len(slice_result) == 1

            if len(worklist.names) > 1:
                one_item = worklist[[0]]
                assert isinstance(one_item, Sample)

            found = worklist[worklist.names[0].name]
            assert isinstance(found, Sample)

    def test_repr_and_get(self, get_worklist):
        for worklist in get_worklist():
            repr_text = repr(worklist)
            assert worklist.name in repr_text
            assert "[0]" in repr_text

            pattern = worklist[0].name.split("_")[0]
            filtered = worklist.get(pattern)
            assert isinstance(filtered, Worklist)
            assert filtered.names

    def test_append_and_add(self, get_worklist):
        for worklist in get_worklist():
            sample = worklist[0]
            extension = worklist + sample
            assert isinstance(extension, Worklist)
            assert len(extension) == len(worklist.names) + 1

            original_length = len(worklist)
            worklist.append(sample)
            assert len(worklist) == original_length + 1
            assert worklist.names[-1] == sample

    def test_pop_and_iter_and_len(self, get_worklist):
        for worklist in get_worklist():
            original = len(worklist)
            popped = worklist.pop(0)
            assert isinstance(popped, Sample)
            assert len(worklist) == original - 1
            assert list(iter(worklist)) == worklist.names

    def test_to_pickle_roundtrip(self, get_worklist):
        for worklist in get_worklist():
            worklist.to_pickle()
            loaded = Worklist.from_pickle(worklist.name)
            assert isinstance(loaded, Worklist)
            assert loaded.name == worklist.name
            output_file = PATHS["output"] / f"{worklist.name}.pkl"
            if output_file.exists():
                output_file.unlink()

    def test_plot_fit_without_results_returns_none(self, get_worklist):
        for worklist in get_worklist():
            assert worklist.plot("fit", save=False) is None

    def test_results_property_handles_empty_and_nonempty(self, get_worklist):
        for worklist in get_worklist():
            sample = deepcopy(worklist.names[0])
            sample.results["fit"] = {
                "test_reference": pd.DataFrame(
                    {
                        "gas": ["CO2"],
                        "reference": ["test_reference"],
                        "mmol_per_mg": [0.5],
                    }
                )
            }
            worklist_unit = Worklist([sample], name="unit_test")
            worklist_unit._results["robustness"] = pd.DataFrame(
                {"value": [1.0]},
                index=pd.MultiIndex.from_tuples([(1, 2)], names=["foo", "bar"]),
            )

            results = worklist_unit.results
            assert isinstance(results["fit"], pd.DataFrame)
            assert "mmol_per_mg" in results["robustness"].columns
            assert set(results["robustness"].index.names) == {"foo", "bar", "varied"}

    def test_getitem_with_list_returns_worklist_or_sample(self, get_worklist):
        for worklist in get_worklist():
            if len(worklist.names) > 1:
                multiple = worklist[[0, 1]]
                assert isinstance(multiple, Worklist)
            else:
                single = worklist[[0]]
                assert isinstance(single, Sample)

    