# High Priority

| Category      | Task                                                                           | First Noted |
| ------------- |------------------------------------------------------------------------------ | ----------- |
| General | Use `.get()`-method where feasible with sensible defaults (config)                               |             |
| General | Better error messages (Daniels issues)                            |             |
| Robustness    | `.robustness()` for single `Sample`                                            |             |
| Plotting      | Move all plotting to `.plot(...)`-method                                       |             |
| Testing       | Add automated testing for initialization, plotting, fitting                    |             |
| Testing       | Add automated testing for different input data                                 |             |
| Testing       | Test all possible args and kwargs                                              |             |
| Documentation | Add docstrings and signatures for every function                               |             |
| General | Replace references to `Sample.tga` and `Sample.ega` where possible with `Sample.data`  |   Oct. 26          |
|Import| `mode = [auto, pickle, construct]`. After `init` auto-save pickle if not exists, else load pickle ||
|Worklist| Rename attributes `.names`, `.samples`, `.aliases` and/or make property| |
|Import| rework profile prompt| |
|Import| update profiles| |
|`get_value`| FTIR-Info not by default, only with `get_value`, automatic calibration<br>more `which`-values| |



# Medium Priority

| Category      | Task                                                                                   | First Noted |
| ------------- | -------------------------------------------------------------------------------------- | ----------- |
| Errors        | `.corr(...)` if `plot=dict`, fill missing keys                                         |             |
| Config        | Make config a function that reads from / writes to each time                           |             |
| Import        | Move device specific settings to profile (e.g. mass resolution, savgol-settings)       |             |
| Import        | Select profile on init and ask to save as default in settings                          |             |
| Import        | Restructure import profile to distinguish required and optional fields                 |             |
| Import        | Calculate `_info`, `.dry_weight` on the fly in `.info`-property                        |             |
| Import        | Merge init operations to separate methods (e.g. `load_sample`)                         |             |
| Import        | Take correction info from file header                                                  |             |
| Fitting       | `.fit()` plot improvements (remove error x-axis, move SQERR, move SOG-names)           |             |
| Fitting       | Integration bounds (baseline not constant due to window width)                         |             |
| Fitting       | Parallelize computations                                                       |             |
| Fitting       | Use linearized gaussian model to enhance fitting performance                   |             |
| Fitting       | Move saving of results to `Sample` level                   |             |
| Plotting      | Add corrected plot (Baseline, raw and corrected data)                                  |             |
| Plotting      | Add DTG to mass-steps (optional)                                                       |             |
| Plotting      | For `Worklist` use Seaborn to plot all samples in one figure                           |             |
| Plotting      | For `Worklist` use Sample.plot and add all to same axis                                |             |
| General       | Add `**kwargs` to every function for flexibility                                       |             |
| Installation       | Convert `setup.py` to `pyproject.toml`                                                 |             |
| Calibration   | Check unit for calibration methods other than "max"                                    |             |
| Calibration   | Add option to load specific calibration (not newest)                                   |             |
| Installation  | Add auxiliary files to project definition for download during install                  |             |
| Documentation | Add example folder                                                                     |             |
| Installation | Add example folder                                                                     |             |
| Documentation | Add Calibration method max explanation                                                 |             |
| Documentation | Explain robustness                                                                     |             |
| Documentation | Make `to_...` methods default for saving                                               |             |
| Other         | `dry_weight(step_temp=..., mass_steps=..., step_time=...)` as arguments                |             |
|     README   | `Worklist.get_value`   |             |
|      README  | List functions available for both `Sample` and `Worklist`  |             |
|     README   | Import Profile examples with comments  |             |
|   README     |   Explain inner working of functions more (with collapsable fields?) |             |
|    README    |  Check if every `arg`, `*args` and `**kwargs` are documented  |             |
|     Plotting   | `dry_weight` add plot  |             |
|      Plotting  |  `Sample.plot(ylim = ...)` subset gases to show  |             |




# Low Priority

| Category      | Task                                                                                   | First Noted |
| ------------- | -------------------------------------------------------------------------------------- | ----------- |
| Config        | Encoding of `settings.ini` (requires *ANSI* atm)                                       |             |
| Import        | Merge log messages from loops for fewer outputs                                        |             |
| Import        | Calculate file hash of import profile to check for updates in remote repo              |             |
| Import        | Add more gases/groups to Netzsch (aromatics, acetic acid, carbonyls, etc.)             |             |
| Fitting       | Remove extra markers from fit-plot                                                     |             |
| Fitting       | Set default `plot=False`                                                               |             |
| Fitting       | plot: Make arrows narrower and labels outside box area                                       |             |
| Fitting       | Calculate statistical values and append to `.results["fit"]`                           |             |
| General       | Don't log `Baseline`                                                                   |             |
| General       | Use root-folder names as profiles or set folder patterns in profile definition?         |             |
| General       | Code formatting                                                                        |             |
| General       | Settings for `units` int_ega = ega                                                     |             |
| General       | Auto-save objects after changes (global setting)                                       |             |
| General       | Unify kwargs for plotting and saving (make_dir / directory → parent_dir)               |             |
| General       | Allow passing of worklist as positional argument                                       |             |
| General       | Polyorder = 2, rel_window_length = .01 (settings)                                      |             |
| General       | Add option for automatic correction (buoyancy reference ≈ a·log x + b)                 |             |
| Calibration   | Add sample labels to points (if specified)                                             |             |
| Documentation | Reorganize/reorder Class-specific methods                                              |             |

# Glamour

| Category      | Task                                                                                   | First Noted |
| ------------- | -------------------------------------------------------------------------------------- | ----------- |
| General       | Use [rich logging](https://rich.readthedocs.io/en/stable/logging.html)                 |             |
| General       | rename package in `setuptools.setup(name=)` _e.g. **tgega**_          |             |
| General       | Reduce amount of logging (adjust level to DEBUG)                                       |             |
| General       | Remove unused functions                                                                |             |
| Calibration      | Integration plot: add legend and reduce padding between subfigures                     |             |
| Plotting      | Make plot interactive (measure mass-, temperature- or time-differences)                |             |
| Import        | Helper to determine `kwargs` for `read_csv` (auto encoding, delimiter detection, etc.) |             |
|Import.Helper| get `ext` from sample names||
||1. look for files<br>2. list suffixes<br>3. detect \<suffix> <br> 4. first 100 rows, determine `sep` and `decimal` <br>5. skiprows||

