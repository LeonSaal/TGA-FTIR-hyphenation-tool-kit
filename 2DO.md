| Category      | Task                                                                                   | Priority | First Noted |
| ------------- | -------------------------------------------------------------------------------------- | -------- | ----------- |
| Errors        | `.corr(...)` if `plot=dict`, fill missing keys                                         | Medium   |             |
| Config        | Make config a function that reads from / writes to each time                           | Medium   |             |
| Config        | Encoding of `settings.ini` (requires *ANSI* atm)                                       | Low      |             |
| Import        | Helper to determine `kwargs` for `read_csv` (auto encoding, delimiter detection, etc.) | High     |             |
| Import        | Move device specific settings to profile (e.g. mass resolution, savgol-settings)       | Medium   |             |
| Import        | Select profile on init and ask to save as default in settings                          | Medium   |             |
| Import        | Restructure import profile to distinguish required and optional fields                 | Medium   |             |
| Import        | Merge log messages from loops for fewer outputs                                        | Low      |             |
| Import        | Calculate `_info`, `.dry_weight` on the fly in `.info`-property                        | Medium   |             |
| Import        | Calculate file hash of import profile to check for updates in remote repo              | Low      |             |
| Import        | Add more gases/groups to Netzsch (aromatics, acetic acid, carbonyls, etc.)             | Low      |             |
| Import        | Merge init operations to separate methods (e.g. `load_sample`)                         | Medium   |             |
| Import        | Take correction info from file header                                                  | Medium   |             |
| Import        | In Worklist, don't ask for profile each time on first init                             | Medium   |             |
| Fitting       | `.fit()` plot improvements (remove error x-axis, move SQERR, move SOG-names)           | Medium   |             |
| Fitting       | Remove extra markers from fit-plot                                                     | Low      |             |
| Fitting       | Parallelize computations                                                               | High     |             |
| Fitting       | Calculate statistical values and append to `.results["fit"]`                           | Medium   |             |
| Fitting       | Set default `plot=False`                                                               | Low      |             |
| Fitting       | Make arrows narrower and labels outside box area                                       | Low      |             |
| Fitting       | Use linearized gaussian model to enhance fitting performance                           | High     |             |
| Robustness    | `.robustness()` for single `Sample`                                                    | High     |             |
| Plotting      | Add corrected plot (Baseline, raw and corrected data)                                  | Medium   |             |
| Plotting      | Make plot interactive (measure mass-, temperature- or time-differences)                | Low      |             |
| Plotting      | Add DTG to mass-steps (optional)                                                       | Medium   |             |
| Plotting      | DTG for worklists (adhere to README.md)                                                | Medium   |             |
| Plotting      | Move all plotting to `.plot(...)`-method                                               | High     |             |
| Plotting      | For `Worklist` use Seaborn to plot all samples in one figure                           | Medium   |             |
| Plotting      | For `Worklist` use Sample.plot and add all to same axis                                | Medium   |             |
| General       | Reduce amount of logging (adjust level to DEBUG)                                       | Medium   |             |
| General       | Don't log `Baseline`                                                                   | Low      |             |
| General       | Add `**kwargs` to every function for flexibility                                       | Medium   |             |
| General       | Use root-folder names as profiles or set folder patterns in profile definition         | Low      |             |
| General       | Convert `setup.py` to `pyproject.toml`                                                 | Medium   |             |
| General       | Code formatting                                                                        | Low      |             |
| General       | Use [rich logging](https://rich.readthedocs.io/en/stable/logging.html)                 | Low      |             |
| General       | Settings for `units` int_ega = ega                                                     | Low      |             |
| General       | Auto-save objects after changes (global setting)                                       | Low      |             |
| General       | Remove unused functions                                                                | Medium   |             |
| General       | Unify kwargs for plotting and saving (make_dir / directory → parent_dir)               | Low      |             |
| Testing       | Add automated testing for initialization, plotting, fitting                            | High     |             |
| Testing       | Add automated testing for different input data                                         | High     |             |
| Testing       | Test all possible args and kwargs                                                      | High     |             |
| Other         | `dry_weight(step_temp=..., mass_steps=..., step_time=...)` as arguments                | Medium   |             |
| Calibration   | Add sample labels to points (if specified)                                             | Low      |             |
| Calibration   | Check unit for calibration methods other than "max"                                    | Medium   |             |
| Calibration   | Add option to pass gases to calibrate explicitly                                       | Medium   |             |
| Calibration   | Add option to load specific calibration (not newest)                                   | Medium   |             |
| Plotting      | Integration plot: add legend and reduce padding between subfigures                     | Medium   |             |
| General       | Add date column                                                                        | Low      |             |
| General       | Allow passing of worklist as positional argument                                       | Low      |             |
| General       | Polyorder = 2, rel_window_length = .01 (settings)                                      | Low      |             |
| Fitting       | Integration bounds (baseline not constant due to window width)                         | Medium   |             |
| General       | Add option for automatic correction (buoyancy reference ≈ a·log x + b)                 | Medium   |             |
| Installation  | Add auxiliary files to project definition for download during install                  | Medium   |             |
| Documentation | Add docstrings and signatures for every function                                       | High     |             |
| Documentation | Add example folder                                                                     | Medium   |             |
| Documentation | Add Calibration method max explanation                                                 | Medium   |             |
| Documentation | Explain robustness                                                                     | Medium   |             |
| Documentation | Reorganize/reorder Class-specific methods                                              | Low      |             |
| Documentation | Make `to_...` methods default for saving                                               | Medium   |             |
| Documentation | Escape `$` in README.md under settings.ini File section                                | Low      |             |