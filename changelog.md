# 2026-02
- add debug statements to `read_data` and fix logging level issue.
- adjusted default savgol `window_length_rel = .05` (previously `0.01`)
- bugfix in baseline correction (always subtracted self, when not specified)
- more logging for calibration

# 2025-11
- use `skiprows: str` in `load_data` to skip rows starting with given string
- fix bug in calibration without plotting

# 2025-09
- added option to `samplelog` to read/write specific sheets
- added `Worklist.from_samplelog` to create worklist from specific sheet in samplelog
- fixed bug in ``Sample.get_values`` when no units are present [#87](https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit/issues/87)
- update import profile to be more flexible
- moved `baseline_als` to `corrections.py`
- add ``corrections`` to import profile
- restructured `Sample` -> 
  - kwarg `mode=pickle` is replaced by ``Sample.from_pickle()``
  - `.save(how=...)` is replaced by `.to_...()`
- removed `DEGA`-plot
- restructuring of ``README.md``
- add `Worklist.calibrate()` to run calibration.
- add `min_r2` to calibration to warn about bad fits
- Improved logic in determining dry point in `dry_weight()`
- made correction more flexible by allowing custom functions to be passed to `.corr()`
- cleanup of code and logging
- update testing
- updated `settings.ini`
- changed `SQERR` to `RMSE`
- reorganise testing

# 2025-08-25
- improved integration plot from calibration
- fixed bug in calibration when molecular formula is invalid
- updated Netzsch import profile

# 2025-08-21
- added `plot_results` to `Worklist.plot`
- added `plot_results` to `Sample.plot`
- fixed ``Sample`` and `Worklist`.results

# 2025-08-20
- new import profile for Netzsch
- made _usecols_ in import profile more flexible
  - pass integer indices of columns
  - pass slices as strings _e.g._ "1:-4"
  - pass names of columns
- seperated renaming of columns and mapping of names to internal names
- improved plotting

# 2025-08-07
- bugfix testing
- fixed unit conversions
- reorganized calibration output

# 2025-08-02
- checked ``.fit``
- changed fit visuals
- trying to make worklist functions concurrent

## 2025-08-01
- made ``mass_step`` more robust with new technique for finding of bounds
- improved visualisation
- fixed bug when updating calibration file
- fixed saving of worklist to samplelog

## 2025-07-31
- added testing routines with `pytest`
- added profile keyword to `Worklist`
- fixed ``mass_step``
- adjusted TGA3 import profile
- made *DTG*-smoothing more flexible and dependant on input data
- estimated reasonable values for *DTG*-smoothing and -peak detection in ``mass_step``
- bugfixes and cleanup

## 2025-07-30
- make info derived from raw data more flexible
- added unit support for info derived from file
- make calibration compliant with units
- begin refactoring of ``calibrate``

## 2025-07-29
### 1
- initial setup for testing
- refactor ``dry_weight``
- refactor ``Sample.info``
- refactor ``Sample.__post_init__``
- restructured import profile
### 2
- remove unmaintained chempy from dependencies and replace with molmass
- fix calculation of total amount per element
- clean up ``config.py``
- added unit support in import_profile
- added ``pint_pandas`` to handly unit annotation

## 2025-07-28
- mostly renamed *IR* to *EGA* (`Sample.ir`, "IR" as arg)
- moved calibration plots to separate files
- add option 'calibration' to `Sample.plot`
- add possibility to initialize `Worklist` with name or list of names in addition to list of `Sample`
- add possibility to add `Sample` + `Sample` | `Worklist` 
- cleaned up `Sample.plot`
- cleaned up dependencies

## 2022-10-13
- fixed bug in linked groups
- changed output path of worklist.fit to ~/worklist.name
- fixed bug in plotting of fit
- improved plotting overall
- added reference to sample.\_\_repr__()
- refactoring of mass_step
- added plot option "mass_steps"

## 2022-10-06
- update README
- update Worklist.\_\_add__()

## 2022-10-05
- minor bugfixes regarding Worklist.save and Sample.fit
- updated README

## 2022-10-04
- fixed bug introduced by the use of pathlib.Path
- cleaned up settings.ini
- fixed encoding error

## 2022-09-28
- improved download of supplementary data

## 2022-09-26
- added changelog
- added Sample.raw to access raw-data
- removed chemical formulas from settings -> use of chempy.Substance.latex_name
- added Sample.baseline to access baseline-data
- deleted old IO-functions
- updated paths to use Path from pathlib