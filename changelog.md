## 2025-07-28
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