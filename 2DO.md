# Errors
- `dry_weight`
- `.corr(...)` if `plot=dict`, fill missing keys
- saving of `Worklist` doesn't seem to update the samplelog

# Further Improvements
## Config
- make config a function that reads from / writes to each time
    - use it everytime

## Import
- add `Worklist.info` (
    - `dict` of `dicts`
    -  `pd.DataFrame`)
-  get unit from header with `re`
- initialize ``Worklist`` from list of names to pass to ``Sample``

## Fitting
- `.robustness()` for single `Sample`
    - Workaround: Make Worklist from Single `Sample` to perform `.robustness()`
- `.fit()` plot 
    - remove error x-axis and add tick-marks at the top of the subplot
    - move **SQERR** to bottom 
    - move _SOG_-names to legend 
- remove extra markers from _fit-plot_
- parallelize computations

## Plotting
- `.plot()` add corrected plot. `Baseline`, `.raw` and corrected data are all available
- add better summary plots for fitting and robustness-results (e.g. with tidy data in combination with `seaborn`)
- make plot interactive
  - measure mass-, temperature- or time-differences
- plotting of calibration data
- add DTG to mass-stops (optional)

## General
- reduce amount of logging
    - adjust level for initialization substeps from *INFO* to *DEBUG*
- don't log `Baseline`
- add automated testing
- let `pint` and `molmass` handle units
    - `pint-pandas` for unit inside `pd.DataFrame`
- use minutes internally
- rename *IR* to *EGA*
- add `**kwargs` to every function for flexibility
- remove functions / REs specific to BAM-devices
- use root-folder names as profiles?

## Other
- `dry_weight(step_temp=..., mass_steps=..., step_time=...)` as arguments
- ? shift `Baseline` 

# Documentation
- for _import profiles_
- combining of objects
    - `Sample`+`Sample`=`Worklist`
    - `Sample`+`Worklist`=`Worklist`
    - `Worklist`+`Worklist`=`Worklist`
- calibration