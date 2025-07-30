# Errors
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
- helper to determine `kwargs` for `read_csv`
	- auto 
		- encoding from chardet, python-magic
		- read first n lines
		- count delimiter, quotes, \\t
			- guess delim and skiprows
		- count numbers and seps to guess header length
			- after skiprows guess first numerical row
			- merge multiindex header to single row afterwards
	- CLI?
        - first ask for required kwargs
        - only if not successful ask for other kwargs 
		- update preview or errors
    - combination of above
  - move device specific settings to profile (e.g. mass resolution)
    - savgol-settings?
    - move fields specific to ``Sample``-level to combined profile

## Correction
- add synthetic baseline
  - ALS, linear, const
- make correction of $CO_2$ more flexible / less specific to first device

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
- make plot interactive?
  - measure mass-, temperature- or time-differences
- add DTG to mass-stops (optional)

## General
- reduce amount of logging
    - adjust level for initialization substeps from *INFO* to *DEBUG*
- don't log `Baseline`
- let `pint` handle units
    - add unit annotations where possible
- add `**kwargs` to every function for flexibility
- remove functions / REs specific to BAM-devices
- use root-folder names as profiles?
- ``setup.py`` $\rightarrow$ ``pyproject.toml``
- testing with ``pytest`` (+`tox`?)
- clean up `Sample`, `calibration.py`
  - reduce if statements
  - combine code blocks to and or move functions to separate .py-files
- formatting of code 

## Testing
- add automated testing for
  - inititialization
  - plotting
  - correction
  - fitting
- for different input data
- all possible args, kwargs

## Other
- `dry_weight(step_temp=..., mass_steps=..., step_time=...)` as arguments
- ? shift `Baseline` 

## Calibration
- add sample labels to points (if specified)
- check unit for calibration methods other than "max"

# Documentation
- for _import profiles_
- combining of objects
    - `Sample`+`Sample`=`Worklist`
    - `Sample`+`Worklist`=`Worklist`
    - `Worklist`+`Worklist`=`Worklist`
- calibration
- add docstrings, signatures for every function