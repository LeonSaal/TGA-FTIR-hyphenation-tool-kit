# Errors
- `.corr(...)` if `plot=dict`, fill missing keys

# Further Improvements
## Config
- make config a function that reads from / writes to each time
    - use it everytime

## Import
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
  - select profile on init?
  - ask to save as default in settings
  - or move init operations to separate methods (e.g. `load_sample`)?
  - what if not all data is available for all sample for a certain profile? 
    - $\rightarrow$ number of columns doesn't match number of supplied names
    - $\rightarrow$ don't allow lists, only dicts for clear assignment
  - warning if profile doesn't match data
  - determine gases in ega more flexible (currently by index)
  - restructure import profile to better distinguish required and optional fields

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
- calculate statistical values and append to ``.results["fit"]`` or add as new key:value-pair $\rightarrow$ maybe not necessary if plotting with `seaborn`
- default ``plot=False``

## Plotting
- `.plot()` add corrected plot. `Baseline`, `.raw` and corrected data are all available
- make plot interactive?
  - measure mass-, temperature- or time-differences
- add DTG to mass-stops (optional)

## General
- reduce amount of logging
    - adjust level for initialization substeps from *INFO* to *DEBUG*
- don't log `Baseline`
- let `pint` handle units
    - or convert values to base units and handle bare numbers internally
- add `**kwargs` to every function for flexibility
- remove functions / REs specific to BAM-devices
- use root-folder names as profiles?
  - or set folder pattern(s) in profile definition (top-level)
  - warn if none found
- ``setup.py`` $\rightarrow$ ``pyproject.toml``
- clean up `calibration.py`
  - reduce if statements
  - combine code blocks to and or move functions to separate .py-files
- formatting of code 
- use [rich logging](https://rich.readthedocs.io/en/stable/logging.html)

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
- can only pass initialized worklist guard
- add date column

# Documentation
- calibration
- add docstrings, signatures for every function
- add example folder