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
    - or have global default overwritten by profile
  - select profile on init?
  - ask to save as default in settings
  - or move init operations to separate methods (e.g. `load_sample`)?
  - restructure import profile to better distinguish required and optional fields
  - merge log messages from loops for fewer outputs
  - calculate ``_info``, ``.dry_weight`` usw. on the fly in ``.info``-property 

## Correction
- make correction of $CO_2$ more flexible / less specific to first device
  - expect corrected profiles?
- pass custom function with given signature to `.corr()`

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
- arrows narrower and labels outside box area

## Plotting
- `.plot()` add corrected plot. `Baseline`, `.raw` and corrected data are all available
- make plot interactive?
  - measure mass-, temperature- or time-differences
- add DTG to mass-steps (optional)
- DTG for worklists (adhere to `README.md`)
- move all pltting to `.plot(...)`-method

## General
- reduce amount of logging
    - adjust level for initialization substeps from *INFO* to *DEBUG*
- don't log `Baseline`
- add `**kwargs` to every function for flexibility
- use root-folder names as profiles?
  - or set folder pattern(s) in profile definition (top-level)
  - warn if none found
- ``setup.py`` $\rightarrow$ ``pyproject.toml``
- formatting of code 
- use [rich logging](https://rich.readthedocs.io/en/stable/logging.html)
- ``settings[units] int_ega = ega``
- auto-save objects after changes -> global setting
- remove functions that are not used anymore

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

## Calibration
- add sample labels to points (if specified)
- check unit for calibration methods other than "max"
- add date column
- passing of worklist as positional argument? 
- polyorder = 2, rel_window_lenght = .01
- integration bounds (baseline not constant due to window width)

# Documentation
- add docstrings, signatures for every function
- add example folder