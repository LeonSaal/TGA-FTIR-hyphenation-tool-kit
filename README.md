# TGA-FTIR Hyphenation Tool Kit

A Python package for handling hyphenated TGA and EGA data, including plotting and advanced deconvolution of EGA data.
- [TGA-FTIR Hyphenation Tool Kit](#tga-ftir-hyphenation-tool-kit)
  - [Installation](#installation)
  - [Quick Start](#quick-start)
    - [Home Environment](#home-environment)
    - [Core Classes](#core-classes)
    - [Basic Usage](#basic-usage)
  - [API Reference](#api-reference)
    - [Initialization](#initialization)
      - [Sample](#sample)
      - [Worklist](#worklist)
    - [Loading](#loading)
    - [Saving](#saving)
    - [Common Dunder Methods](#common-dunder-methods)
    - [Main Methods](#main-methods)
      - [Correction](#correction)
      - [Calibration](#calibration)
      - [Plotting](#plotting)
      - [Fitting](#fitting)
      - [Robustness Analysis](#robustness-analysis)
    - [Class-Specific Methods](#class-specific-methods)
      - [Sample Methods](#sample-methods)
      - [Worklist Methods](#worklist-methods)
  - [Additional Information](#additional-information)
    - [Data Loading Profile](#data-loading-profile)
      - [Data](#data)
      - [Supplementary](#supplementary)
      - [Corrections](#corrections)
    - [`settings.ini` File](#settingsini-file)
    - [`Fitting_parameter.xlsx` File](#fitting_parameterxlsx-file)

---

## Installation

Create a new environment and install via ``pip``:

```sh
conda create -n <YOUR_ENV> python
conda activate <YOUR_ENV>
pip install git+https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit.git
```
---

## Quick Start
### Home Environment
On first import of the package this will set your home directory to the current working directory.
A fresh [`settings.ini`](#settingsini-file) file as well as a file with fit preferences will be copied here. You can modify both files to your needs.

Then, you are prompted to set a data directory. This will set the root folder under which all your raw data files are stored. 
To initialize a ``Sample`` or ``Worklist``, it suffices to provide just the filename(s) (without path). The package will search for the files under the data directory.

To load the raw data into the respective objects, you need to provide a [data loading profile](#data-loading-profile). If none is provided during initialization, an interactive prompt will guide you through the creation of one. Profiles are saved in the installation folder under *settings/import_profiles/profiles*.

### Core Classes

- **Sample**: Represents a single sample.
- **Worklist**: Holds multiple `Sample` objects.

### Basic Usage
```Python
from TGA_FTIR_tools import Sample, Worklist
s = Sample("example_sample", profile="example_profile")  # Initialize a Sample
w = Worklist(["example_sample", "another_sample"], profile="example_profile")  # Initialize a Worklist
```
Then you can access the data via the class attributes and use the main methods described below, e.g.:
```Python
ega = s.ega # Access EGA data
print(s.info) # Print all Sample information
s.plot("TG") # Plot TG data
...
```
---

## API Reference
### Initialization
#### Sample

```python
Sample(name, alias=None, mode="construct", profile=None)
```

**Parameters**

| Name    | Type       | Default | Description                                                                              |
| ------- | ---------- | ------- | ---------------------------------------------------------------------------------------- |
| name    | str        | —       | Name of the sample (usually filename)                                                    |
| profile | str        | None    | Data loading profile, if None is supplied, launched interactive prompt to construct one. |
|         |            |         |                                                                                          |
| alias   | str        | None    | Alias for plots and grouping                                                             |
| sample  | str        | alias   | Sample type for grouping                                                                 |
| run     | str \| int | alias   | Sample run for grouping replicates                                                       |


**Sample Attributes**

| Attribute | Type             | Default                   | Description                            |
| --------- | ---------------- | ------------------------- | -------------------------------------- |
| name      | str              | —                         | Name of sample                         |
| alias     | str              | name                      | Alias for sample                       |
| sample    | str              | alias                     | Sample type for grouping               |
| run       | str \| int       | alias                     | Sample run for grouping replicates     |
| info      | dict             | —                         | Sample information (e.g. initial_mass) |
| tga       | pd.DataFrame     | None                      | TGA data                               |
| ega       | pd.DataFrame     | None                      | EGA data                               |
| linreg    | pd.dataFrame     | None                      | Calibration data                       |
| results   | dict             | "fit":{}, "robustness":{} | Fit and robustness results             |
| raw       | Sample           | None                      | Raw data from initialization           |
| baseline  | Baseline(Sample) | None                      | Baseline data after correction         |



---

#### Worklist

```python
Worklist(samples, name=None, profiles=[], aliases=[])
```
**Parameters**

| Name     | Type                                 | Default | Description                                                     |
| -------- | ------------------------------------ | ------- | --------------------------------------------------------------- |
| samples  | list[str \| Sample] \| Sample \| str | —       | List of Sample objects or names or single sample name or Sample |
| name     | str                                  | None    | Name of the worklist                                            |
| profiles | list \| str                                 | []      | Data loading profile(s)                                           |
| aliases  | list                                 | []      | Sample aliases                                                  |

Or from samplelog:
```python
Worklist.from_samplelog(sheet_name=0)
```
**Parameters**

| Name       | Type       | Default | Description                                                                      |
| ---------- | ---------- | ------- | -------------------------------------------------------------------------------- |
| sheet_name | str \| int | 0       | Name or integer number of sheet in `Samplelog.xlsx` to initialize worklist from. |


**Worklist Attributes**

| Attribute | Description                                                                                              |
| --------- | -------------------------------------------------------------------------------------------------------- |
| samples   | List of Sample objects                                                                                   |
| name      | Name of worklist                                                                                         |
| profiles  | Data loading profiles                                                                                    |
| aliases   | Sample aliases from initialization. For current alias use `print(Worklist)` or access `Sample` directly. |
| info      | Combined Sample.info for all samples in Worklist.           |
| results   | Combined Sample.results for all samples in Worklist.             |

---

### Loading

Load from [pickle file](#saving).

```python
Sample.from_pickle(name)
Worklist.from_pickle(fname)
```

---

### Saving

Save sample or worklist data to output directory (specified  under `settings.ini/[paths]/output`).

```python
Sample.save(how="samplelog"|"excel"|"pickle")
Worklist.save(how="samplelog"|"pickle")
```

| `how`     | shortcut           | Description                                                                                        |
| --------- | ------------------ | -------------------------------------------------------------------------------------------------- |
| samplelog | object.save()      | Writes information from object to `Samplelog.xlsx`.                                                |
| pickle    | object.to_pickle() | Saves object as pickle-file.                                                                       |
| excel     | Sample.to_excel()  | Writes `Sample.info`,  `Sample.tga` and `Sample.ega` from sample to separate sheets in excel-file. |

---


### Common Dunder Methods

Both `Sample` and `Worklist` implement several Python "dunder" (double underscore) methods to provide intuitive behavior:

- `print(obj)`: Returns a readable string representation of the object, useful for debugging and interactive sessions.
- `obj1 + obj2`: Allows combining objects using the `+` operator. For example, `Sample + Sample` or `Sample + Worklist` returns a new `Worklist` containing both.
- `len(obj)`: Returns the number of contained samples (`1` for `Sample`, number of samples for `Worklist`).
- `for sample in obj:`: Enables iteration over contained samples (a single-item iterator for `Sample`, all samples for `Worklist`).
- `obj[]` (`Worklist` only): Supports indexing and slicing to access individual samples or subsets. Works with single integer index as well as with a list of indices and also with a single sample name.

These methods make it easy to work with `Sample` and `Worklist` objects in a Pythonic way, supporting iteration, combination, and inspection.

---

### Main Methods

#### Correction

Correct TG and EGA data for baseline or blank measurements. For
correction on import see [data loading profile](#corrections).

```python
Sample.corr(baseline, corrs, **kwargs)
Worklist.corr(baselines, corrs, **kwargs)
```

**Parameters**

| Name      | Type                                           | Default | Description                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| --------- | ---------------------------------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| baseline  | Sample \| Baseline                             |         | Baseline to use for correction.                                                                                                                                                                                                                                                                                                                                                                                                                     |
| baselines | Sample \| Baseline \| list[Sample \| Baseline] |         | Baselines to use for correction.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| corrs     | dict                                           |         | Dict of structure ``{"ega": {signal: fun}, "tga": ...}``. ``fun`` is function with signature ``fun(x, y, **kwargs) -> z`` that takes in two numeric vectors  and returns one vector, all of the same length. ``x`` and ``y`` are the respective signal columns from the sample to be corrected and the associated baseline. <br>*E.g.* `{"ega":{"CO": lambda x, y, offs=0: x - y}}` performs simple baseline subtraction on the $\text{CO}$-signal. |
| kwargs    | dict                                           |         | Additional keyword-arguments passed on to custom functions. <br> *E.g.* `offs=1` to shift the baseline from the above row.                                                                                                                                                                                                                                                                                                                          |

---

#### Calibration
Quantitive information on evolved gases can be obtained by calibrating the EGA signal with samples of known mass and thermogravimetric behaviour. 

```python
Worklist.calibrate()
```

**Parameters**

| Name               | Type     | Default               | Description                                                                                                                                                                      |
| ------------------ | -------- | --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| molecular_formulas | dict     | {}                    | If EGA signal name is no valid molecular formula *e.g.* `"mz_44"`, the corresponding one can be passed here *e.g.* `{"mz_44": "CO2"}`.                                           |
| method             | str      | "max"                 | Calibration method. *max* refers to                                                                                                                                              |
| width_T            | np.array | np.array([0, np.inf]) | Width range of peaks in DTG to determine mass steps.                                                                                                                             |
| min_rel_height     | float    | 0.2                   | Minimum relative height of DTG peak to qualify mass step.                                                                                                                        |
| corr_baseline      | str      | "linear"              | Baseline correction for integration of EGA peak. <ul><li>linear = Linear baseline connecting edges of peak</li><li>const = Constant baseline at minimum intensity</li></ul>      |
| min_r2             | float    | 0.95                  | Minimum R² below which a warning message is triggered (Also at R²=1).                                                                                                            |
| plot               | bool     | False                 | Show integration and calibration plots.                                                                                                                                          |
| **fig_args         | dict     | {}                    | Extra arguments passed to [``matplotlib.pyplot.subplots``](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html#matplotlib.pyplot.subplots) if `plot=True`. |

#### Plotting

Plot sample or worklist data.

```python
Sample.plot(option, ax=None, save=False, **kwargs)
Worklist.plot(option, **kwargs)
```

**Options**

| Option        | Description            |
| ------------- | ---------------------- |
| "TG"          | TG data                |
| "mass_steps"  | Mass steps             |
| "heat_flow"   | Heat flow              |
| "EGA"         | EGA data               |
| "cumsum"      | Cumulative sum  of EGA       |
| "EGA_to_DTG"  | EGA to DTG             |
| "fit"         | Fit results            |
| "calibration" | Calibration data       |
| "results"     | Fit/robustness results |

**Common Parameters**

| Name   | Type | Default       | Description             |
| ------ | ---- | ------------- | ----------------------- |
| ax     | axis | None          | Matplotlib axis to plot on       |
| save   | bool | False         | Save figure             |
| x_axis | str  | "sample_temp" | "time" or "sample_temp" |
| y_axis | str  | "orig"        | "orig" or "rel"         |
| legend | bool | True          | Show legend             |
| title  | bool | True          | Show title              |

See function docstrings for option-specific parameters.

---

#### Fitting

Fit sample or worklist data.

```python
Sample.fit(reference_name, T_max=None, T_max_tol=50, save=True, plot=True, presets=None, mod_sample=True, overwrite=False)
Worklist.fit(reference_name, ...)
```

**Parameters**

| Name           | Type  | Default | Description                  |
| -------------- | ----- | ------- | ---------------------------- |
| reference_name | str   | —       | Reference for fitting. Refers to first column in [`Fitting_parameter.xlsx`](#fitting_parameterxlsx-file).        |
| T_max          | float | None    | Limit fit interval           |
| T_max_tol      | float | 50      | Tolerance for center value   |
| save           | bool  | True    | Save results                 |
| plot           | bool  | True    | Plot results                 |
| presets        | dict  | None    | Custom presets               |
| mod_sample     | bool  | True    | Modify object during fitting |
| overwrite      | bool  | False   | Overwrite existing fit       |

Results are stored in `results["fit"][reference_name]`.

---

#### Robustness Analysis

Assess sensitivity of fit to parameter variations.

```python
Worklist.robustness(reference_name, T_max=None, save=True, var_T=10, var_rel=0.3)
```

**Parameters**

| Name    | Type  | Default | Description           |
| ------- | ----- | ------- | --------------------- |
| var_T   | float | 10      | Absolute variance (K) |
| var_rel | float | 0.3     | Relative variance     |

---

### Class-Specific Methods

#### Sample Methods

| Method                      | Signature                                                                 | Description                                                                                                                                                                                         |
| --------------------------- | ------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `reference_mass` (property) | `Sample.reference_mass`                                                   | Returns reference mass value.                                                                                                                                                                       |
| `step_data`                 | `Sample.step_data(ref_mass_name=None)`                                    | Returns DataFrame with step data. Default steps include initial and final mass but can be extended with the dry mass or other steps as determined by `Sample.dry_weight()` or `Sample.mass_step()`. |
| `get_value`                 | `Sample.get_value(*values, which="sample_mass", at="sample_temp")`        | Extracts values from TG data at specified points.                                                                                                                                                   |
| `dry_weight`                | `Sample.dry_weight(plot=False, **kwargs)`                                 | Determines dry point and dry mass of sample. By default, if $\text{H}_2 \text{O}$-data is present, a peak within the minimum temperature and 200 °C is used to determine dryness, alternatively a peak in the DTG in the same temeperature range. |
| `mass_step`                 | `Sample.mass_step(ax=None, plot=True, min_rel_height=0.2, width_T=[0,∞])` | Detects mass steps in TGA data. ``min_rel_height`` sets lower limit of peaks in DTG, ``width_T`` controls peak width in K                                                                           |


---


#### Worklist Methods

| Method      | Signature                             | Description                                                                                                                                                                                                                                                 |
| ----------- | ------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `get`       | `Worklist.get(pattern, attr="name")`  | Returns a Worklist filtered with regular expression on attribute.                                                                                                                                                                                           |
| `append`    | `Worklist.append(Sample \| Worklist)` | Adds a Sample or another Worklist. In contrast to `+` operation modifies inplace.                                                                                                                                                                            |
| `pop`       | `Worklist.pop(i)`                     | Removes and returns sample at index `i`.                                                                                                                                                                                                                    |                                                                      |
| `load`      | `Worklist.load(fname)`                | Loads a Worklist from pickle file.                                                                                                                                                                                                                          |
| `calibrate` | `Worklist.calibrate(**kwargs)`        | Calibrates EGA data with Worklist. Worklist should contain measurements of the thermal degradation of different masses from a substance with known thermogravimetric behaviour and consecutive release of different gases. See [Calibration](#calibration). |

---

## Additional Information

- For fit references, see [`Fitting_parameter.xlsx`](#fitting_parameterxlsx-file) via `TGA_FTIR_tools.fit_preferences()`.
- For more details, see function docstrings and examples.

### Data Loading Profile
The hyphenation is achieved with a `.json`-file in the ``settings/import_profiles/profiles``-folder of the following structure:
```json
{
    "data": {
        "tga": "path to tga profile.json",
        "ega": "path to ega profile.json"
    },
    "supplementary": {
        "name_pattern": "",
        "mass_resolution_ug": ""
    },
    "corrections": {
        "ega": {
            "H2O": {
                "method": "arpls",
                "kwargs": {
                    ...
                }
            }
        }
    }
}
```
#### Data
For the *data*-field, for each of the hyphenated devices a path leading to a file of the following structure has to be provided in the respective *tga* and *ega* folders.

Each file must contain the two topmost and can contain other otional key-value-pairs.

| key          | Type            | Default                                                             | Description                                                                                                                                                                                                                                                                                                                                                                         |
| ------------ | --------------- | ------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| spec         | dict            | "manufacturer": str, "device": str, "software": str, "version": str | Metainformation on the device used to collect data.                                                                                                                                                                                                                                                                                                                                 |
| ext          | str             | —                                                                   | Suffix of raw data files. Can be either just the file extension *e.g.* ".csv" or a regular expression containing the capture group *suffix e.g.*: "_?(?P\<suffix\>.*).txt". If a column in the raw data gets assigned the name "suffix", this is later replaced with the captured suffix.                                                                                           |
|              |                 |                                                                     |                                                                                                                                                                                                                                                                                                                                                                                     |
| kwargs       | dict            | —                                                                   | Extra arguments passed on to [`pandas.read_csv`](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html#pandas.read_csv). Currently only data readable by this function is supported. For the `skiprows` argument, you can pass a string marking the start of the last row you want to skip. This gives more flexibility for different data exports from the same device.                                                                                                                                                                                   |
| map_suffix   | dict            | —                                                                   | Map values of capture group "suffix" (see *ext*) to custom name.                                                                                                                                                                                                                                                                                                                    |
| info_pattern | dict            | —                                                                   | Key-regular expressions (regex) pairs to extract additional information from raw data *e.g.* from the header or footer of the file. The regex must have named capture groups: "name", "value" and "unit" (can be empty). If the regex matches once, the key from the dict is used, otherwise the name from the capture group.                                                       |
| usecols      | dict            | —                                                                   | Subset and rename columns in combined data. Dict keys can contain a combination of integer indices, column names and column name patterns. The corresponding value (if not `null`) is the replacement value. For patterns, the replacement is done using [`re.sub`](https://docs.python.org/3/library/re.html#re.sub), so capture groups can be accessed via *e.g.* number or name. |
| units        | str, list, dict | —                                                                   | Assign units to columns. If str, regex with capture-group "unit" to extract unit from column header. If list, unit for respective column (must have same length). If dict, a mapping of column-name: unit.                                                                                                                                                                          |

For examples, see [settings/import_profiles](https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit/tree/main/TGA_FTIR_tools/settings/import_profiles). The default profile is saved in the [settings-File](#settingsini-file).


#### Supplementary
The *supplementary* field is currently not used, but should incorporate settings specific to the respective hyphenation in the future.

#### Corrections
The *corrections*-field can contain correction methods to be applied during initialization of a `Sample` object. By default, only correction of water in EGA data is performed due to commonly occuring drift issues. The `kwargs` field can be used to adjust method-specific parameters.

| Method  | Description                                                                                                                 | kwargs with defaults                   | kwargs-Description                                                                                                                                                                                      |
| ------- | --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `arpls` | [Baseline correction using asymmetrically reweighted penalized least squares smoothing](https://doi.org/10.1039/C4AN01061B) | lam=1e5<br> ratio=1e-6<br> iter_max=50 | <ul><li>lam = smoothness parameter (should be varied by orders of magnitude)</li> <li> ratio = early termination criterion of weight-change </li><li>iter_max = maximum number of iterations.</li></ul> |
| `als`   | Baseline correction using asymmetrically least squares smoothing (see link above).                                          | lam=1e6<br> p=0.01<br> iter_max=10     | <ul><li>p = parameter to set asymmetry in weights for least squares (recommended range 0.1 - 0.001).</li></ul>                                                                                          |
| `min`   | Using minimum value as constant baseline.                                                                                   | —                                      |                                                                                                                                                                                                         |
| `const` | Estimate noise level by subtracting Savitzky-Golay-smoothed from original signal and constructing constant baseline.        | spread=1                               | <ul><li>spread = factor to scale noise estimate.</li></ul>                                                                                                                                              |
| `none`  | Baseline is set to zero.                                                                                                    | —                                      |                                                                                                                                                                                                         |


### `settings.ini` File
The `settings.ini` file is created in the home directory upon first import of the package. It contains paths and default settings for various functions. You can modify these settings to your needs.

| Section    | key               | default               | Description                        |   
| ---------- | ----------------- | --------------------- | ---------------------------------- | 
| **paths**      |                   |                       |  Directories for organisation of data.                                  |  
|       |                   |                       |                                    |  
|            | home              |                       | Working directory.                 |  
|            | data              |                       | Data directory.                    |
|            | calibration       | %(home)s/Calibration  | Directory for calibration results. |
|            | output            | %(home)s/Output       | Directory for file output.         |
|            | plots             | %(output)s/Plots      |                                    |
|            | fitting           | %(output)s/Fitting    |                                    |
|            | robustness        | %(output)s/Robustness |                                    |
|       |                   |                       |                                    |  
| **savgol**     |                   |                       |  Smoothing parameters for DTG                                |  
|       |                   |                       |                                    |  
|            | window_length_rel | 0.01                   |     Relative smoothing window length. For the absolute window length, the length of data is considered.                               |
|            | polyorder         | 2                     |        Polynomial order of filter.                            |
|       |                   |                       |                                    |  
| **defaults**   |                   |                       |                                    |
|       |                   |                       |                                    |  
|            | logging_level     | INFO                  |     Set [logging level](https://docs.python.org/3/library/logging.html#logging-levels).                                |
|            | profile           |                       |         Default [import profile](#data-loading-profile)                           |
|       |                   |                       |                                    |  
| **correction** |                   |                       |                                    |
|       |                   |                       |                                    |  
|            | tga_plot          | True                  |                                    |
|            | dry_weight_plot   | True                  |                                    |
|            | ega_plot          | True                  |                                    |
|       |                   |                       |                                    |  
| **plotting**   |                   |                       |         Plot appearance defaults.                           |
|       |                   |                       |                                    |  
|            | mpl-style         | default               |           [matplotlib style](https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html)                         |
|            | x_axis            | sample_temp           |              Default x_axis, must be present in data. Another option can be "time".                      |
|            | y_axis            | orig                  |             Default x_axis. <br><ul><li>orig = data in original units</li><li>rel = data in relative values (%)</li></ul>                       |
|            | xlim              | None,None             |     Set x-limits of plot.                               |
|            | ylim              | auto                  |     Set y-limits of plot.                                          |
|            | save              | False                 |         Save plots.                           |
|            | title             | False                 |         Show title.                           |
|            | legend            | True                  |         Show legend.                           |
|       |                   |                       |                                    |  
| **labels**     |                   |                       |       Labels for parameters found in data. Styling with $\LaTeX$ possible *e.g.* \$m_0\$ $\rightarrow$ $m_0$                             |
|       |                   |                       |                                    |  
|            | sample_mass       | mass                  |                                    |
|            | time              | t                     |                                    |
|            | sample_temp       | T                     |                                    |
|            | dtg               | mass loss rate        |                                    |
|            | molar_amount      | n                     |                                    |
|            | heat_flow         | $\dot{Q}$             |                                    |
|            | ega               | EGA                   |                                    |
|            | dry_mass          | $m_{dry}$             |                                    |
|            | initial_mass      | $m_0$                 |                                    |
|            | final_mass        | $m_{end}$             |                                    |
|       |                   |                       |                                    |  
| **units**      |                   |                       |     Units for plots. **Currently not unit-safe!**                               |
|       |                   |                       |                                    |  
|            | sep               | /                     |        Seperator between parameter and unit.                            |
|            | int_ega           | A                     |           EGA-intensity.                         |
|            | sample_mass       | mg                    |                                    |
|            | time              | min                   |                                    |
|            | sample_temp       | °C                    |                                    |
|            | molar_amount      | mmol                  |                                    |
|            | heat_flow         | mW                    |                                    |
|            | dtg               | mg\,min^{{-1}}        |                                    |
|       |                   |                       |                                    |  
| **fitting**    |                   |                       |       Default fitting parameter if not specified in [`Fitting_parameter.xlsx`](#fitting_parameterxlsx-file).                             |
|       |                   |                       |                                    |  
|            | tol_center        | 30.0                  |  Allow center to be moved $\pm 30 K$ arount initial value.                                  |
|            | hwhm_0            | 47.5                  |   Initial HWHM.                                 |
|            | hwhm_min          | 0                     |       Lower bound HWHM.                             |
|            | hwhm_max          | 95                    |          Upper bound HWHM.                          |
|            | height_0          | 0.5                   |      Initial relative height $\frac{h}{\max(h)}$                          |
|            | height_min        | 0                     |   Lower bound relative height.                                 |
|            | height_max        | 1                     |    Upper bound relative height.                                 |

### `Fitting_parameter.xlsx` File

Structure of the `Fitting_parameter.xlsx` file to define fitting preferences. The first column, `group`, contains the name of the reference, which is passed to the `.fit()` or `.robustness()` method. The other columns contain the names of the signals to be fitted and their initial center values (in °C). In the second row ,`gas`, the gas for the respective group is specified. This must correspond to a column name in the EGA data (*e.g.* `CO2` or `mz_44`).

Additional sheets can contain preset values for `hwhm_0`, `hwhm_min`, `hwhm_max`, `height_0`, `height_min` and `height_max`. If no preset is given, the default from the [settings file, section **fitting**](#settingsini-file) is used.

| group<br>gas | adsorbed CO2<br>CO2 | carboxylic acids_1<br>CO2 |  ...   |
| ------------ | ------------------- | ------------------------- | --- |
| reference0  |                     |                        |   ...  |
| reference1 |                     |                        |   ...  |
| ... |            ...         |...|   ...  |

See [default file in repository](https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit/blob/main/TGA_FTIR_tools/settings/Fitting_parameter.xlsx).