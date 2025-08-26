# TGA-FTIR Hyphenation Tool Kit

A Python package for handling hyphenated TGA and EGA data, including plotting and advanced deconvolution of EGA data.
- [TGA-FTIR Hyphenation Tool Kit](#tga-ftir-hyphenation-tool-kit)
  - [Installation](#installation)
  - [Quick Start](#quick-start)
    - [Home Environment](#home-environment)
    - [Core Classes](#core-classes)
    - [Basic Usage](#basic-usage)
  - [API Reference](#api-reference)
    - [Sample](#sample)
      - [Initialization](#initialization)
    - [Worklist](#worklist)
      - [Initialization](#initialization-1)
  - [Main Methods](#main-methods)
    - [Correction](#correction)
    - [Fitting](#fitting)
    - [Robustness Analysis](#robustness-analysis)
    - [Plotting](#plotting)
    - [Saving](#saving)
    - [Loading](#loading)
  - [Common Dunder Methods](#common-dunder-methods)
  - [Class-Specific Methods](#class-specific-methods)
    - [Sample Methods](#sample-methods)
    - [Worklist Methods](#worklist-methods)
  - [Additional Information](#additional-information)
    - [Data Loading Profile](#data-loading-profile)

---

## Installation

Create a new environment and install via pip:

```sh
conda create -n <YOUR_ENV> python
conda activate <YOUR_ENV>
pip install git+https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit.git
```
---

## Quick Start
### Home Environment
On first import of the package this will set your home directory to the current working directory.
A fresh `settings.ini` file as well as a file with fit preferences will be copied here. You can modify both files to your needs.

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

### Sample

#### Initialization

```python
Sample(name, alias=None, mode="construct", profile=None)
```

**Parameters**

| Name    | Type | Default     | Description                                                                              |
| ------- | ---- | ----------- | ---------------------------------------------------------------------------------------- |
| name    | str  | —           | Name of the sample (usually filename)                                                    |
| alias   | str  | None        | Alias for plots and grouping                                                             |
| mode    | str  | "construct" | "construct" (from raw data) or "pickle" (from [pickle](#saving))                         |
| profile | str  | None        | Data loading profile, if None is supplied, launched interactive prompt to construct one. |

**Attributes**

| Attribute | Type             | Default                     | Description                            |
| --------- | ---------------- | --------------------------- | -------------------------------------- |
| name      | str              | —                           | Name of sample                         |
| alias     | str              | name                        | Alias for sample                       |
| sample    | str              | —                           | Sample type for grouping               |
| info      | dict             | —                           | Sample information (e.g. initial_mass) |
| tga       | pd.DataFrame     | None                        | TGA data                               |
| ega       | pd.DataFrame     | None                        | EGA data                               |
| linreg    | pd.dataFrame     | None                        | Calibration data                       |
| results   | dict             | {"fit":{}, "robustness":{}} | Fit and robustness results             |
| raw       | Sample           | None                        | Raw data from initialization           |
| baseline  | Baseline(Sample) | None                        | Baseline data after correction         |

---

### Worklist

#### Initialization

```python
Worklist(samples, name=None, profile=None)
```

**Parameters**

| Name    | Type        | Default | Description                     |
| ------- | ----------- | ------- | ------------------------------- |
| samples | list/Sample | —       | List of Sample objects or names |
| name    | str         | None    | Name of the worklist            |
| profile | str         | None    | Data loading profile            |

**Attributes**

| Attribute | Description            |
| --------- | ---------------------- |
| samples   | List of Sample objects |
| name      | Name of worklist       |
| profile   | Data loading profile   |

---

## Main Methods

### Correction

Correct TG and EGA data for baseline or blank measurements.

```python
Sample.corr(baseline_name=None, plot=False)
Worklist.corr(baseline_names=None, plot=False)
```

**Parameters**

| Name           | Type          | Default | Description             |
| -------------- | ------------- | ------- | ----------------------- |
| baseline_name  | str/None      | None    | Baseline name or None   |
| baseline_names | str/list(str) | None    | Baseline name/s or None |
| plot           | bool          | False   | Plot corrections        |

---

### Fitting

Fit sample or worklist data.

```python
Sample.fit(reference_name, T_max=None, T_max_tol=50, save=True, plot=True, presets=None, mod_sample=True, overwrite=False)
Worklist.fit(reference_name, ...)
```

**Parameters**

| Name           | Type  | Default | Description                  |
| -------------- | ----- | ------- | ---------------------------- |
| reference_name | str   | —       | Reference for fitting        |
| T_max          | float | None    | Limit fit interval           |
| T_max_tol      | float | 50      | Tolerance for center value   |
| save           | bool  | True    | Save results                 |
| plot           | bool  | True    | Plot results                 |
| presets        | dict  | None    | Custom presets               |
| mod_sample     | bool  | True    | Modify object during fitting |
| overwrite      | bool  | False   | Overwrite existing fit       |

Results are stored in `results["fit"][reference_name]`.

---

### Robustness Analysis

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

### Plotting

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
| "DEGA"        | DEGA data              |
| "cumsum"      | Cumulative sum         |
| "EGA_to_DTG"  | EGA to DTG             |
| "fit"         | Fit results            |
| "calibration" | Calibration data       |
| "results"     | Fit/robustness results |

**Common Parameters**

| Name   | Type | Default       | Description             |
| ------ | ---- | ------------- | ----------------------- |
| ax     | axis | None          | Matplotlib axis         |
| save   | bool | False         | Save figure             |
| x_axis | str  | "sample_temp" | "time" or "sample_temp" |
| y_axis | str  | "orig"        | "orig" or "rel"         |
| legend | bool | True          | Show legend             |
| title  | bool | True          | Show title              |

See function docstrings for option-specific parameters.

---

### Saving

Save sample or worklist data.

```python
Sample.save(how="samplelog"|"excel"|"pickle")
Worklist.save(fname=None)
```

---

### Loading

Load worklist from pickle file.

```python
Worklist.load(fname)
```


---

## Common Dunder Methods

Both `Sample` and `Worklist` implement several Python "dunder" (double underscore) methods to provide intuitive behavior:

- `print(obj)`: Returns a readable string representation of the object, useful for debugging and interactive sessions.
- `obj1 + obj2`: Allows combining objects using the `+` operator. For example, `Sample + Sample` or `Sample + Worklist` returns a new `Worklist` containing both.
- `len(obj)`: Returns the number of contained samples (`1` for `Sample`, number of samples for `Worklist`).
- `for sample in obj:`: Enables iteration over contained samples (a single-item iterator for `Sample`, all samples for `Worklist`).
- `obj[]` (`Worklist` only): Supports indexing and slicing to access individual samples or subsets. Works with single integer index as well as with a list of indices and also with a single sample name.

These methods make it easy to work with `Sample` and `Worklist` objects in a Pythonic way, supporting iteration, combination, and inspection.

## Class-Specific Methods

### Sample Methods

| Method                      | Signature                                                                 | Description                                                                                                                                                                                                                |
| --------------------------- | ------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `info` (property)           | `Sample.info`                                                             | Returns sample information.                                                                                                                                                                                                |
| `reference_mass` (property) | `Sample.reference_mass`                                                   | Returns reference mass value.                                                                                                                                                                                              |
| `step_data`                 | `Sample.step_data(ref_mass_name=None)`                                    | Returns DataFrame with step data. Default steps include initial and final mass but can be extended with the dry mass or other steps as determined by `Sample.dry_weight()` or `Sample.mass_step()`.                        |
| `get_value`                 | `Sample.get_value(*values, which="sample_mass", at="sample_temp")`        | Extracts values from TG data at specified points.                                                                                                                                                                          |
| `dry_weight`                | `Sample.dry_weight(plot=False, **kwargs)`                                 | Determines dry point and mass of sample.                                                                                                                                                                                   |
| `mass_step`                 | `Sample.mass_step(ax=None, plot=True, min_rel_height=0.2, width_T=[0,∞])` | Detects mass steps in TGA data. ``min_rel_height`` sets lower limit of peaks in DTG, ``width_T`` controls peak width in K                                                                                                  |
| `calibrate`                 | `Sample.calibrate(worklist = Worklist,profile=None, **kwargs)`            | Calibrates EGA data with Worklist. Worklist should contain measurements of the thermal degradation of different masses from a substance with known thermogravimetric behaviour and consecutive release of different gases. |

---


### Worklist Methods

| Method   | Signature                             | Description                                                                      |
| -------- | ------------------------------------- | -------------------------------------------------------------------------------- |
| `get`    | `Worklist.get(pattern, attr="name")`  | Returns a Worklist filtered with regular expression on attribute.                |
| `append` | `Worklist.append(Sample \| Worklist)` | Adds a Sample or another Worklist.In contrast to `+` operation modifies inplace. |
| `pop`    | `Worklist.pop(i)`                     | Removes and returns sample at index `i`.                                         |
| `info`   | `Worklist.info(type="df")`            | Returns info for all samples as DataFrame or dict.                               |
| `load`   | `Worklist.load(fname)`                | Loads a Worklist from pickle file.                                               |

---

## Additional Information

- For fit references, see `Fitting_parameter.xlsx` via `TGA_FTIR_tools.fit_preferences()`.
- For more details, see function docstrings and examples.

### Data Loading Profile

For each of the hyphenated devices a separate file of the following structure has to be provided in the respective *tga* and *ega* folders.

The file must contain two required and can contain other otional key-value-pairs.

|    | key  | Type | Default | Description     |
| --- | ---- | ---- | ------- | --------------- |
| required    | spec   | dict | {"manufacturer": str, "device": str, "software": str, "version": str}   | Metainformation on the device used to collect data. |
|  required   | ext | str |  —  | Suffix of raw data files. Can be either just the file extension *e.g.* ".csv" or a regular expression containing the capture group *suffix e.g.*: "_?(?P\<suffix\>.*).txt". If a column in the raw data gets assigned the name "suffix", this is later replaced with the captured suffix.    |
|optional|kwargs|dict|—|Extra arguments passed on to [`pandas.read_csv`](https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html#pandas.read_csv). Currently only data readable by this function is supported.|
|optional|map_suffix|dict|—|Map values of capture group "suffix" to custom name.|
|optional|info_pattern|dict|—|Key-regular expressions (regex) pairs to extract additional information from raw data *e.g.* from the header or footer of the file. The regex must have named capture groups: "name", "value" and "unit" (can be empty). If the regex matches once, the key from the dict is used, otherwise the name from the capture group.|
|optional|usecols|list|—|Subset columns in combined data. List can contain a combination of integer indices, column names (after renaming) and slices as string *e.g.* "1:3"|
|optional|units|str, list, dict|—|Assign units to columns. If str, regex with capture "unit" to extract unit from column header. If list, unit for respective column (must have same length). If dict, a mapping of column-name: unit. |
|optional|rename|str, list, dict|—|enaming of columns. If str, assumes function to be evaluated *e.g.* ``lambda x: x.upper()``. If list, new name for respective column (must have same length). If dict, a mapping of old-name:new-name|
|optional|name_mapping|dict|—|A mapping of old-name:new-name (old-name after rename operation)|

For examples, see [settings/import_profiles](https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit/tree/main/TGA_FTIR_tools/settings/import_profiles).

The hyphenation is achieved with a separate `.json` of the following structure:
```json
{
    "data": {
        "tga": "path to tga profile.json",
        "ega": "path to ega profile.json"
    },
    "supplementary": {
        "name_pattern": "",
        "mass_resolution_ug": ""
    }
}
```
The supplementary field is currently not used, but should incorporate settings specific to the respective hyphenation in the future.