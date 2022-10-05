# TGA-FTIR-hyphenation-tool-kit

A package for handling hyphenated TGA and FTIR data. Includes basic plotting as well as as advanced deconvolution of FTIR data.
For installation with conda and pip:

First create fresh environment with name <YOUR_ENV>:

`conda create -n <YOUR_ENV> python`

`conda activate <YOUR_ENV>`

Then run:

`pip install git+https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit.git`

To run GUI, download `run.py`, navigate to folder in your command prompt and run:

`python run.py`

---
---
## Main Usage
Core classes of the package are `Sample` and `Worklist`, where a `Worklist` object holds multiple `Sample` objects.

[1. Initialization](#1-initialization)\
[2. Correction](#2-correction)\
[3. Fitting](#3-fitting)\
[4. Robustness of fit](#4-robustness-of-fit)\
[5. Plotting](#5-plotting)\
[6. Saving](#6-saving)\
[7. Loading](#7-loading)

---
### 1. Initialization
`obj = Sample(<SAMPLE_NAME>)`
Optional arguments:

    alias           alias to use in e.g. plots
    mode            one of ["construct", "pickle"] to either:
                    - construct Sample from raw data
                    - load Sample from pickle
    profile         name of data loading profile. Can be set in
                    settings.ini

For more info on how to save as pickle, see [6. Saving](#6-saving).

Class attributes:

    name            name of Sample, default is filename used for initialization
    alias           alias of Sample, defaults to Sample.name
    sample          sample type used for grouping, defauls to 
                    Sample.name
    info = None     contains sample information e.g. initial_mass
    tga = None      contains tga data
    ir = None       contains ir data
    linreg = None   contains calibration data
    raw             conatins raw data from first initialization,
                    remains unchanged after e.g. correction. Is of 
                    type Sample    
---
after correction:

    baseline        contains data of baseline used for correction. 
                    Is of type Baseline which inherits from Sample

`wl = Worklist([obj], name = <WORKLIST_NAME>)`

Class attributes: 

    samples         List of Sample objects
    name            name of Worklist
    
For overview on Worklist use print(Worklist).


---
### 2. Correction
`obj.corr(<BASELINE_NAME> | None)`

`wl.corr(<BASELINE_NAME> | [<BASELINE_NAMES>] | None)`

If argument = None, it is lookes for a refernece given in the samplelog.

Optional arguments:

    plot = False        plot corrections

---
### 3. Fitting
`obj.fit(<REFERENCE_NAME>)`

`wl.fit(<REFERENCE_NAME>)`

Valid values for `<REFERENCE_NAME>` can be found and/ or set in `Fitting_parameter.xlsx` which is accessible trough ``TGA_FTIR_tools.fit_preferences()``.

Optional arguments:

    T_max = None        limit fit interval
    T_max_tol = 50      tolerance of center value
    save = True         save results
    plot = True         plot results
    presets = None      pass presets other than those specified by <REFERENCE_NAME>
    mod_sample = True   modify object during fitting

---
### 4. Robustness of fit

`wl.robustness(<REFERENCE_NAME>)`

Optional arguments:

    T_max = None        limit fit interval
    save = True         save results
    var_T = 10          absolute variance of HWHM_max and center_0
    var_rel = 0.3       relative variance of HWHM_0 and height_0

---
### 5. Plotting

`obj.plot(<OPTION>)` 
options = ["TG", "heat_flow", "IR", "DIR", "cumsum", "IR_to_DTG", "fit"]

`wl.plot(<OPTION>)` options = ["TG", "heat_flow", "IR", "DIR", "cumsum", "IR_to_DTG", "fit", "robustness", "results"]

---
### 6. Saving
`obj.save()`

 Optional arguments:

    how        one of ["samplelog", "excel", "pickle"] to save    
                - obj.info in the samplelog
                - obj.info, obj.ir, obj.tga to excel-file
                - obj as pickle-file

`wl.save()`

Optional arguments:

    fname       filename, if None: wl.name

---
### 7. Loading
`wl.load(<FNAME>)`

Load Worklist from pickle-file as produced by [6. Saving](#6-saving).

For Sample, see [1. Initialization](#1-initialization) and the use of the `mode`-argument.