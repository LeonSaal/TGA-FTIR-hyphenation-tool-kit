# TGA-FTIR-hyphenation-tool-kit

A package for handling hyphenated TGA and FTIR data. Includes basic plotting as well as as advanced deconvolution of FTIR data.
For installation with conda and pip:

First create fresh environment with name <YOUR_ENV>:

``conda create -n <YOUR_ENV> python``

``conda activate <YOUR_ENV>``

Then run:

``pip install git+https://github.com/LeonSaal/TGA-FTIR-hyphenation-tool-kit.git``

To run GUI, download ``run.py``, navigate to folder in your command prompt and run:

``python run.py``

---

---
## Main Usage
Core classes of the package are ``Sample`` and ``Worklist``, where a ``Worklist`` object holds multiple ``Sample`` objects.

---
### Initialization
``obj = Sample(<SAMPLE_NAME>, alias = <ALIAS>)``

Class attributes 

    info = None     contains sample information e.g. name
    tga = None      contains tga data
    ir = None       contains ir data
    linreg = None   conains calibration data

``wl = Worklist([obj], name = <WORKLIST_NAME>)``


---
### Correction
``obj.corr(<BASELINE_NAME>)``

``wl.corr(<BASELINE_NAME> | [<BASELINE_NAMES>])``

Optional arguments:

    plot = False        plot corrections

---
### Fitting
``obj.fit(<REFERENCE_NAME>)``

``wl.fit(<REFERENCE_NAME>)``

Optional arguments:

    T_max = None        limit fit interval
    T_max_tol = 50      tolerance of center value
    save = True         save results
    plot = True         plot results
    presets = None      pass presets other than those specified by <REFERENCE_NAME>
    mod_sample = True   modify object during fitting

---
### Robustness of fit

``wl.robustness(<REFERENCE_NAME>)``

Optional arguments:

    T_max = None        limit fit interval
    save = True         save results
    var_T = 10          absolute variance of HWHM_max and center_0
    var_rel = 0.3       relative variance of HWHM_0 and height_0

---
### Plotting

``obj.plot(<OPTION>)`` 
options = ["TG", "heat_flow", "IR", "DIR", "cumsum", "IR_to_DTG", "fit"]

``wl.plot(<OPTION>)`` options = ["TG", "heat_flow", "IR", "DIR", "cumsum", "IR_to_DTG", "fit", "robustness", "results"]
