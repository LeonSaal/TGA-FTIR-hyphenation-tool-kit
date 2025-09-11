from .calibration import calibrate
from .classes import Baseline, Sample, Worklist, SampleInfo
#from .config import PLOTTING, fit_references
from .fitting import bar_plot_results, concatenate, robustness, summarize
from .input_output import samplelog
from .plotting import plots
import pint_pandas
import pint
pint_pandas.PintType.ureg.formatter.default_format = "P~"
ureg = pint.get_application_registry()
ureg.autoconvert_offset_to_baseunit =True
pint.set_application_registry(ureg)
ureg.setup_matplotlib()
ureg.mpl_formatter = "{:~P}"
