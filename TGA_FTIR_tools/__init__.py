from .classes import TG_IR
from .config import PLOTTING
from .calibration import calibrate
from .input_output import samplelog
from .input_output import overview
from .plotting import plots
from .fitting import robustness
from .fitting import fits

# plot settings
import matplotlib as plt
plt.rcParams.update({'font.size': PLOTTING.getint('font_size')})
plt.rcParams['figure.figsize'] = PLOTTING.getfloat('figure_width')/2.54,PLOTTING.getfloat('figure_height')/2.54

# try loading calibration data
try:
    linreg,stats=calibrate()
except:
    pass
