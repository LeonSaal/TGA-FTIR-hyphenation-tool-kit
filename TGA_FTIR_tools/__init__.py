from TGA_FTIR_tools.classes import TG_IR
from TGA_FTIR_tools.config import PLOTTING
from TGA_FTIR_tools.calibration import calibrate
from TGA_FTIR_tools.input_output.samplelog import samplelog
from TGA_FTIR_tools.input_output.general import overview
from TGA_FTIR_tools.plotting import plots

import matplotlib as plt
plt.rcParams.update({'font.size': PLOTTING.getint('font_size')})
plt.rcParams['figure.figsize'] = PLOTTING.getfloat('figure_width')/2.54,PLOTTING.getfloat('figure_height')/2.54
