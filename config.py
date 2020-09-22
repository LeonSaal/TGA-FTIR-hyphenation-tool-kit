import configparser
config = configparser.ConfigParser()
config.read('settings.ini')

PATHS=config['paths']

UNITS=config['units']
keys=['sample_mass','time','sample_temp','molar_amount','heat_flow']
units=['mg','min','°C','mmol','mW']
for key, val in zip(keys,units):
    UNITS[key]=val
    
SEP=UNITS['sep']

PARAMS=config['parameters']

MOLAR_MASS=config['molar_mass']

PLOTTING=config['plotting']

DPI=PLOTTING.getint('dpi')

LABELS=config['labels']

COUPLING=config['coupling']

SAVGOL=config['savgol']

BOUNDS=config['fitting']

