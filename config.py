import configparser
config = configparser.ConfigParser()
config.read('settings.ini')
PATHS=config['paths']
UNITS=config['units']
SEP=UNITS['sep']
PARAMS=config['parameters']
MOLAR_MASS=config['molar_mass']
DPI=config['plotting'].getint('dpi')
LABELS=config['labels']
COUPLING=config['coupling']
SAVGOL=config['savgol']
BOUNDS=config['fitting']