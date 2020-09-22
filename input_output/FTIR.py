import re
import pandas as pd
import numpy as np

from .general import find_files
from ..config import PATHS


def read_FTIR(file_name):
    files=find_files(file_name,'.csv',PATHS['dir_data'])
    #find the gases by looking at the suffix of the files
    gases=[]
    paths=[]
    for file in files:
        gases.append(file[file.rfind('_')+1:file.rfind('.')].lower())
        paths.append(file) 

    if gases==[]:
        return
    else:
        #make DataFrame with the first gas, keeping the time column
        data=pd.read_csv(files[0], delimiter=';', decimal=',', names=['time',gases[0]])
    
        #append the IR data from the other gases as new columns
        for i in range(1,len(gases)):
            data[gases[i]]=pd.read_csv(files[i], delimiter=';', decimal=',', names=['time',gases[i]], usecols=[gases[i]])
        
        #convert time column for minutes to seconds 
        data['time']=(data['time']*60).astype(int)
        return data.dropna()

def FTIR_info(TG_IR):
    info={}
    for gas in TG_IR.info['gases']:
        info['area_{}'.format(gas)]=np.sum(TG_IR.ir[gas])
    for gas in TG_IR.info['gases']:
        try:
            info['mmol_{}'.format(gas)]=(info['area_{}'.format(gas)]-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]
        except:  
            pass
    try:
        elems=list(set(re.sub('\d', '',''.join([gas for gas in TG_IR.info['gases'] if len(gas)<5]))))
        for elem in elems:
            info['mmol_{}'.format(elem)]=0
            for gas in TG_IR.linreg.index:
                if elem in gas:
                    if re.search('(?<='+elem+')\d',gas) != None:
                        info['mmol_{}'.format(elem)]+=int(re.search('(?<='+elem+')\d',gas).group())*info['mmol_{}'.format(gas)]
                    else:
                        info['mmol_{}'.format(elem)]+=info['mmol_{}'.format(gas)]
    except:
        pass
     
    return info