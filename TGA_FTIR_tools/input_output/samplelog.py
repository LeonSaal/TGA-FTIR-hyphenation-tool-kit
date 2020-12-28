import os
import pandas as pd
import numpy as np

from ..config import PATHS

def samplelog(info=None,overwrite=False):
    path=os.path.join(PATHS['dir_home'],'Samplelog.xlsx')
    
    try:
        samplelog=pd.read_excel(path,index_col=0)
    except:
        print('> \'Samplelog.xlsx\' was not found. New file was created under \'{}\'.'.format(path))
        samplelog=pd.DataFrame(columns=['alias','reference'])
        samplelog.to_excel(path)
        
    if info!=None:
        name=info['name']
        data=pd.DataFrame.from_dict(info,orient='index',columns=[name]).T.drop(['name'],1)
        
        for key in data.columns:
            if key not in samplelog.columns:
                samplelog[key]=np.nan
        if name in samplelog.index:
            if overwrite==False:
                samplelog=samplelog.fillna(data)
            else:
                samplelog.loc[[name]]=data
        else:
            samplelog=samplelog.append(data)
    
        try:
            samplelog.to_excel(path)
            print('Successfully updated \'Samplelog.xlsx\'.')
        except:
            print('Unable to write on \'Samplelog.xlsx\'. Please close file and try again!')
        
    return samplelog
