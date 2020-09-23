import os
import pandas as pd
import numpy as np

from ..config import PATHS

def samplelog(info=None,overwrite=False):
    path=os.path.join(PATHS['dir_home'],'Samplelog.xlsx')
    
    try:
        samplelog=pd.read_excel(path,sep=';',decimal=',',index_col=0)
    except:
        print('> Samplelog konnte nicht gefunden werden. Leerer Samplelog wird angelegt.')
        samplelog=pd.DataFrame(columns=['alias','reference'])
        
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
    
        samplelog.to_excel(os.path.join(PATHS['dir_home'],'Samplelog.xlsx'))
    return samplelog