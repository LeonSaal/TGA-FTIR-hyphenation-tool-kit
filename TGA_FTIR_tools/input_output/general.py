import datetime as dt
import os
import re
import pandas as pd
    
def time():
    return str(dt.datetime.now().date())+'_'+str(dt.datetime.now().hour)+'-'+str(dt.datetime.now().minute).zfill(2)+'-'+str(dt.datetime.now().second).zfill(2)+'_'    

def find_path(file,parent_dir):
    for dirpath,dirnames,filenames in os.walk(parent_dir):
        for filename in filenames:
            if filename[:filename.rfind('.')].find(file)!=-1:
                return dirpath
            
def find_files(file,suffix,parent_dir):
    files=[]
    for dirpath,dirnames,filenames in os.walk(parent_dir):
        for filename in filenames:
            if (re.search(re.escape(file),filename)!=None) and (re.search(suffix+'$',filename,flags=re.IGNORECASE)):
                files.append(os.path.join(dirpath,filename))   
    return files

def overview(objs):
    out=pd.DataFrame(columns=['name','alias'])
    out.index.name='index'
    for attr in ['name','alias']:
        out[attr]=[obj.info[attr] for obj in objs]
    print(out)
    return
