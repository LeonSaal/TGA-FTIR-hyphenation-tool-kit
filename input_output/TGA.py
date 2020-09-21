import re
import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from .general import find_files
from ..config import PATHS, COUPLING, DPI
from ..plotting import get_label



def read_TGA(file,profile='Otto'):
    #open file from TGA in given directory and make a DataFrame from it
    path=find_files(file,'.txt',PATHS['dir_data'])[0]
    
    if profile=='Otto':
        skiprows=13
    elif profile=='Falk':
        skiprows=11
        
    try:
        data=pd.read_csv(path, delim_whitespace=True,decimal=',' ,names=['Index','t','Ts','Tr','mass'],skiprows=skiprows, skipfooter=11,converters={'mass':lambda x: float(x)},engine='python').drop(columns='Index')
        
    except:
        return
    
    #check if there is heat flow information and append it 
    try:
        path_mW=find_files(file,'_mW.txt',PATHS['dir_data'])[0]
        data['heat_flow']=pd.read_csv(path_mW, delim_whitespace=True,decimal=',' ,names=['Index','t','Ts','Tr','heat_flow'],skiprows=skiprows, skipfooter=11, usecols=['heat_flow'],engine='python')
    except:
        pass
    
    return data


def TGA_info(file,TGA,profile='Otto'):
    #open file from TGA in given directory and make a DataFrame from it
    path=find_files(file,'.txt',PATHS['dir_data'])[0]
    
    if profile=='Otto':
        skipheader=7
        
    if profile=='Falk':
        skipheader=6
        
    #extract information on the measurement from the header and footer of the TGA file
    footer=pd.read_table(path,encoding='ansi',skipfooter=2,index_col=False,names=[0],engine='python').tail(3)
    header=str(pd.read_table(path,encoding='ansi',skiprows=skipheader,nrows=1,index_col=False,names=[0],engine='python').iloc[0,0])
    info={}
    info['name']=re.search('\S+(?=,)',header).group()
    info['date']=re.search('(?<=,\s)\S+(?=\s)',header).group()#header[header.find(',')+2:header.rfind(' ')].strip()
    info['time']=re.search('(?<=\s)\S+$',header).group()#header[header.rfind(' ')+1:].strip()
    method_name=str(footer.iloc[2,0]).strip()
    info['method_name']=method_name
    info['sample_mass']=pd.to_numeric(re.search('(?<=\s)\S+(?=\smg)',footer.iloc[0,0]).group().replace(',','.'))#str(footer).strip()
    #if the sample wasn't weighed in automatically, the mass at t=0 is used instead
    if info['sample_mass']==0:
        info['sample_mass']=TGA['mass'][0]
    info['reference_mass']='sample_mass'

    # extract method info from method
    last_i=0
    values=[]
    parameters=[]
    for i in range(len(method_name)):
        if method_name[i]=='=':
            parameters.append('background_state')
            
        elif method_name[i]=='<':
            parameters.append('lower_temp')
            
        elif (method_name[i]=='>') or (method_name[i]=='/'):
            parameters.append('high_temp')
        
        elif (method_name[i]==')') and (method_name[last_i-1]=='('):
            parameters.append('method_gas')
            
        elif (method_name[i]=='(') and (method_name[last_i-1]=='/'):
            parameters.append('gradient')
            
        elif (method_name[i]=='/') and (method_name[last_i-1]=='<'):
            parameters.append('high_temp')
            
        if (method_name[i] in '=<>()/')==True:
            val=method_name[last_i:i]
            if val.isnumeric()==True:
                val=int(val)
            values.append(val)
            last_i=i+1
    
    parameters.append('crucible')
    values.append(method_name[method_name.rfind('_')+1:])
    try:
        info['background_delay']=int(re.search('^\d+(?==)',method_name).group())
    except:
        info['background_delay']=COUPLING.getint('background_delay')
    info['method_gases']=[values[index].lower() for index in range(len(parameters)) if parameters[index]=='method_gas']
    info['switch_temp']=[values[index] for index in range(len(parameters)) if parameters[index]=='high_temp']
    
    return info

def dry_weight(TG_IR,how='h2o',plot=False,ref_mass=None,save=False):
    if type(how)==type(None):
        dry_point=0
    elif type(how)!=str:
        dry_point=TG_IR.tga['t'][TG_IR.tga['Ts']>=how].iloc[0]
    else :
        if how=='h2o':
            try:
                ref=TG_IR.ir.filter(items=['Ts','h2o'])
                #ref['h2o']/=TG_IR.linreg['slope']['h2o']*18/15/TG_IR.info['sample_mass']
                ylabel=get_label(how)
            except:
                print('No water signal found. Falling back to DTG.')
                how='mass'
        
        if how=='mass':
            ref=TG_IR.tga.filter(items=['Ts','mass'])
            ref['mass']=-sp.signal.savgol_filter(TG_IR.tga['mass']/TG_IR.info['sample_mass'], 51, 3, deriv=1)
            ylabel='DTG'
        min_T=ref['Ts'][ref[how]>=max(ref[how][(ref['Ts']>50) & (ref['Ts']<200)])].values[0]
        max_T=min_T+50

        x=ref['Ts'][(ref['Ts']>min_T) & (ref['Ts']<max_T)]
        y=ref[how][(ref['Ts']>min_T) & (ref['Ts']<max_T)]
        slope,intercept,r_val,p_val,std_err =sp.stats.linregress(x,y)

        dry_point=ref['Ts'][ref['Ts'] >=-intercept/slope].index[0]
    
    #getting the dry_mass at the dry_point as well as the final weight and calculating the relative
    #mass-loss and the water content from it
    dry_temp=TG_IR.tga['Ts'][dry_point]
    names=['dry']+TG_IR.info['method_gases']
    info={}
    if (how=='h2o') or (how=='mass'):
        times=[0,dry_point]+list(TG_IR.tga.index[TG_IR.tga['Tr'].isin(TG_IR.info['switch_temp'])])
        names=['dry']+TG_IR.info['method_gases']
        info['dry_mass']=TG_IR.tga['mass'][dry_point]
        info['reference_mass']='dry_mass'
        info['dry_temp']=dry_temp
        info['dry_time']=dry_point
    elif how==None:
        info['reference_mass']='sample_mass'
        times=[0]+list(TG_IR.tga.index[TG_IR.tga['Tr'].isin(TG_IR.info['switch_temp'])])
        names=TG_IR.info['method_gases']   
    if ref_mass!=None:
        info['reference_mass']=ref_mass
    weights=TG_IR.tga['mass'][TG_IR.tga.index.isin(times)].values
    mass_loss=abs(np.diff(weights))
    
    info['final_mass']=TG_IR.tga['mass'][len(TG_IR.tga)-1]
    TG_IR.info.update(info)
    for name,ml in zip(names,mass_loss):
        info['ml_'+name]=ml
        info['rel_ml_'+name]=ml/TG_IR.info[TG_IR.info['reference_mass']]
    TG_IR.info.update(info)
    
    #plotting
    if plot==True:
        fig=plt.figure()
        x=TG_IR.tga['Ts']
        y=TG_IR.tga['mass']
        plt.plot(x,y,label='TGA')
        
        for i in range(len(times)-1):
            plt.annotate(s='', xy=(x[times[i]],y[times[i]]), xytext=(x[times[i]],y[times[i+1]]), arrowprops=dict(arrowstyle='<->'))
            plt.text(x[times[i]]+20,y[times[i]],'{:.2f} mg @ {:.2f} °C'.format(weights[i],x[times[i]]))
            plt.text(x[times[i]]+20,(y[times[i]]+ y[times[i+1]])/2,'$ML$ {}: {:.2f} mg ({:.1f} %)'.format(get_label(names[i]),mass_loss[i],mass_loss[i]/TG_IR.info[TG_IR.info['reference_mass']]*100))
        
        plt.scatter(x[times],y[times],c='r')
        plt.hlines(weights[1:],x[times[:-1]],x[times[1:]],linestyle='dashed')
        plt.ylabel('TGA /mg')
        plt.xlabel('T /°C')
        if type(how)==str:
            ax2=plt.twinx()
            ax2.plot(ref['Ts'],ref[how],linestyle='dashed')
            ax2.set_ylabel(ylabel)
        plt.show()   
                
        if save:
            path_plots=PATHS['dir_plots']
            if os.path.exists(path_plots)==False:
                os.makedirs(path_plots)
            fig.savefig(os.path.join(path_plots,'{}_mass_steps.png'.format(TG_IR.info['name'])), bbox_inches='tight',dpi=DPI)
        
