import re
import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from .general import find_files
from ..config import PATHS, COUPLING, DPI, PARAMS, UNITS, SEP, SAVGOL
from ..plotting import get_label

WINDOW_LENGTH=SAVGOL.getint('window_length')
POLYORDER=SAVGOL.getint('POLYORDER')

def read_profiles(profile):
    profiles=pd.read_excel(os.path.join(PATHS['dir_home'],'TGA_import_profiles.xlsx'),index_col=0)
    
    return profiles.loc[profile,:].to_dict() 
    
def read_TGA(file,profile=COUPLING['profile']):
    "load TG data from file"
    #open file from TGA in given directory and make a DataFrame from it
    try:
        path=find_files(file,'.txt',PATHS['dir_data'])[0]
    except:
        print('No TG data for {} was found'.format(file))
        return
    
    profiles=read_profiles(profile)
    names=profiles['names'].split(',')
    names_heatflow=profiles['names_heatflow']
    skipfooter=profiles['skipfooter']
    skiprows=profiles['skiprows']
    sep=profiles['sep']
    drop=list(profiles['drop'])
        
    try:
        data=pd.read_csv(path, delim_whitespace=True,decimal=sep ,names=names,skiprows=skiprows, skipfooter=skipfooter,
                         converters={'sample_mass':lambda x: float(x.replace(sep,'.'))},engine='python').drop(columns=drop,errors='ignore')
        
    except:
        print('Failed to read TG-data from {}'.format(file))
        return
    
    #check if there is heat flow information and append it 
    try:
        path_mW=find_files(file,'_mW.txt',PATHS['dir_data'])[0]
        data['heat_flow']=pd.read_csv(path_mW, delim_whitespace=True,decimal=sep ,names=names_heatflow,skiprows=skiprows, skipfooter=skipfooter, converters={'sample_mass': lambda x: float(x.replace(sep,'.'))}, usecols=['heat_flow'],engine='python')
    except:
        pass
    
    return data


def TGA_info(file,TGA,profile='Otto'):
    "extract TG info e.g. measurement time, initial mass... from TG file"
    # open file from TGA in given directory and make a DataFrame from it
    path=find_files(file,'.txt',PATHS['dir_data'])[0]
    
    if profile=='Otto':
        skipheader=7
        
    if profile=='Falk':
        skipheader=6

    # extract information on the measurement from the header and footer of the TGA file
    footer = pd.read_table(path, encoding='ansi', skipfooter=2, index_col=False, names=[0], engine='python').tail(3)
    header = str(pd.read_table(path, encoding='ansi', skiprows=skipheader, nrows=1, index_col=False, names=[0], engine='python').iloc[0,0])
    info={}
    info['name']=file
    info['date']=re.search('\d{2}\.\d{2}\.\d{4}',header).group()
    info['time']=re.search('\d{2}:\d{2}:\d{2}',header).group()
    method=str(footer.iloc[2,0]).strip()
    info['method'] = method
    info['initial_mass']=pd.to_numeric(re.search('(?<=\s)\S+(?=\smg)',footer.iloc[0,0]).group().replace(',','.'))#str(footer).strip()
    
    # if the sample wasn't weighed in automatically, the mass at t=0 is used instead
    if info['initial_mass']==0:
        info['initial_mass']=TGA['sample_mass'][0]
    info['reference_mass']='initial_mass'
    
    # extract method info from method
    last_i=0
    values=[]
    parameters=[]
    for i in range(len(method)):
        if method[i]=='=':
            parameters.append('background_state')
            
        elif method[i]=='<':
            parameters.append('lower_temp')
            
        elif (method[i]=='>') or (method[i]=='('):
            parameters.append('high_temp')
            
        elif (method[i]=='/') and (method[last_i-1]=='<'):
            parameters.append('high_temp')  
        
        elif (method[i]==')') and (method[last_i-1]=='('):
            parameters.append('method_gas')
            
        elif (method[i]=='/'):
            parameters.append('gradient')

        if method[i] in '=<>()/_' and i-last_i>=1:
            val=method[last_i:i]
            if val.isnumeric()==True:
                val=int(val)
            values.append(val)
            last_i=i+1
            if method[i]=='_':
                break
        if i-last_i==0 and method[i] in '=<>()/_' and method[last_i] in '=<>()/_':
            last_i=i+1

    parameters.append('crucible')
    values.append(method[last_i:])
    try:
        info['background_delay']=int(re.search('^\d+(?==)',method).group())
    except:
        info['background_delay']=COUPLING.getint('background_delay')

    info['mass_steps']=[values[index].upper() for index in range(len(parameters)) if parameters[index]=='method_gas']
    info['step_temp']=[values[index] for index in range(len(parameters)) if parameters[index]=='high_temp']
    
    return info

def dry_weight(TG_IR, how_dry='H2O', plot=False, ref_mass='dry_mass', save=False, xlim=[None,None], ylim=[None,None]):
    "determine dry point and mass from TG data"
    # make h2o to H2O
    if how_dry=='h2o':
        how_dry=how_dry.upper()
        
    # if how_dry is None, no dry point is determined
    if how_dry==None:
        dry_point = 1
    
    # if how_dry is a number, the dry point is set to that temperature
    elif type(how_dry)!=str:
        # check if value is within the temperature range
        if (how_dry > max(TG_IR.tga['sample_temp'])) or (how_dry <= min(TG_IR.tga['sample_temp'])):
            how_dry = None
            dry_point = 1
            print('Supplied value is out of range. \'how_dry\' is set to None. Be aware that \'reference_mass\' is also set to \'initial_mass\'.')
        else:
            dry_point = TG_IR.tga['time'][TG_IR.tga['sample_temp']>=how_dry].values[0]
            # check if how_dry value could not be found (very scarce, but possible..)
            if (dry_point == 0):
                how_dry = None
                dry_point = 1
                print('Supplied value is out of range. \'how_dry\' is set to None. Be aware that \'reference_mass\' is also set to \'initial_mass\'.')
        
    # if how_dry is 'H2O' or 'sample_mass', the dry point is determined from the respective data
    else :
        if how_dry=='H2O':
            try:
                ref=TG_IR.ir.filter(items=['sample_temp','H2O'])
                ylabel=get_label(how_dry)
            except:
                how_dry='sample_mass'
        
        if how_dry=='sample_mass':
            ref=TG_IR.tga.filter(items=['sample_temp','sample_mass'])
            ref['sample_mass']=-sp.signal.savgol_filter(TG_IR.tga['sample_mass']/TG_IR.info['initial_mass'], WINDOW_LENGTH, POLYORDER, deriv=1)
            ylabel='DTG'
        min_T=ref['sample_temp'][ref[how_dry]>=max(ref[how_dry][(ref['sample_temp']>50) & (ref['sample_temp']<200)])].values[0]
        max_T=min_T+50

        x=ref['sample_temp'][(ref['sample_temp']>min_T) & (ref['sample_temp']<max_T)]
        y=ref[how_dry][(ref['sample_temp']>min_T) & (ref['sample_temp']<max_T)]
        slope,intercept,r_val,p_val,std_err =sp.stats.linregress(x,y)

        dry_point=ref['sample_temp'][ref['sample_temp'] >=-intercept/slope].index[0]
    
    # getting the dry_mass at the dry_point as well as the final weight and calculating the relative
    # mass-loss and the water content from it
    dry_temp=TG_IR.tga['sample_temp'][dry_point]
    info={}
    
    if how_dry == None:
        info['reference_mass']='initial_mass'
        times=[0]+list(TG_IR.tga.index[TG_IR.tga['reference_temp'].isin(TG_IR.info['step_temp'])])
        names=TG_IR.info['mass_steps']
        
    elif (how_dry == 'H2O') or (how_dry == 'sample_mass') or (type(how_dry) != str):   # turns to be true everytime
        times=[0,dry_point]+list(TG_IR.tga.index[TG_IR.tga['reference_temp'].isin(TG_IR.info['step_temp'])])
        names=['dry']+TG_IR.info['mass_steps']
        info['dry_mass']=TG_IR.tga['sample_mass'][dry_point]
        info['reference_mass']='dry_mass'
        info['dry_temp']=dry_temp
        info['dry_time']=dry_point
        
    if ref_mass!='dry_mass':
        info['reference_mass']=ref_mass
    weights=TG_IR.tga['sample_mass'][TG_IR.tga.index.isin(times)].values
    mass_loss=abs(np.diff(weights))
    
    info['final_mass']=TG_IR.tga['sample_mass'][len(TG_IR.tga)-1]
    TG_IR.info.update(info)
    for name,ml in zip(names,mass_loss):
        info['ml_'+name]=ml
        info['rel_ml_'+name]=ml/TG_IR.info[TG_IR.info['reference_mass']]
    TG_IR.info.update(info)
    
    # plotting
    if plot:
        fig, ax = plt.subplots()
        x=TG_IR.tga['sample_temp']
        y=TG_IR.tga['sample_mass']
        ax.plot(x,y,label='TGA')
        
        for i in range(len(times)-1):
            ax.annotate(text='', xy=(x[times[i]],y[times[i]]), xytext=(x[times[i]],y[times[i+1]]), arrowprops=dict(arrowstyle='<->'))
            #plt.text(x[times[i]]+20,y[times[i]],'{:.2f} mg @ {:.2f} °C'.format(weights[i],x[times[i]]))
            ax.text(x[times[i]]+20,(y[times[i]]+ y[times[i+1]])/2,'$ML$ {}: {:.2f} mg ({:.1f} %)'.format(get_label(names[i]),mass_loss[i],mass_loss[i]/TG_IR.info[TG_IR.info['reference_mass']]*100))
        
        ax.scatter(x[times],y[times],c='r')
        ax.hlines(weights[1:],x[times[:-1]],x[times[1:]],linestyle='dashed')
        ax.set_ylabel('{} {} {}'.format(PARAMS['sample_mass'],SEP,UNITS['sample_mass']))
        ax.set_xlabel('{} {} {}'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        ax.set_ylim(ylim)
        if type(how_dry)==str:
            ax2=plt.twinx()
            ax2.plot(ref['sample_temp'],ref[how_dry],linestyle='dashed',label=ylabel)
            ax2.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())  # switch on minor ticks on each axis
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        
        ax.set(title = 'Dry mass and mass steps determination')
        ax.legend()
        plt.show()   
                
        if save:
            path_plots=PATHS['dir_plots']
            if os.path.exists(path_plots)==False:
                os.makedirs(path_plots)
            fig.savefig(os.path.join(path_plots,'{}_mass_steps.png'.format(TG_IR.info['name'])), bbox_inches='tight',dpi=DPI)
        
