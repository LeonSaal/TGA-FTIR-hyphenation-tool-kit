import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import copy
import os
import re
import pickle

from .config import SAVGOL, PATHS
from .input_output import corrections, TGA, FTIR, samplelog, general
from .calibration import calibrate
from .plotting import plot_FTIR, plot_TGA, FTIR_to_DTG
from .fitting import fitting as fit

WINDOW_LENGTH=SAVGOL.getint('window_length')
POLYORDER=SAVGOL.getint('POLYORDER')


class TG_IR:
    try:
        linreg,stats = calibrate(mode='load')
    except:
        pass

    def __init__(self,name,mode='construct',profile='Otto',alias='load'):
        if mode=='construct':
            try:
                self.tga=TGA.read_TGA(name,profile=profile)
                self.info=TGA.TGA_info(name,self.tga,profile=profile)
                self.tga['dtg']=-savgol_filter(self.tga['mass'],WINDOW_LENGTH,POLYORDER,deriv=1)
            except:
                print('No TG data for {} was found'.format(name))
            
            try:
                self.ir=FTIR.read_FTIR(name)
                self.info['gases']=self.ir.columns[1:].to_list()
            except:
                print('No IR data for {} was found'.format(name))
            self.info.update(FTIR.FTIR_info(self))
            try:
                self.ir['t']=self.ir['t']+60*self.info['background_delay']
                self.ir=pd.merge(self.tga.filter(['t','Ts'],axis=1),self.ir, how='left',on='t').dropna(axis=0)
            except:
                pass
            if alias=='load':
                try:
                    alias=samplelog().loc[name,'alias']
                except:
                    alias=np.nan
                if type(alias)!=str:
                    self.info['alias']=self.info['name']
                else:
                    self.info['alias']=alias
                    
            elif type(alias)==str:
                self.info['alias']=alias
            else:
                self.info['alias']=self.info['name']
        if mode=='pickle':
            with open(os.path.join(PATHS['dir_output'],name+'.pkl'), 'rb') as inp:
                obj = pickle.load(inp)
            for key in obj.__dict__: 
                self.__dict__[key]=obj.__dict__[key]
               
        
    def corr(self,reference,**kwargs):
        if 'reference' in self.info:
            print('Sample has already been corrected! Re-initialise object for correction.')
            return
        
        # correction of data
        self.info['reference']=reference
        try:
            self.tga=corrections.corr_TGA(self.tga,reference)
            self.tga['dtg']=-savgol_filter(self.tga['mass'],WINDOW_LENGTH,POLYORDER,deriv=1)
        except:
            print('Failed to correct TG data.')
            
        try:
            self.ir=corrections.corr_FTIR(self.ir,reference)
        except:
            print('Failed to correct IR data.')
            
        # filling TG_IR.info
        try:
            TGA.dry_weight(self,**kwargs)
        except:
            print('Failed to derive TG info.')
        try:
            self.info.update(FTIR.FTIR_info(self))
        except:
            print('Failed to derive IR info.')
                
    def plot(self,which,**kwargs):
        if ('ir' not in self.__dict__) and (which in ['IR','DIR','cumsum','IR_to_DTG']):
            print('Option unavailable without IR data.')
            return
        else:
            if which=='IR':
                plot_FTIR(self,**kwargs)
            
            if which=='DIR':
                temp=copy.deepcopy(self)
                temp.ir.update(self.ir.filter(self.info['gases'],axis=1).diff().ewm(span = 10).mean())
                plot_FTIR(temp,**kwargs)
                
            if which=='cumsum':
                temp=copy.deepcopy(self)
                temp.ir.update(self.ir.filter(self.info['gases'],axis=1).cumsum())
                plot_FTIR(temp,**kwargs)
                
        if ('tga' not in self.__dict__) and (which in ['TG','heat_flow','IR_to_DTG']):
            print('Option unavailable without TGA data.')
            return
        else:            
            if which=='TG':
                plot_TGA(self,'mass',**kwargs)
            
            if which=='heat_flow':
                if 'heat_flow' in self.tga.columns:
                    plot_TGA(self,which,**kwargs)
                else:
                    print('No heat flow data available!')
        
        if ('linreg' in TG_IR.__dict__) and (which == 'IR_to_DTG'):
            FTIR_to_DTG(self,**kwargs)
        elif ('linreg' not in TG_IR.__dict__) and (which == 'IR_to_DTG'):
            print('Option unavailable without calibration!')
            return
        
        
    def fit(self,reference,tol_c=30,tol_hwhm=95,T_max=None,plot=True,func=fit.multi_gauss,y_axis='orig',save=True):
        if T_max==None:
            T_max=max(self.tga['Ts'])
        elif T_max>max(self.tga['Ts']):
            print('$T_{max}$ exceeds maximum temperature of data')
            T_max=max(self.tga['Ts'])
            
        ###extracting initial values for fitting from reference
        references=pd.read_excel(os.path.join(PATHS['dir_home'],'Initial_center.xlsx'),index_col=0,header=None)
        init_params=references.loc[['gas',reference]].dropna(axis=1)

        #find needed gases
        gases=list(set(init_params.loc['gas']))

        #initial center
        temps=pd.DataFrame(columns=gases,index=range(len(init_params.columns)))
        
        #labels for the center
        labels=pd.DataFrame(columns=gases,index=range(len(init_params.columns)))
        for column in init_params.columns:
            temps[init_params[column]['gas']][temps[init_params[column]['gas']].count()]=init_params[column][reference]
            labels[init_params[column]['gas']][labels[init_params[column]['gas']].count()]=references[column]['group']
        temps=temps.where(temps<T_max+tol_hwhm).dropna(thresh=1)
        labels=labels.where(temps<T_max+tol_hwhm).dropna(thresh=1)
        
        if save==True:
            path=os.path.join(PATHS['dir_fitting'],general.time()+reference+'_'+self.info['name'])
            os.makedirs(path)
            os.chdir(path)
            
        temp=copy.deepcopy(self)
        temp.tga=temp.tga[temp.tga['Ts']<T_max]
        temp.ir=temp.ir[temp.ir['Ts']<T_max]
        peaks, sumsqerr=fit.fitting(temp,temps,labels,func,tol_center=tol_c,max_hwhm=tol_hwhm,plot=plot,y_axis=y_axis,save=save)
        
        os.chdir(PATHS['dir_home'])
        return peaks, sumsqerr
    
    def save(self,how='pickle',**kwargs):
        samplelog(self.info,**kwargs)
        path_output=PATHS['dir_output']
        if os.path.exists(path_output)==False:
            os.makedirs(path_output)
        if how=='pickle':
            with open(os.path.join(path_output,self.info['name']+'.pkl'),'wb') as output:
                pickle.dump(self,output,pickle.HIGHEST_PROTOCOL)
        elif how=='excel':
            with pd.ExcelWriter(os.path.join(path_output,self.info['name']+'.xlsx')) as writer:
                try:
                    pd.DataFrame.from_dict(self.info,orient='index').to_excel(writer,sheet_name='info')
                except:
                    pass
                for key in ['tga','ir']:
                    try:
                        self.__dict__[key].to_excel(writer,sheet_name=key)
                    except:
                        pass
    def calibrate(self,**kwargs):
        try:
            self.linreg,self.stats=calibrate(**kwargs)
        except:
            pass
                

def fits(*TG_IR,reference=None,tol_c=30,tol_hwhm=95,T_max=None,plot=True,y_axis='orig',init_offs=[0,0,0],save=True,labels=None,temps=None):
    if reference!=None:
        references=pd.read_excel(os.path.join(PATHS['dir_home'],'Initial_center.xlsx'),index_col=0,header=None)
        if T_max==None:
            T_max=min([max(obj.tga['Ts']) for obj in TG_IR])

        elif T_max>min([max(obj.tga['Ts']) for obj in TG_IR]):
            T_max=T_max=min([max(obj.tga['Ts']) for obj in TG_IR])
            print('$T_{max}$ exceeds maximum temperature of data')


        init_params=references.loc[['gas',reference]].dropna(axis=1)
        #find needed gases
        gases=list(set(init_params.loc['gas']))

        ###extracting initial values for fitting from reference
        #initial center
        temps=pd.DataFrame(columns=gases,index=range(len(init_params.columns)))
        #labels for the center
        labels=pd.DataFrame(columns=gases,index=range(len(init_params.columns)))
        for column in init_params.columns:
            temps[init_params[column]['gas']][temps[init_params[column]['gas']].count()]=init_params[column][reference]
            labels[init_params[column]['gas']][labels[init_params[column]['gas']].count()]=references[column]['group'] 
        temps=temps.where(temps<T_max+tol_hwhm).dropna(thresh=1)
        labels=labels.where(temps<T_max+tol_hwhm).dropna(thresh=1)

    else:
        gases=labels.columns
    
    #initializing of output DataFrames
    col_labels=[group+'_'+gas.upper() for gas in gases for group in labels[gas].dropna()]+[gas.upper() for gas in gases]
    err=pd.DataFrame(columns=gases)
    names=['center','height','hwhm','area','mmol','mmol_per_mg']
    res=dict()
    for name in names:
        res[name]=pd.DataFrame(columns=col_labels)
    
    #make subdirectory to save data
    if save==True:
        path=os.path.join(PATHS['dir_fitting'],general.time()+reference+'_'+'_'.join(list(set([str(obj.info['sample']) for obj in TG_IR]))))
        os.makedirs(path)
        os.chdir(path)
    
    #cycling through samples
    for obj in TG_IR:
        #fitting of the sample and calculating the amount of functional groups
        if T_max!=None:
            temp=copy.deepcopy(obj)
            temp.tga=temp.tga[temp.tga['Ts']<T_max]
            temp.ir=temp.ir[temp.ir['Ts']<T_max]
        else:
            temp=obj
        peaks,sumsqerr=fit.fitting(temp,temps,labels,fit.multi_gauss,tol_center=tol_c,max_hwhm=tol_hwhm,plot=plot,y_axis=y_axis,save=save,init_offs=init_offs)

        #writing data to output DataFrames
        for key in res:
            res[key]=res[key].append(peaks[key].rename(obj.info['name']).T)  
        err=err.append(sumsqerr)

    # calculate statistical values
    dm=1e-6
    for key in res:
        samples=list(set([plt.get_label(re.search('(?<=_)\d{5}(?=_\d{2,3})',index).group()) for index in res[key].index]))
        stddev=pd.DataFrame(columns=res[key].columns,index=[sample+'_stddev' for sample in samples])
        mean=pd.DataFrame(columns=res[key].columns,index=[sample+'_mean' for sample in samples])
        for sample in samples:
            for column in res[key].columns:
                gas=column[column.rfind('_')+1:].lower()
                indices=[index for index in res[key].index if plt.get_label(re.search('(?<=_)\d{5}(?=_\d{2,3})',index).group())==sample]
                subset=res[key][column].loc[indices]
                if key=='mmol_per_mg':
                    mmol=res['mmol'][gas.upper()].loc[indices]#res['mmol'][column].loc[indices]
                    g=mmol/subset
                    lod=TG_IR[0].stats['x_NG'][gas]
                    dmmolg_i=np.power(np.power(lod/mmol,2)+np.power(dm/g,2),0.5)*subset
                    dmmol=np.power(np.sum(np.power(dmmolg_i,2)),0.5)
                    stddev[column][sample+'_stddev']=dmmol
                else:
                    stddev[column][sample+'_stddev']=np.std(subset)
                mean[column][sample+'_mean']=np.mean(subset)
        res[key]=res[key].append(mean)
        res[key]=res[key].append(stddev)        
    
    #exporting data
    if save==True:
        with pd.ExcelWriter('summary.xlsx') as writer:
            for key in res:
                res[key].dropna(axis=1).to_excel(writer,sheet_name=key)
            err.to_excel(writer,sheet_name='sum_squerr')
        os.chdir(PATHS['dir_home'])
    return res