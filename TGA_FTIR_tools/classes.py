import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import copy
import os

import pickle
from .config import SAVGOL, PATHS, COUPLING
from .input_output import corrections, TGA, FTIR, general
from .input_output import samplelog
from .calibration import calibrate
from .plotting import plot_FTIR, plot_TGA, FTIR_to_DTG
from .fitting import fitting, get_presets 

WINDOW_LENGTH=SAVGOL.getint('window_length')
POLYORDER=SAVGOL.getint('POLYORDER')


class TG_IR:
    def __init__(self,name,mode='construct',profile='Otto',alias='load',**kwargs):
        if mode=='construct':
            try:
                self.tga=TGA.read_TGA(name,profile=profile)
                self.tga['dtg']=-savgol_filter(self.tga['sample_mass'],WINDOW_LENGTH,POLYORDER,deriv=1)

                try:
                    self.info=TGA.TGA_info(name,self.tga,profile=profile)
                except:
                    print('Failed to derive TG info. Using default values.')
                    self.info=dict()
                    self.info['name']=name
                    self.info['initial_mass']=self.tga.loc[0,'sample_mass']
                    self.info['reference_mass']='initial_mass'
                    self.info['background_delay']=COUPLING.getint('background_delay')
                    self.info['switch_temp']=[max(self.tga['reference_temp'])]
                    self.info['method_gases']=['? method_gas ?']
                try:
                    TGA.dry_weight(self,**kwargs)
                except:
                    pass
                print('TGA data was found.')
        
            except:
                del self.tga
        
            try:
                self.ir=FTIR.read_FTIR(name)
                self.info['gases']=self.ir.columns[1:].to_list()
                try:
                    self.linreg,self.stats = calibrate(mode='load')
                except:
                    pass
                print('IR data found{} for gases {}.'.format(' and calibrated (*)' if 'linreg' in self.__dict__ else '',', '.join([gas+('*' if 'linreg' in self.__dict__ and gas in self.linreg.index else '') for gas in self.info['gases']])))
                try:
                    self.info.update(FTIR.FTIR_info(self))
                except:
                    pass
                try:
                    self.ir['time']=self.ir['time']+60*self.info['background_delay']
                    self.ir=pd.merge(self.tga.filter(['time','sample_temp','reference_temp'],axis=1),self.ir, how='left',on='time').dropna(axis=0)
                except:
                    pass
            except:
                del self.ir
                print('No IR data was found.')
            
            if 'ir' not in self.__dict__ and 'tga' not in self.__dict__:
                return
            else:
                print('>> \'{}\' successfully initialiazed.\n'.format(name))
            
            if alias=='load':
                try:
                    alias=samplelog().loc[name,'alias']
                except:
                    alias=np.nan
                if type(alias)!=str:
                    self.info['alias']=name
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
               
        
    def corr(self,reference='load',plot=False,**kwargs):
        if 'reference' in self.info:
            print('Sample has already been corrected! Re-initialise object for correction.')
            return
        
        if reference=='load':
                try:
                    reference=samplelog().loc[self.info['name'],'reference']
                except:
                    print('No reference found in Samplelog. Please supply \'reference = \'')
                    return

        # correction of data
        self.info['reference']=reference
        try:
            self.tga=corrections.corr_TGA(self.tga,reference,plot=plot)
            self.tga['dtg']=-savgol_filter(self.tga['sample_mass'],WINDOW_LENGTH,POLYORDER,deriv=1)
        except:
            print('Failed to correct TG data.')
            
        try:
            self.ir.update(corrections.corr_FTIR(self.ir,reference,plot=plot))
        except:
            print('Failed to correct IR data.')
            
        # filling TG_IR.info
        try:
            TGA.dry_weight(self,plot=plot,**kwargs)
            print('\'TG_IR.info\' was updated. To store these in Samplelog.xlsx run \'TG_IR.save()\'')
            success=True
        except:
            print('Failed to derive TG info.')

        try:
            self.info.update(FTIR.FTIR_info(self))
            if not success:
                print('\'TG_IR.info\' was updated. To store these in Samplelog.xlsx run \'TG_IR.save()\'')
        except:
            print('Failed to derive IR info.')
            
    def get_value(self,*values, which='sample_mass', at='sample_temp'):
        out = pd.DataFrame(index=[which],columns=pd.Index(values,name=at))
        for value in values:
            out.loc[which,value]=self.tga[which][self.tga[at]>=value].values[0]
    
        return out
    def dry_weight(self,**kwargs):
        try:
            TGA.dry_weight(self,**kwargs)
            print('\'TG_IR.info\' was updated. To store these in Samplelog.xlsx run \'TG_IR.save()\'')

        except:
            print('Failed to derive TG info.')
                
    def plot(self,which,**kwargs):
        options=['TG', 'heat_flow', 'IR', 'DIR', 'cumsum', 'IR_to_DTG']
        if which not in options:
            print('\'TG_IR.plot\' supports {} as input for \'which\' figure to plot.'.format(', '.join(['\''+option+'\'' for option in options])))
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
                plot_TGA(self,'sample_mass',**kwargs)
            
            if which=='heat_flow':
                if 'heat_flow' in self.tga.columns:
                    plot_TGA(self,which,**kwargs)
                else:
                    print('No heat flow data available!')
        
        if ('linreg' in self.__dict__) and (which == 'IR_to_DTG'):
            FTIR_to_DTG(self,**kwargs)
        elif ('linreg' not in self.__dict__) and (which == 'IR_to_DTG'):
            print('Option unavailable without calibration!')
            return
        
        
    def fit(self,reference,T_max=None,save=True,plot=True,presets=None,**kwargs):
        if T_max==None:
            T_max=max(self.tga['sample_temp'])
        elif T_max>max(self.tga['sample_temp']):
            print('$T_{max}$ exceeds maximum temperature of data')
            T_max=max(self.tga['sample_temp'])
        
        if presets==None:
            presets=get_presets(PATHS['dir_home'], reference)

        if save:
            path=os.path.join(PATHS['dir_fitting'],general.time()+reference+'_'+self.info['name']).replace(os.sep,os.altsep)
            os.makedirs(path)
            os.chdir(path)
            
        print('Fitting according to \'{}\' in Fitting_parameters.xlsx is in progress ...'.format(reference))
        temp=copy.deepcopy(self)
        temp.tga=temp.tga[temp.tga['sample_temp']<T_max]
        temp.ir=temp.ir[temp.ir['sample_temp']<T_max]
        peaks, sumsqerr=fitting(temp,presets,plot=plot,save=save,**kwargs)
        if save:
            print('Fitting finished! Plots and results are saved in \'{}\'.'.format(path))
            os.chdir(PATHS['dir_home'])
        return peaks, sumsqerr;
    
    
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
