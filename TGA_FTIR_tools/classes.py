import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import copy
import os

import pickle
from .config import SAVGOL, PATHS, COUPLING
from .input_output import corrections, TGA, FTIR, general, samplelog
from .calibration import calibrate
from .plotting import plot_FTIR, plot_TGA, FTIR_to_DTG
from .fitting import fitting, get_presets 

WINDOW_LENGTH=SAVGOL.getint('window_length')
POLYORDER=SAVGOL.getint('POLYORDER')


class TG_IR:
    def __init__(self, name, mode='construct', profile=COUPLING['profile'], alias='load',**kwargs):
        if mode=='construct':
            try:
                # load TG data
                self.tga=TGA.read_TGA(name,profile=profile)
                self.tga['dtg']=-savgol_filter(self.tga['sample_mass'],WINDOW_LENGTH,POLYORDER,deriv=1)
                
                # deriving TG info
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
                print('TGA data found.')
        
            except:
                del self.tga
        
            try:
                # load IR data
                self.ir=FTIR.read_FTIR(name)
                self.info['gases']=self.ir.columns[1:].to_list()
                try:
                    # load calibration
                    self.linreg,self.stats = calibrate(mode='load')
                except:
                    pass
                print('IR data found{} for gases {}.'.format(' and calibrated (*)' if 'linreg' in self.__dict__ else '',', '.join([gas+('*' if 'linreg' in self.__dict__ and gas in self.linreg.index else '') for gas in self.info['gases']])))
                #derivinf IR info
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
                print('No IR data found.')
            
            if 'ir' not in self.__dict__ and 'tga' not in self.__dict__:
                return
            else:
                print('>> \'{}\' successfully initialiazed.'.format(name))
            
            # assigning alias
            if alias=='load':
                try:
                    # from samplelog
                    alias=samplelog(create=False).loc[name,'alias']
                except:
                    alias=np.nan
                if type(alias)!=str:
                    # default name
                    self.info['alias']=name
                else:
                    self.info['alias']=alias
                    
            elif type(alias)==str:
                self.info['alias']=alias
            else:
                self.info['alias']=self.info['name']
        
        # initialize object from pickle file
        if mode=='pickle':
            with open(os.path.join(PATHS['dir_output'],name+'.pkl'), 'rb') as inp:
                obj = pickle.load(inp)
            for key in obj.__dict__: 
                self.__dict__[key]=obj.__dict__[key]
               
        
    def corr(self,reference='load',plot=False,**kwargs):
        "correction of TG and IR data"
        
        if 'reference' in self.info:
            print('Sample has already been corrected! Re-initialise object for correction.')
            return
        
        if reference=='load':
            # try to load reference from samplelog if none is supplied
            try:
                reference=samplelog(create=False).loc[self.info['name'],'reference']
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
        
        if hasattr(self, 'ir'):
            try:
                self.ir.update(corrections.corr_FTIR(self.ir,reference,plot=plot))
            except:
                print('Failed to correct IR data.')
            
        # filling TG_IR.info
        try:
            TGA.dry_weight(self,plot=plot,**kwargs)
            print('\'TG_IR.info\' of {} was updated. To store these in Samplelog.xlsx run \'TG_IR.save()\''.format(self.info['name']))
            success=True
        except:
            print('Failed to derive TG info.')

        if hasattr(self, 'ir'):
            try:
                self.info.update(FTIR.FTIR_info(self))
                if not success:
                    print('\'TG_IR.info\' was updated. To store these in Samplelog.xlsx run \'TG_IR.save()\'')
            except:
                print('Failed to derive IR info.')
            
    def get_value(self, *values, which='sample_mass', at='sample_temp'):
        "extract values from TG data at e.g. certain temperatures"
        
        out = pd.DataFrame(index=[which],columns=pd.Index(values,name=at))
        for value in values:
            out.loc[which,value]=self.tga[which][self.tga[at]>=value].values[0]
    
        return out
    
    def dry_weight(self, **kwargs):
        "determine dry point and mass of sample"
        
        try:
            TGA.dry_weight(self,**kwargs)
            print('\'TG_IR.info\' was updated. To store these in Samplelog.xlsx run \'TG_IR.save()\'')

        except:
            print('Failed to derive TG info.')
                
    def plot(self, plot, **kwargs):
        "plotting TG and or IR data"
        
        options=['TG', 'heat_flow', 'IR', 'DIR', 'cumsum', 'IR_to_DTG']
        if plot not in options:
            print('\'TG_IR.plot\' supports {} as input to plot a respective figure.'.format(', '.join(['\''+option+'\'' for option in options])))
            
        if ('ir' not in self.__dict__) and (plot in ['IR','DIR','cumsum','IR_to_DTG']):
            print('Option unavailable without IR data.')
            return
        else:
            if plot == 'IR':
                plot_FTIR(self,**kwargs)
            
            if plot == 'DIR':
                temp=copy.deepcopy(self)
                temp.ir.update(self.ir.filter(self.info['gases'],axis=1).diff().ewm(span = 10).mean())
                plot_FTIR(temp,**kwargs)
                
            if plot == 'cumsum':
                temp=copy.deepcopy(self)
                temp.ir.update(self.ir.filter(self.info['gases'],axis=1).cumsum())
                plot_FTIR(temp,**kwargs)
                
        if ('tga' not in self.__dict__) and (plot in ['TG','heat_flow','IR_to_DTG']):
            print('Option unavailable without TGA data.')
            return
        else:            
            if plot == 'TG':
                plot_TGA(self,'sample_mass',**kwargs)
            
            if plot == 'heat_flow':
                if 'heat_flow' in self.tga.columns:
                    plot_TGA(self, plot, **kwargs)
                else:
                    print('No heat flow data available!')
        
        if ('linreg' in self.__dict__) and (plot == 'IR_to_DTG'):
            FTIR_to_DTG(self,**kwargs)
        elif ('linreg' not in self.__dict__) and (plot == 'IR_to_DTG'):
            print('Option unavailable without calibration!')
            return
        
        
    def fit(self,reference,T_max=None,T_max_tol=50,save=True,plot=True,presets=None,**kwargs):
        "deconvolution of IR data"
        if 'linreg' not in self.__dict__:
            print('Option unavailable without calibration!')
            return
            
        # setting upper limit for data 
        if T_max==None:
            T_max=max(self.tga['sample_temp'])
        elif T_max>max(self.tga['sample_temp']):
            print('T_max exceeds maximum temperature of data')
            T_max=max(self.tga['sample_temp'])
        
        # load presets for deconvolution
        if presets==None:
            presets=get_presets(PATHS['dir_home'], reference)
            
        for gas in presets:
            presets[gas]=presets[gas].drop(presets[gas].index[presets[gas].loc[:,'center_0']>T_max+T_max_tol])
        
        # setting up output directory
        if save:
            path=os.path.join(PATHS['dir_fitting'],general.time()+reference+'_'+self.info['name']).replace(os.sep,os.altsep)
            os.makedirs(path)
            os.chdir(path)
        
        # fitting
        print('Fitting according to \'{}\' in Fitting_parameters.xlsx is in progress ...'.format(reference))
        temp=copy.deepcopy(self)
        temp.tga=temp.tga[temp.tga['sample_temp']<T_max]
        temp.ir=temp.ir[temp.ir['sample_temp']<T_max]
        peaks, sumsqerr=fitting(temp,presets,plot=plot,save=save,**kwargs)
        if save:
            print('Fitting finished! Plots and results are saved in \'{}\'.'.format(path))
            os.chdir(PATHS['dir_home'])
        return peaks, sumsqerr;
    
    
    def save(self,how='samplelog',**kwargs):
        "save object or its contents as pickle file or excel"
        
        # update samplelog
        samplelog(self.info, create=True, **kwargs)
        path_output=PATHS['dir_output']
        if how=='samplelog':
            return
        
        if os.path.exists(path_output)==False:
            os.makedirs(path_output)
        # save object
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
        "calibrate object"
        try:
            self.linreg,self.stats=calibrate(**kwargs)
        except:
            pass
