import pandas as pd
import numpy as np
from .fitting import fits, get_presets
from ..input_output.general import time
from ..plotting import get_label
import os
import re
import matplotlib.pyplot as plt
from ..config import PATHS, BOUNDS, UNITS, DPI
import copy

def robustness(*TG_IR,reference,T_max=None,save=True,var_T=10,var_rel=0.3,ylim=[0,None],**kwargs):
    presets=get_presets(PATHS['dir_home'], reference, TG_IR[0].ir)
    params=['center_0','tol_c','hwhm_max','height_0','hwhm_0']
    results=dict()
    variance=dict(zip(params,[var_T,var_T,var_T,var_rel,var_rel]))
    default=dict(zip(params,[0,BOUNDS.getfloat('tol_center'),BOUNDS.getfloat('hwhm_max'),BOUNDS.getfloat('height_0'),BOUNDS.getfloat('hwhm_0')]))
    
    #Init values
    print('Initial results:')
    res=fits(*TG_IR,reference=reference,save=False,T_max=T_max,presets=presets,**kwargs)
    for key in params:
        results[key+'_init']=res['mmol_per_mg'].drop(['CO2'],axis=1)
    
    results['center_0_minus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
    results['center_0_plus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
    gases=[key for key in presets]

    for key in params:
        for i in  [-1,1]:
            print('{}\n{}: {}'.format('_'*15,key,variance[key]*i))
            temp_presets=copy.deepcopy(presets)
            for gas in gases:
                if key =='center_0':
                    for group in temp_presets[gas].index:
                        print(group)
                        temp_presets[gas].loc[group]=presets[gas].loc[group]+i*np.array([variance[key], 0, 0, variance[key], 0, 0, variance[key], 0, 0])
                        res=fits(*TG_IR,reference=reference,save=False,plot=False,presets=temp_presets,**kwargs)
                        
                        col=group+'_'+gas.upper()
                        if i==-1:
                            results[key+'_minus'][col]=res['mmol_per_mg'][col]
                        elif i==1:
                            results[key+'_plus'][col]=res['mmol_per_mg'][col] 
                
                else: 
                    if key=='hwhm_max':
                        temp_presets[gas]+=i*np.array([0, 0, 0, 0, 0, 0, 0, variance[key], 0])
                    elif key=='tol_c':
                        temp_presets[gas]+=i*np.array([0, 0, 0, -variance[key], 0, 0, variance[key], 0, 0])
                    elif key in ['height_0','hwhm_0']:
                        temp_presets[gas][key]=temp_presets[gas][key[:key.rfind('_')]+'_max']*(default[key]+i*variance[key])
                        
        if key !='center_0':
            res=fits(*TG_IR,reference=reference,save=False,T_max=T_max,presets=temp_presets, **kwargs)
            if i==-1:
                results[key+'_minus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
            elif i==1:
                results[key+'_plus']=res['mmol_per_mg'].drop(['CO2'],axis=1)

    #make subdirectory to save data
    if save:
        path=os.path.join(PATHS['dir_home'],time()+reference+'_{}_{}'.format(var_T,var_rel))
        os.makedirs(path)
        os.chdir(path)
    with pd.ExcelWriter('robustness.xlsx') as writer:
        for key in results:
            results[key].to_excel(writer,sheet_name=key)
    samples=[re.search('^.+(?=_mean)',index).group() for index in res['mmol_per_mg'].index if re.search('^.+(?=_mean)',index)!=None]     
    
    for sample in samples:
        x=['$T_m$','$\Delta T_m$','$HWHM_{max}$','$h_0$','$HWHM_0$']
        x=dict(zip(params,x))
        print(sample)
        
        for param in params:
            data=dict()
            data['minus']=dict({'x':[],'y':[],'yerr':[]})
            data['init']=dict({'x':[],'y':[],'yerr':[]})
            data['plus']=dict({'x':[],'y':[],'yerr':[]})
        
            for key in results:
                for column in res['mmol_per_mg'].drop(['CO2'],axis=1).columns:#'adsorbed_CO2',
                    if key.find(param)!=-1:
                        data[key[key.rfind('_')+1:]].setdefault('x', []).append(column)
                        data[key[key.rfind('_')+1:]].setdefault('y', []).append(results[key].loc[sample+'_mean',column])
                        data[key[key.rfind('_')+1:]].setdefault('yerr', []).append(results[key].loc[sample+'_stddev',column])
            fig=plt.figure()
            plt.title('{}: {}'.format(sample,x[param]))
            for key in data:
                xticks=['{} {}'.format(group[:group.rfind('_')].capitalize(),get_label(group[group.rfind('_')+1:].lower())) for group in data[key]['x']]
                plt.errorbar(xticks,data[key]['y'],yerr=data[key]['yerr'],label='{} {}'.format(key,variance[param] if key!='init' else default[param]),marker='x',capsize=10,ls='none')
            plt.ylim(ylim)
            plt.ylabel('${}\,{}^{{-1}}$'.format(UNITS['molar_amount'],UNITS['sample_mass']))
            plt.legend()
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.show()
            path=os.path.join(PATHS['dir_home'],'Robustness')
            if not os.path.exists(path):
                os.makedirs(path)
            fig.savefig(os.path.join(path,sample+'_'+param+'.png'), bbox_inches='tight', dpi=DPI)
    os.chdir(PATHS['dir_home'])
return
