import pandas as pd
import numpy as np
from .fitting import fits, get_presets
from ..input_output.general import time
from ..plotting import get_label
import os
import matplotlib.pyplot as plt
from ..config import PATHS, BOUNDS, UNITS, DPI
import copy
import time as tm

def robustness(objs,reference,T_max=None,save=True,var_T=10,var_rel=0.3,ylim=[0,None],**kwargs):
    "perform robustness test on multiple TG_IR objects"
    
    # load default presets
    presets_rob=get_presets(PATHS['dir_home'], reference)
    
    # setting up output DataFrames and lists of parameters to test
    params=['center_0','tolerance_center','hwhm_max','height_0','hwhm_0']
    results=dict()
    results['summary']=pd.DataFrame()
    variance=dict(zip(params,[var_T,var_T,var_T,var_rel,var_rel]))
    default=dict(zip(params,[0,BOUNDS.getfloat('tol_center'),BOUNDS.getfloat('hwhm_max'),BOUNDS.getfloat('height_0'),BOUNDS.getfloat('hwhm_0')]))
    
    # Calculating initial results and timing initial fit for estimation of remaining time
    print('Initial results:')
    start=tm.time()
    res=fits(objs,reference=reference,plot=False, save=False, T_max=T_max, presets=presets_rob, **kwargs)
    length=tm.time()-start
    for key in params:
        results[key+'_init']=res['mmol_per_mg']
    del res
    
    gases=[key for key in presets_rob]
    print('\nVarying fitting parameters...\nApproximate remaining time: {:.1f} min'.format((length*2*(4+sum([len(presets_rob[gas]) for gas in presets_rob]))/60)*1.1))
    
    # cycling through parameters to vary
    for key in params:
        for i,suffix in  zip([-1,1],['_minus','_plus']):
            print('{0}\n{0}\n{1} {2:+}:'.format('_'*90,key,variance[key]*i))
            temp_presets=copy.deepcopy(presets_rob)
            
            if key =='center_0':
                results[key+suffix]=pd.DataFrame(columns=results[key+'_init'].columns)
                for gas in gases: 
                    # vary center of each group
                    for group in temp_presets[gas].index:
                        cols=temp_presets[gas].drop('link',axis=1).columns
                        temp_presets=copy.deepcopy(presets_rob)
                        print('\n{} {}:'.format(group.capitalize(),gas))
                        temp_presets[gas].loc[group,cols]=presets_rob[gas].loc[group,cols]+i*np.array([variance[key], 0, 0,variance[key], 0, 0,variance[key], 0, 0])
                        col=group+'_'+gas
                        res=fits(objs,reference=reference,save=False,plot=False,presets=temp_presets,**kwargs)
                        results[key+suffix][col]=res['mmol_per_mg'][col]
                        del res
            else:
                for gas in gases:
                    cols=temp_presets[gas].drop('link',axis=1).columns
                    if key=='hwhm_max':
                        temp_presets[gas].loc[:,cols]+=i*np.array([0, 0, 0, 0, 0, 0, 0,variance[key], 0])
                    elif key=='tolerance_center':
                        temp_presets[gas].loc[:,cols]+=i*np.array([0, 0, 0,-variance[key], 0, 0,variance[key], 0, 0])
                    elif key in ['height_0','hwhm_0']:
                        temp_presets[gas][key]=temp_presets[gas][key[:key.rfind('_')]+'_max']*(default[key]+i*variance[key])
                
                res=fits(objs,reference=reference,plot=False, save=False, T_max=T_max,presets=temp_presets, **kwargs)
                results[key+suffix]=res['mmol_per_mg']
            results[key+suffix].sort_index().sort_index(axis=1,inplace=True)

    # make subdirectory to save data
    if save:
        path=os.path.join(PATHS['dir_home'],'Robustness',time()+reference+'_{}_{}'.format(var_T,var_rel))
        os.makedirs(path)
        os.chdir(path)
        
    # sort and summarize results for each sample
    samples=results['center_0_init'].index.levels[0]
    print('{0}\n{0}\nResults:\n{0}'.format('_'*30)) 
    for sample in samples:
        print(sample)
        # labels and units for plots
        labels=['$center$','$tolerance\,center$','$HWHM_{max}$','$height_0$','$HWHM_0$']
        units=['°C','°C','°C','$height_{max}$','$HWHM_{max}$']
        
        # exclude 'mean' and 'sum' columns from plots
        drop_cols=[gas for gas in gases]
        drop_cols_plot=[col for col in results['center_0_init'].columns if ('_sum' in col) or ('_mean' in col)]
        x=results['center_0_init'].columns.drop(drop_cols+drop_cols_plot)
        
        data=dict()
        data['mean']=pd.DataFrame(columns=x)
        data['all']=pd.DataFrame(columns=x)
        
        # cycle through varied parameters
        for param,label,unit in zip(params,labels,units):
            
            # setup figure
            fig=plt.figure()
            plt.title('{}: {}'.format(sample,label))
            
            # cycle through upper, initial and lower bound of parameter
            for run,i in zip(['plus','init','minus'],[-1,0,1]):
                index='_'.join([param,run])
                
                # make list of mean value and deviation
                y=results[index].loc[sample,'mean',:].drop(drop_cols,axis=1)
                yall=results[index].loc[sample,results[index].index.levels[1].drop(['mean','stddev','dev'],errors='ignore'),:].drop(drop_cols,axis=1)
                yerr=results[index].loc[sample,'dev',:].drop(drop_cols,axis=1)
                
                # plot errorbar
                xticks=['{} {}'.format(group[:group.rfind('_')].capitalize() if group.rfind('_')!=-1 else '',get_label(group[group.rfind('_')+1:].lower())) for group in x]
                plt.errorbar(xticks,y.drop(drop_cols_plot,axis=1).values[0],yerr=yerr.drop(drop_cols_plot,axis=1).values[0],label='{} {}'.format(default[param]+i*variance[param] if param !='center_0' else '{:+}'.format(default[param]+i*variance[param]),unit),marker='x',capsize=10,ls='none')
                
                # collect mean as well as individual results for further statistical evalutaion
                data['mean']=data['mean'].append(y)
                data['all']=data['all'].append(yall)
                
            # draw and save plot
            plt.ylim(ylim)
            plt.ylabel('${}\,{}^{{-1}}$'.format(UNITS['molar_amount'],UNITS['sample_mass']))
            plt.legend()
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.show()
            fig.savefig(sample+'_'+param+'.png', bbox_inches='tight', dpi=DPI)
            
        # make further statistical summary
        results['summary']=results['summary'].append(pd.concat({sample:pd.DataFrame(data['mean'].mean(axis=0).rename('mean')).T}, names=['samples','run']))
        results['summary']=results['summary'].append(pd.concat({sample:pd.DataFrame(data['mean'].std(axis=0).rename('meanstddev')).T}, names=['samples','run']))
        results['summary']=results['summary'].append(pd.concat({sample:pd.DataFrame(data['all'].std(axis=0).rename('stddev')).T}, names=['samples','run']))
        results['summary']=results['summary'].append(pd.concat({sample:pd.DataFrame(data['all'].min(axis=0).rename('min')).T}, names=['samples','run']))
        results['summary']=results['summary'].append(pd.concat({sample:pd.DataFrame(data['all'].max(axis=0).rename('max')).T}, names=['samples','run']))
     
    # save results to excel file
    if save:
        with pd.ExcelWriter('robustness.xlsx') as writer:
            for key in results:  
                results[key].drop(drop_cols,axis=1,errors='ignore').to_excel(writer,sheet_name=key)
    os.chdir(PATHS['dir_home'])
    return 
