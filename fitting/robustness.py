import pandas as pd
import numpy as np
from TG_IR import fits
from Functions.general import time
from Functions.plotting import get_label
import os
import re
import configparser
import matplotlib as plt
config = configparser.ConfigParser()
config.read('settings.ini')


def robustness(*TG_IR,reference,tol_c=30,tol_hwhm_0=95,T_max=None,save=True):
    var_T=10
    var_rel=.3
    params=['offs_c','tol_c','tol_hwhm','offs_h','offs_hwhm']
    x=['$T_m$','$\Delta T_m$','$HWHM_{max}$','$h_0$','$HWHM_0$']
    x=dict(zip(params,x))
    variance=[var_T,var_T,var_T,var_rel,var_rel]
    values=dict(zip(params,variance))
    
    units=dict(zip(params,['°C','°C','°C','% $h_{max}$','% $HWHM_{max}$']))
    #glabelpre=dict(zip(['minus','init','plus'],['-','initial','+']))
    
    initials=dict(zip(params,['initial',tol_c,tol_hwhm_0,50,50]))
    minus=dict(zip(params,[-values['offs_c'],tol_c-values['tol_c'],tol_hwhm_0-values['tol_hwhm'],int(np.ceil(initials['offs_h']-100*values['offs_h'])),int(np.ceil(initials['offs_hwhm']-100*values['offs_hwhm']))]))
    plus=dict(zip(params,['+'+str(values['offs_c']),tol_c+values['tol_c'],tol_hwhm_0+values['tol_hwhm'],int(np.ceil(initials['offs_h']+100*values['offs_h'])),int(np.ceil(initials['offs_hwhm']+100*values['offs_hwhm']))]))
    
    glabels=dict(zip(['minus','init','plus'],[minus,initials,plus]))
    results=dict()
    
    #Init values
    res=fits(*TG_IR,reference=reference,tol_c=tol_c,tol_hwhm=tol_hwhm_0,T_max=T_max,init_offs=[0,0,0],save=False,plot=False)
    for key in params:
        results[key+'_init']=res['mmol_per_mg'].drop(['CO2'],axis=1)
    
    results['offs_c_minus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
    results['offs_c_plus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
    
    for key in ['offs_h']:#params:
        for i in  [-1,1]:
            tol_center=tol_c
            tol_hwhm=tol_hwhm_0
            offs_c=0
            offs_h=0
            offs_hwhm=0
            if key=='offs_c':
                
                print('tol_center {}, tol_hwhm {}, offs_c {}, offs_h {}, offs_hwhm {}'.format(tol_center,tol_hwhm,values[key]*i,offs_h,offs_hwhm))
                references=pd.read_excel(os.path.join(config['paths']['dir_configuration'],'surfgr_C.xlsx'),index_col=0,header=None)
                if T_max==None:
                    T_max=min([max(obj.tga['Ts']) for obj in TG_IR])

                elif T_max>min([max(obj.tga['Ts']) for obj in TG_IR]):
                    T_max=T_max=min([max(obj.tga['Ts']) for obj in TG_IR])
                    print('$T_{max}$ exceeds maximum temperature of data')
                
                for col in references.loc[['gas',reference]].dropna(axis=1):
                    print(col)
                    init_params=references.loc[['gas',reference]].dropna(axis=1)
                    #find needed gases
                    gases=list(set(init_params.loc['gas']))
                    init_params[col][reference]+=values[key]*i
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
                    
                    res=fits(*TG_IR,tol_c=tol_center,tol_hwhm=tol_hwhm,T_max=T_max,init_offs=[offs_c,offs_h,offs_hwhm],save=False,plot=False,temps=temps,labels=labels)
                    gas=init_params[col]['gas']
                    group=references[col]['group']+'_'+gas.upper()
                    
                    if i==-1:
                        results[key+'_plus'][group]=res['mmol_per_mg'][group]
                    elif i==1:
                        results[key+'_minus'][group]=res['mmol_per_mg'][group] 
                    break
                    
            else: 
                if key=='tol_hwhm':
                    tol_hwhm+=values[key]*i
                elif key=='tol_c':
                    tol_center+=values[key]*i
                elif key=='offs_h':
                    offs_h=values[key]*i
                elif key=='offs_hwhm':
                    offs_hwhm=values[key]*i
                print('tol_center {}, tol_hwhm {}, offs_c {}, offs_h {}, offs_hwhm {}'.format(tol_center,tol_hwhm,offs_c,offs_h,offs_hwhm))
                res=fits(*TG_IR,reference=reference,tol_c=tol_center,tol_hwhm=tol_hwhm,T_max=T_max,init_offs=[offs_c,offs_h,offs_hwhm],save=True,plot=False)
                if i==-1:
                    results[key+'_plus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
                elif i==1:
                    results[key+'_minus']=res['mmol_per_mg'].drop(['CO2'],axis=1)
            #break
        #break
    #make subdirectory to save data
    if save==True:
        path=os.path.join(config['paths']['dir_configuration'],time()+reference+'_{}_{}'.format(tol_c,tol_hwhm_0))
        os.makedirs(path)
        os.chdir(path)
    with pd.ExcelWriter('robustness.xlsx') as writer:
        for key in results:
            results[key].to_excel(writer,sheet_name=key)
        #err.to_excel(writer,sheet_name='sum_squerr')
    #print(results)
    #print(res['mmol_per_mg'].index)
    samples=[re.search('^.+(?=_mean)',index).group() for index in res['mmol_per_mg'].index if re.search('^.+(?=_mean)',index)!=None]
    
##    for sample in samples:
#        print(sample)
#        for column in res['mmol_per_mg'].columns:
#            fig=plt.figure()
#            plt.title('{}: {} {}'.format(sample,column[:column.rfind('_')].capitalize(),get_label(column[column.rfind('_')+1:].lower())))
#            data=dict()
#            data['minus']=dict({'x':[],'y':[],'yerr':[]})
#            data['init']=dict({'x':[],'y':[],'yerr':[]})
#            data['plus']=dict({'x':[],'y':[],'yerr':[]})
#            for key in results:
#                #print(key[key.rfind('_')+1:])
#                data[key[key.rfind('_')+1:]].setdefault('x', []).append(re.search('^.+(?=_p|_m|_init)',key).group())
#                data[key[key.rfind('_')+1:]].setdefault('y', []).append(results[key].loc[sample+'_mean',column])
#                data[key[key.rfind('_')+1:]].setdefault('yerr', []).append(results[key].loc[sample+'_stddev',column])
#            for key in data:
#                plt.errorbar([x[param] for param in data[key]['x']],data[key]['y'],yerr=data[key]['yerr'],label=labels[key],marker='x',capsize=10,ls='none')
#            plt.ylim(0,None)
 ##           plt.ylabel('$mmol/g$')
 #           plt.legend()
 #           plt.show()
#            fig.savefig(os.path.join(dir_plots,'Robustness',sample+'_'+column+'.png'), bbox_inches='tight', dpi=300)
                #plt.errorbar(re.search('^.+(?=_p|_m|_init)',key).group(0),1,yerr=.5,label=label)
            #plt.errorbar(['test'],[1],yerr=[.5],label=label)       
    
    for sample in samples:
        print(sample)
        
        for param in params:
            data=dict()
            data['minus']=dict({'x':[],'y':[],'yerr':[]})
            data['init']=dict({'x':[],'y':[],'yerr':[]})
            data['plus']=dict({'x':[],'y':[],'yerr':[]})
        
            for key in results:
                #print(key)
                for column in res['mmol_per_mg'].drop(['CO2'],axis=1).columns:#'adsorbed_CO2',
                    if key.find(param)!=-1:
                        data[key[key.rfind('_')+1:]].setdefault('x', []).append(column)
                        data[key[key.rfind('_')+1:]].setdefault('y', []).append(results[key].loc[sample+'_mean',column])
                        data[key[key.rfind('_')+1:]].setdefault('yerr', []).append(results[key].loc[sample+'_stddev',column])
                   # print(data)
            fig=plt.figure()
            plt.title('{}: {}'.format(sample,x[param]))
            #print(data)
            for key in data:
                #print(data[key])
                plt.errorbar(['{} {}'.format(group[:group.rfind('_')].capitalize(),get_label(group[group.rfind('_')+1:].lower())) for group in data[key]['x']],data[key]['y'],yerr=data[key]['yerr'],label='{} {}'.format(glabels[key][param],units[param]),marker='x',capsize=10,ls='none')
            plt.ylim(0,None)
            plt.ylabel('$mmol/g$')
            plt.legend()
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.show()
            path=os.path.join(config['paths']['dir_home'],'Robustness')
            if not os.exists(path):
                os.makedirs(path)
            fig.savefig(os.path.join(path,sample+'_'+param+'.png'), bbox_inches='tight', dpi=300)
            #plt.errorbar(re.search('^.+(?=_p|_m|_init)',key).group(0),1,yerr=.5,label=label)
        #plt.errorbar(['test'],[1],yerr=[.5],label=label)
    os.chdir(config['paths']['dir_home'])
