import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from ..plotting import get_label
from ..config import UNITS, SEP, DPI, BOUNDS, PATHS, COUPLING
from ..input_output.general import time
import os
import copy
import re


def gaussian(x,height,center,hwhm):
    return height*np.exp(-np.log(2)*np.power((x-center)/hwhm,2))

def multi_gauss(x,*args):
    n=int(len(args)/3)
    heights=args[:n]
    centers=args[n:2*n]
    hwhms=args[2*n:len(args)]
    
    s=0
    for i in range(n):
        s=s+gaussian(x,heights[i],centers[i],hwhms[i])
    return s

def baseline_als(y, lam=1e6, p=0.01, niter=10): #https://stackoverflow.com/questions/29156532/python-baseline-correction-library
    L = len(y)
    D = sp.sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sp.sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = sp.sparse.linalg.spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def fitting(TG_IR,presets,func=multi_gauss,y_axis='orig',plot=False,save=True,predef_tol=0.01):
    data=[]
    for key in presets:
        data+=[presets[key].loc[:,'link'].rename(key)]
    links=pd.concat(data,axis=1)
    gas_links=links.replace('0',np.nan).dropna(thresh=1).dropna(thresh=1,axis=1)
    if gas_links.dropna(axis=1).empty and not links.replace('0',np.nan).dropna(thresh=1).empty:
        print('You cannnot predefine fitting parameters for all supplied gases!')
        gas_links=pd.DataFrame()
        links=pd.DataFrame()
    gases=[key for key in presets]
    for gas in gas_links.columns:
        gases.remove(gas)
        gases.append(gas)
    #thresholds for fit parameters
    ref_mass=TG_IR.info['reference_mass']
    
    #initializing output DataFrame
    peaks=pd.DataFrame(columns=['center','height','hwhm','area','mmol','mmol_per_mg'],index=[group+'_'+gas for gas in gases for group in presets[gas].index]+[gas for gas in gases])
    sumsqerr=pd.DataFrame(index=[TG_IR.info['name']],columns=gases)

    #cycling through gases
    FTIR=TG_IR.ir.copy()
    for gas in gases:
        # correction of water-signal drift
        if gas=='H2O':
            FTIR[gas]-=baseline_als(FTIR[gas])
            
        ## molar desorption
        tot_area=np.sum(TG_IR.ir[gas])
    
        if gas == 'H2O':
            # total area of water is calculated above dry-point
            tot_area=np.sum(TG_IR.ir[gas][TG_IR.ir['sample_temp']>TG_IR.info['dry_temp']])
            
        tot_mol=(tot_area-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]
        peaks['area'][gas]=tot_area
        peaks['mmol'][gas]=tot_mol
        peaks['mmol_per_mg'][gas]=tot_mol/TG_IR.info[ref_mass]
        
        if y_axis=='rel':
            FTIR.update(FTIR[gas]/tot_area*tot_mol)
        
        # predefining guesses and bounds
        if gas in gas_links:
            df=links.dropna(thresh=2).replace('0',np.nan).dropna(thresh=1)
            for group in df.index:
                other=links.loc[group,:].index[links.loc[group]=='0'].values[0]
                for letter in df.loc[group,gas]:
                    if letter=='c':
                        param='center'
                    if letter=='w':
                        param='hwhm'
                    if letter=='h':
                        param='height'
                    if param=='height':
                        preset=peaks.loc['{}_{}'.format(group,other),param]/TG_IR.linreg['slope'][other]*TG_IR.linreg['slope'][gas]
                    else:
                        preset=peaks.loc['{}_{}'.format(group,other),param]
                    presets[gas].loc[group,'{}_0'.format(param)]=preset 
                    presets[gas].loc[group,'{}_min'.format(param)]=preset*(1-predef_tol)
                    presets[gas].loc[group,'{}_max'.format(param)]=preset*(1+predef_tol)

        #guesses 
        params_0=np.concatenate(([presets[gas].loc[:,key+'_0'] for key in ['height', 'center', 'hwhm']])) 
        # ...and bounds
        params_min=np.concatenate(([presets[gas].loc[:,key+'_min'] for key in ['height', 'center', 'hwhm']])) 
        params_max=np.concatenate(([presets[gas].loc[:,key+'_max'] for key in ['height', 'center', 'hwhm']])) 
        

        
        #actual fitting
        x=FTIR['sample_temp']
        try:
            popt,pcov=sp.optimize.curve_fit(func,x,FTIR[gas],p0=params_0,bounds=(params_min,params_max))
        except:
            print('Failed to fit {} signal'.format(gas))
            break
        
        #return values
        num_curves=len(presets[gas])
        for i in range(num_curves):
            group=presets[gas].index[i]+'_'+gas
            peaks['height'][group]=popt[i]
            peaks['center'][group]=popt[i+num_curves]
            peaks['hwhm'][group]=popt[i+2*num_curves]
            if y_axis=='orig':
                peaks['area'][group]=np.sum(gaussian(x,popt[i],popt[i+num_curves],popt[i+2*num_curves]))
                peaks['mmol'][group]=peaks['area'][group]/tot_area*tot_mol
                peaks['mmol_per_mg'][group]=peaks['mmol'][group]/TG_IR.info[ref_mass]
            elif y_axis=='rel':
                peaks['mmol'][group]=peaks['area'][group]/tot_area*tot_mol
                peaks['mmol_per_mg'][group]=peaks['mmol'][group]/TG_IR.info[ref_mass]
        
        ###plotting
        
        profiles=pd.DataFrame()
        data=FTIR[gas]
        fit=multi_gauss(x,*popt)
        diff=data-fit
        sumsqerr[gas][TG_IR.info['name']]=np.sum(np.power(diff,2))
        profiles['sample_temp']=x
        profiles['data']=data
        profiles['fit']=fit
        profiles['diff']=diff
            
        if plot:
            #setup plot
            fig=plt.figure(constrained_layout=True)
            gs = fig.add_gridspec(8, 1)
            fitting = fig.add_subplot(gs[:-1, 0])
            fitting.set_title('{}, {:.2f} mg'.format(TG_IR.info['alias'],TG_IR.info[ref_mass]))
            error = fig.add_subplot(gs[-1,0],sharex=fitting)
            #fitting.xaxis.set_ticks(np.arange(0, 1000, 50))
            
            #plotting of fit
            fitting.plot(x,data,label='data',lw=2,zorder=num_curves+1)#,ls='',marker='x',markevery=2,c='cyan')
            fitting.plot(x,fit,label='fit',lw=2,zorder=num_curves+2)
        for i in range(0,num_curves):
            y=gaussian(x,popt[i],popt[i+num_curves],popt[i+2*num_curves])
            profiles[presets[gas].index[i]]=y
            if plot:
                fitting.text(popt[num_curves+i],popt[i],presets[gas].index[i],zorder=num_curves+3+i)
                fitting.plot(x,y,linestyle='dashed',zorder=i)
        if plot:
            fitting.legend()
            fitting.set_xlabel('{} {} ${}$'.format(UNITS['sample_temp'], SEP, UNITS['sample_temp']))
            if y_axis=='orig':
                fitting.set_ylabel('{} {} ${}$'.format(get_label(gas), SEP, UNITS['ir']))
            elif y_axis=='rel':
                fitting.set_ylabel('{} {} ${}\,{}^{{-1}}\,{}^{{-1}}$'.format(get_label(gas), SEP, UNITS['molar_amount'], UNITS['sample_mass'], UNITS['time']))

            #mark center on x-axis
            fitting.scatter(popt[num_curves:2*num_curves],np.zeros(num_curves),marker=7,color='k',s=100,zorder=num_curves+3)

            #plotting of absolute difference
            abs_max=0.05*max(data)
            
            error.text(0,abs_max,'SQERR: {:.2e}'.format(sumsqerr[gas][TG_IR.info['name']]))#,'SQERR: '+'%.2E'% Decimal(sumsqerr[gas][TG_IR.info['name']]),va='bottom')
            error.plot(x,diff)
            error.hlines(0,min(x),max(x),ls='dashed')
            error.set_xlabel('{} {} ${}$'.format(UNITS['sample_temp'], SEP, UNITS['sample_temp']))
            error.set_ylabel('error {} ${}$'.format(SEP, UNITS['ir']))
            error.set_ylim(-abs_max,abs_max)
            plt.show()
            fig.savefig(TG_IR.info['name']+'_'+gas+'.png', bbox_inches='tight', dpi=DPI)
        if save:
            f_name=TG_IR.info['name']+'_'+y_axis+'.xlsx'
            try:
                with pd.ExcelWriter(f_name,engine='openpyxl', mode='a') as writer:
                    profiles.to_excel(writer,sheet_name=gas)
                    presets[gas].to_excel(writer,sheet_name=gas+'_param')
            except:
                with pd.ExcelWriter(f_name,engine='openpyxl') as writer:
                    profiles.to_excel(writer,sheet_name=gas)
                    presets[gas].to_excel(writer,sheet_name=gas+'_param')
                    
    #calculate summarized groups
    groups=list(set([re.split('_| ',group)[0] for group in peaks.index if group not in gases]))
    for group in groups:
        group_set=peaks[['mmol','mmol_per_mg']].loc[peaks.index.map(lambda x: x.startswith(group))]
        if len(group_set)>1:
            peaks=peaks.append(pd.DataFrame(group_set.sum(axis=0).rename(group+'_sum')).T)
            peaks=peaks.append(pd.DataFrame(group_set.mean(axis=0).rename(group+'_mean')).T)
    
    if save:                
        with pd.ExcelWriter(f_name,engine='openpyxl', mode='a') as writer:
                peaks.astype(float).to_excel(writer,sheet_name='summary')
    return peaks.astype(float),sumsqerr

def fits(TG_IR,reference,save=True,presets=None,**kwargs):
    if presets==None:
        presets=get_presets(PATHS['dir_home'], reference,TG_IR[0].ir)

    #initializing of output DataFrames
    err=pd.DataFrame()
    names=['center','height','hwhm','area','mmol','mmol_per_mg']
    res=dict()
    for name in names:
        res[name]=pd.DataFrame()
    
    #make subdirectory to save data
    if save:
        path=os.path.join(PATHS['dir_fitting'],time()+reference+'_'+'_'.join(list(set([str(obj.info['name']) for obj in TG_IR]))))
        os.makedirs(path)
        os.chdir(path)
    
    sample_re='^.+(?=_\d{1,3})'
    num_re='(?<=_)\d{1,3}$'
    #cycling through samples
    for obj in TG_IR:
        #fitting of the sample and calculating the amount of functional groups
        name=obj.info['alias']
        sample=re.search(sample_re,name)
        num=re.search(num_re,name)
        if sample!=None:
            sample=sample.group()
        else:
            sample=name
        if num!=None:
            num=num.group()
        else:
            num=0

        peaks, sumsqerr=obj.fit(reference,presets=presets,**kwargs,save=False)

        #writing data to output DataFrames
        for key in res:
            res[key]=res[key].append(pd.concat({sample:pd.DataFrame(peaks[key].rename(num)).T}, names=['samples','run']))  
        err=err.append(pd.concat({sample:sumsqerr.rename({name:num})}, names=['samples','run']))
        
    # calculate statistical values
    dm=COUPLING.getfloat('mass_resolution')*1e-3
    for key in res:
        samples=res[key].index.levels[0]
        for sample in samples:
            
            if key=='mmol_per_mg':
                drop_cols=[col for col in res[key].columns if ('_sum' in col) or ('_mean' in col)]
                columns=res[key].columns.drop(drop_cols)
                group_gas=[column[column.rfind('_')+1:] for column in columns]

                lod=[TG_IR[0].stats['x_LOD'][gas] for gas in group_gas]
                subset=res['mmol_per_mg'].loc[sample,columns]
                mmol=res['mmol'].loc[sample,columns]#res['mmol'][column].loc[indices]
                g=mmol/subset
                
                dmmolg_i=np.power(np.power(lod/mmol,2)+np.power(dm/g,2),0.5)*subset
                dmmol=np.power(np.sum(np.power(dmmolg_i,2)),0.5)
            
                stddev=pd.concat([pd.DataFrame(dmmol.rename('dev')).T,pd.DataFrame(res[key].loc[sample,drop_cols].std().rename('dev')).T],axis=1)

            else:
                stddev=pd.DataFrame(res[key].loc[sample].std().rename('stddev')).T
            mean=pd.DataFrame(res[key].loc[sample].mean().rename('mean')).T

            res[key]=res[key].append(pd.concat({sample:mean}, names=['samples','run']))
            res[key]=res[key].append(pd.concat({sample:stddev}, names=['samples','run']))
                    
        res[key].sort_index(inplace=True)       
    
    #exporting data
    if save:
        with pd.ExcelWriter('summary.xlsx') as writer:
            for key in res:
                res[key].dropna(axis=1,thresh=1).to_excel(writer,sheet_name=key)
            err.to_excel(writer,sheet_name='sum_squerr')
        os.chdir(PATHS['dir_home'])
    return res

def get_presets(path,reference,FTIR):
    presets=dict()
    references=pd.read_excel(os.path.join(path,'Fitting_parameter.xlsx'),index_col=0,header=None,sheet_name=None)
    gases=list(set(references['center_0'].loc['gas']))
    
    for gas in gases:
        index=[references['center_0'].loc['group',i] for i in references['center_0'].columns if (references['center_0'].loc['gas',i].upper()==gas)]
        data=pd.DataFrame(index=index)
        for key in references:
            data[key]=pd.DataFrame(references[key].loc[reference,:][references[key].loc['gas',:]==gas].T.values,index=index,columns=[key])#.dropna(axis=1)
        presets[gas]=data.dropna(axis=0,how='all')    
        
        params=['height_0', 
            'hwhm_0', 
            'center_min', 
            'hwhm_min',
            'height_min', 
            'center_max',
            'hwhm_max', 
            'height_max',
            'link']
        vals=[pd.Series(presets[gas].loc[:,'height_max']).fillna(max(FTIR[gas]) if BOUNDS['height_max'] == 'max' else BOUNDS.getfloat('height_max'))*BOUNDS.getfloat('height_0'),
              pd.Series(presets[gas].loc[:,'hwhm_max']).fillna(BOUNDS.getfloat('hwhm_max'))*BOUNDS.getfloat('hwhm_0'),
              pd.Series(presets[gas].loc[:,'center_0']-BOUNDS.getfloat('tol_center')), 
              BOUNDS.getfloat('hwhm_min'),
              BOUNDS.getfloat('height_min'), 
              pd.Series(presets[gas].loc[:,'center_0']+BOUNDS.getfloat('tol_center')), 
              BOUNDS.getfloat('hwhm_max'), 
              max(FTIR[gas]) if BOUNDS['height_max'] == 'max' else BOUNDS.getfloat('height_max'),
              '0']
        infill=dict(zip(params,vals))
        presets[gas]=presets[gas].fillna(infill).dropna()
        if presets[gas].empty:
            del presets[gas]
    return presets
