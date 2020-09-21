import os
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from .general import find_files
from ..config import PATHS, DPI
from .FTIR import read_FTIR

def corr_TGA(TGA,file_baseline,plot=False):
    corr_data=TGA.copy()
    path_baseline=find_files(file_baseline,'.txt',PATHS['dir_data'])[0]
    #opens the buoyancy blank value 'baseline' and substracts them from the original data 
    try:
        reference_mass=pd.read_csv(path_baseline, delim_whitespace=True,decimal=',' ,names=['Index','t','Ts','Tr','mass'],skiprows=13, skipfooter=11,converters={'mass': lambda x: float(x.replace(',','.'))},engine='python').drop(columns='Index')
        corr_data['mass']=corr_data['mass'].subtract(reference_mass['mass'][:len(TGA)])
    except:
        print('>',path_baseline,' was not found.')
        return None
    try:
        path_mW=find_files(file_baseline,'_mW.txt',PATHS['dir_data'])[0]
        reference_heat_flow=pd.read_csv(path_mW, delim_whitespace=True,decimal=',' ,names=['Index','t','Ts','Tr','heat_flow'],skiprows=13, skipfooter=11, usecols=['heat_flow'],engine='python')
        corr_data['heat_flow']=corr_data['heat_flow'].subtract(reference_heat_flow['heat_flow'][:len(TGA)])
    except:
        pass
    
    #plotting of data, baseline and corrected value
    if plot==True:
        plt.figure()
        x=TGA['Ts']
        y=TGA['mass']
        plt.plot(x,y,label='data')
        plt.plot(x,reference_mass['mass'][:len(TGA)],label='baseline')
        plt.plot(x,corr_data['mass'],label='corrected')
        plt.xlabel('T /°C')
        plt.ylabel(y.name+' /mg')
        plt.legend()
        plt.title('TGA baseline correction')
        plt.show()
        
        
    return corr_data

def corr_FTIR(FTIR,file_baseline,save=False,plot=False):
    #opens FTIR data of the baseline and takes the 'co2' column
    corr_data=FTIR.copy()
    try:
        baseline=read_FTIR(file_baseline)
        baseline=np.array(baseline['co2'])
    
        #in the baseline the peaks and valleys as well as the amplitude of the baseline are determined
        peaks_baseline,properties_baseline=sp.signal.find_peaks(baseline,height=[None,None])#,height=[tol*min(baseline),tol*max(baseline)])
        valleys_baseline,valley_properties_baseline=sp.signal.find_peaks(-baseline,height=[None,None])#,height=[tol*min(-baseline),tol*max(-baseline)])
        amplitude_baseline=np.mean(properties_baseline['peak_heights'])+np.mean(valley_properties_baseline['peak_heights'])
    
        #in the original data the peaks and valleys that have similar height as the baseline are determined
        tol=1.5
        peaks,properties=sp.signal.find_peaks(FTIR['co2'],height=[-tol*amplitude_baseline,tol*amplitude_baseline])
        valleys,valley_properties=sp.signal.find_peaks(-FTIR['co2'],height=[None,None],prominence=amplitude_baseline*.05) 
    
        #the median distance between between baseline-peaks, the period is determined
        dist_peaks=np.diff(peaks_baseline)
        len_period=int(np.median(dist_peaks))
    
        #determination of the phase shift in x direction by checking if there is also a valley in the baseline in proximity
        #the x shift is calculated as the median of the differences
        dists=[]
        j=0
        for valley in valleys:
            while j<len(valleys_baseline)-1 and (valleys_baseline[j]-valley <=0):
                j=j+1
            if valleys_baseline[j]-valley>=0:
                dists.append(valleys_baseline[j]-valley)
    
        x_shift=int(sp.stats.mode(dists)[0])
    
        ###modifying of the baseline  
        #elongating the baseline by one period
        period=baseline[:len_period]
        co2_baseline=np.concatenate((period,baseline), axis=None)
    
        #shifting the baseline in x direction
        c=[]
        for x_offs in range(-1,len_period%x_shift+1):
            peaks,props=sp.signal.find_peaks(FTIR['co2']-co2_baseline[x_shift+x_offs:len(FTIR)+x_shift+x_offs],height=[None,None],prominence=amplitude_baseline*.02)
            c.append(len(peaks))
        x_offs=np.where(c==np.min(c))[0][0]-1
    
        co2_baseline=co2_baseline[x_shift+x_offs:len(FTIR)+x_shift+x_offs]
        
    except:
        print('No CO2 baseline found.')
        co2_baseline=np.zeros(len(FTIR))

    h2o_baseline=np.zeros(len(FTIR))
    
    ##return value
    
    corr_data['co']=corr_data['co'].subtract(min(corr_data['co']))
    corr_data['co']=corr_data['co'].subtract(const_baseline(corr_data['co'],0.002))
    corr_data['co2']=corr_data['co2'].subtract(co2_baseline)
    corr_data['co2']=corr_data['co2'].subtract(min(corr_data['co2']))
    corr_data['co2']=corr_data['co2'].subtract(const_baseline(corr_data['co2'],0.07))
    corr_data['h2o']=corr_data['h2o'].subtract(h2o_baseline)
    corr_data['h2o']=corr_data['h2o'].subtract(min(corr_data['h2o']))
    corr_data['h2o']=corr_data['h2o'].subtract(const_baseline(corr_data['h2o'],0.001))
    
    
    ###plotting of baseline, data and the corrected data
    if plot==True:
        
        fig=plt.figure()
        try:
            x=FTIR['Ts']
        except:
            x=FTIR['t']
        y=co2_baseline

        
        plt.plot(x,FTIR['co2'],label='data')
        plt.plot(x,baseline[:len(x)], label='baseline')
        plt.plot(x,y,label='corr. baseline')#,alpha=.5)
        plt.plot(x,corr_data['co2'],label='corr. data')
        plt.hlines(0,min(x),max(x),ls='dashed')
        
        plt.vlines(x.iloc[valleys],min(co2_baseline),max(FTIR['co2']-y),linestyle='dashed')
        plt.legend()
        if x.name=='t':    
            plt.xlabel(x.name+' /min')
        elif x.name=='Ts':    
            plt.xlabel('T /°C')
        plt.ylabel('$CO_2$ /AU')
        plt.title('$CO_2$ baseline correction')
        fig.savefig(os.path.join(PATHS['dir_plots'],'IR','IR_corr.png'), bbox_inches='tight',dpi=DPI)
        plt.show()
        
        fig=plt.figure()
        
        y=h2o_baseline
        plt.plot(x,FTIR['h2o'],label='data')
        plt.plot(x,y,label='baseline')
        plt.plot(x,corr_data['h2o'],label='corr. data')
        plt.hlines(0,min(x),max(x),ls='dashed')
        
        plt.legend()
        if x.name=='t':    
            plt.xlabel(x.name+' /min')
        elif x.name=='Ts':    
            plt.xlabel('T /°C')
        plt.ylabel('AU')
        plt.title('$H_2O$ baseline correction')
        plt.show()

    
    if save==True:
        out=pd.DataFrame()
        out['t']=x
        out['data']=FTIR['co2']
        out['baseline']=baseline[:len(x)]
        out['corr_baseline']=co2_baseline
        out['corr_data']=corr_data['co2']
        out.to_excel('baseline.xlsx')
    return corr_data



def const_baseline(data,thres):
    baseline=data[data<thres]
    
    if len(baseline)==0:
        return 0
    else:
        return np.sum(baseline)/len(baseline)

    
