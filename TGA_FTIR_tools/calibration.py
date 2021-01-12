import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os

from .plotting import get_label
from sklearn import linear_model
from .config import SAVGOL, PARAMS, UNITS, SEP, MOLAR_MASS, COUPLING, PATHS
from .input_output import TGA, FTIR, corrections


def mass_step(TGA_data,rel_height=.98,plot=False): #rel_height=.963
    "deriving mass steps via peaks in DTG signal"
    #calculation and smoothing of DTG
    TG=(TGA_data['sample_mass']/TGA_data['sample_mass'][0])
    DTG=sp.signal.savgol_filter(TG, SAVGOL.getint('window_length'), SAVGOL.getint('polyorder'), deriv=1)
  
    #detect mass steps
    peaks,properties=sp.signal.find_peaks(-DTG,height=0.00025,width=0,rel_height=rel_height,prominence=0.00025)
    step_end=properties['right_ips'].astype(np.int, copy=False)
    step_start=properties['left_ips'].astype(np.int, copy=False)
    
    #calculate masssteps
    steps=np.zeros(len(peaks))
    samples=20
    for i in range(len(peaks)):
        steps[i]=np.mean(TGA_data['sample_mass'][step_end[i]:step_end[i]+samples])
    
    #calculate step height
    step_height=np.zeros(len(steps))
    steps=np.insert(steps,0,TGA_data['sample_mass'][0])
        
    step_height=np.zeros(len(step_start))
    for i in range(len(step_start)):
        step_height[i]=np.mean(TGA_data['sample_mass'][step_start[i]-samples:step_start[i]])-np.mean(TGA_data['sample_mass'][step_end[i]:step_end[i]+samples])
    
    rel_step_height=step_height/steps[0]*100
    
    #plotting
    if plot==True:
        #plotting of rel. TG
        x=TGA_data['sample_temp']
        plt.figure()
        rel_steps=steps/steps[0]*100
        plt.hlines(rel_steps[:-1],np.zeros(len(rel_steps)-1),x[step_end],linestyle='dashed')
        plt.vlines(x[step_end],rel_steps[1:],rel_steps[:-1],linestyle='dashed')
        for i in range(len(step_end)):
            plt.text(x[step_end[i]]+5,rel_steps[i+1]+rel_step_height[i]/2,str(round(rel_step_height[i],2))+' %')
        plt.plot(x,TGA_data['sample_mass']/TGA_data['sample_mass'][0]*100)
        plt.text(0.85*max(TGA_data['sample_temp']),100,'sample mass: {:.2f} {}'.format(TGA_data['sample_mass'][0],UNITS['sample_mass']), horizontalalignment='center')
        plt.xlabel('{} {} {}'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        plt.ylabel('{} {} %'.format(PARAMS['sample_mass'],SEP))
        plt.title('TG')
        plt.show()
        
        #plotting of DTG
        plt.figure()
        y=-DTG
        plt.plot(x,y)
        plt.vlines(x[step_end],0,max(y),linestyle='dashed')
        plt.vlines(x[step_start],0,max(y),linestyle='dashed')
        plt.vlines(x[peaks],y[peaks]-properties['peak_heights'],y[peaks])
        plt.hlines(y[peaks]-properties['peak_heights'],x[step_end],x[step_start])
        plt.xlabel('{} {} {}'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        plt.ylabel('{} {} {} ${}^{{-1}}$'.format(PARAMS['dtg'],SEP,UNITS['sample_mass'],UNITS['time']))
        plt.title('DTG')
        plt.show()
        
    return step_height,rel_step_height,step_start,step_end

def integrate_peaks(FTIR_data,step_start,step_end,corr_baseline=None,plot=False,gases=None):
    "integrating IR signal in between given bounds"
    
    # adsjusting step starts according to coupling delay
    step_start=step_start+int(60*COUPLING.getfloat('coupling_delay'))
    step_end=step_end+int(60*COUPLING.getfloat('coupling_delay'))
    integrals=pd.DataFrame(index=range(len(step_start)),columns=gases)
    
    #plotting
    if plot:
        colors =plt.rcParams['axes.prop_cycle'].by_key()['color']
    
        x=FTIR_data['time']/60
        #setup figure and plot first gas
        graph=[None] 
        fig, graph[0]=plt.subplots()
        fig.subplots_adjust(right=.8)
        
        graph[0].set_xlabel('{} {} {}'.format(PARAMS['time'], SEP, UNITS['time']))
        graph[0].set_ylabel('{} {} {}'.format(get_label(gases[0]),SEP,UNITS['ir']))
        graph[0].yaxis.label.set_color(colors[0])
        graph[0].plot(x,FTIR_data[gases[0]])
        
        #append secondary, third... y-axis on right side
        for i,gas in enumerate(gases[1:]):
            y=FTIR_data[gas]
            
            graph.append(graph[0].twinx())
            graph[i+1].spines['right'].set_position(('axes',1+i*.1))
            graph[i+1].plot(x,y, color=colors[i+1])
            graph[i+1].vlines(step_start/60,0,max(y),linestyle='dashed')
            graph[i+1].vlines(step_end/60,0,max(y),linestyle='dashed')
            graph[i+1].set_ylabel('{} {} {}'.format(get_label(gas),SEP,UNITS['ir']))
            graph[i+1].yaxis.label.set_color(colors[i+1])
    
    #integration
    for gas in gases:
        for i in range(len(step_end)):
            subset=FTIR_data[gas][(FTIR_data['time']>=step_start[i]) & (FTIR_data['time']<=step_end[i])]
            #baseline correction
            if corr_baseline=='linear':
                baseline=np.linspace(subset.iloc[0],subset.iloc[len(subset)-1],len(subset))
            elif corr_baseline=='const':
                baseline=min(subset.iloc[0],subset.iloc[len(subset)-1])*np.ones(len(subset))#np.linspace(min(subset),min(subset),len(subset))
            elif corr_baseline==None:
                baseline=np.zeros(len(subset))
                
            integral=sp.integrate.simps(subset-baseline)   
            integrals.loc[i,gas]=integral
            
            if plot==True:
                x=FTIR_data['time'][(FTIR_data['time']>=step_start[i]) & (FTIR_data['time']<=step_end[i])]/60
                graph[gases.index(gas)].plot(x,baseline,color=colors[gases.index(gas)],linestyle='dashed')

    if plot==True:
        plt.show()
        
    return integrals 

def eval_lin(x,slope,intercept):
    "evaluating linear equation with slope and intercept at ponts x"
    return slope*x+intercept

def calibration_stats(x_cali,y_cali,linreg,alpha=.95,beta=None,m=1,k=3):
    "calculating calibration stats according to DIN 32645"
    gases=linreg.index
    
    n=len(x_cali)
    if n==2:
        return pd.DataFrame()
    
    f=n-2
    if beta==None:
        beta=alpha
    
    stats=pd.DataFrame()
    for gas in gases:
        b=linreg['slope'][gas]
        a=linreg['intercept'][gas]
        
        s_yx=np.sqrt(np.sum(np.power(b*x_cali[gas]+a-y_cali[gas],2))/(n-2))
        s_x0=s_yx/b
        x_=np.mean(x_cali[gas])
        Q_x=np.sum(np.power(x_cali[gas]-x_,2))
        
        x_NG=s_x0*sp.stats.t.ppf(alpha,f)*np.sqrt(1/m+1/n+(x_*x_)/Q_x)
        
        #x_EG=x_NG+s_x0*sp.stats.t.ppf(beta,f)*np.sqrt(1/m+1/n+(x_*x_)/Q_x)
    
        x_BG=k*x_NG
                
        stats=stats.append(pd.DataFrame([[s_yx,s_x0,x_NG,x_BG]],index=[gas],columns=['s_yx','s_x0','x_LOD','x_LOQ']))
    return stats
    
def calibrate(plot=False,mode='load',method='max'):
    # check if calibration folder is present
    if os.path.exists(PATHS['dir_calibration'])==False:
        # make directory
        os.makedirs(PATHS['dir_calibration'])
    os.chdir(PATHS['dir_calibration'])
    if mode=='load':
        # try to load saved calibration
        try:
            cali=pd.read_excel('cali.xlsx',sheet_name=None,index_col=0)
            linreg=cali['linreg']
            x_cali=cali['x in {}'.format(UNITS['molar_amount'])]
            y_cali=cali['y in {}'.format(UNITS['int_ir'])]
            stats=cali['stats']
            data=pd.read_excel('cali.xlsx',sheet_name='data',index_col=[0,1])
            gases=linreg.index
        except:
            print('No calibration data found. To obtain quantitative IR data supply an \'Calibration\' folder in the home directory containing cali.xlsx or run TGA_FTIR_tools.calibrate(mode=\'recalibrate\')!')
            os.chdir(PATHS['dir_home'])
            return
    
    #new calibration
    elif mode=='recalibrate':
        #setting up output DataFrames
        x_cali=pd.DataFrame()
        y_cali=pd.DataFrame()
        linreg=pd.DataFrame()
        stats=pd.DataFrame()
        data=pd.DataFrame()
        print('Calibrating...')
        #imporitng sample list for calibration
        try:
            samples=pd.read_csv('Sample_list.txt',delimiter='\t')
        except:
            with open('Sample_list.txt', 'w') as file:
                file.write('Samples\tBaseline')
            print('\'Sample_list.txt\' was created in the \'Calibration\' folder, please fill in calibration measurements and rerun this command.')
            os.chdir(PATHS['dir_home'])
            return
        #calculating mass steps and integrating FTIR_data signals for all samples
        for sample,baseline in zip(samples['Samples'],samples['Baseline']):
            
            #reading and correcting data
            TGA_data=TGA.read_TGA(sample)
            info=TGA.TGA_info(sample,TGA_data)
            TGA_data=corrections.corr_TGA(TGA_data,baseline)
            FTIR_data=FTIR.read_FTIR(sample)
            info['gases']=FTIR_data.columns[1:].to_list()
            gases=info['gases']
            FTIR_data.update(corrections.corr_FTIR(FTIR_data,baseline,plot=plot))
            try:
                FTIR_data['time']+=60*info['background_delay']
            except:
                FTIR_data['time']+=60*COUPLING.getfloat('background_delay')
            print('----------------------------------------------------\n{}'.format(info['name']))
        
            #calculating mass steps and integrating FTIR_data signals 
            [steps,rel_steps,stepstart,stepend]=mass_step(TGA_data,plot=plot)
            integrals=integrate_peaks(FTIR_data,stepstart,stepend,plot=plot,corr_baseline=None,gases=gases)
            
            integrals.insert(loc=0,column='mass loss in {}'.format(UNITS['sample_mass']),value=steps)
            data=data.append(pd.concat({sample: integrals}, names=['samples','step']))
        
        #assigning gases to mass steps
        for sample in data.index.levels[0]:
            release_steps=[]
            for i,step in enumerate(data.loc[sample,'mass loss in {}'.format(UNITS['sample_mass'])]):
                integrals=data.loc[sample].drop(['mass loss in {}'.format(UNITS['sample_mass'])],axis=1)
                norm=integrals.divide(integrals.max(axis=0).values,axis=1).loc[i]
                gas=norm.loc[norm==1].index.values[0]
                release_steps.append(gas)
                if gas not in x_cali.columns:
                    x_cali[gas]=np.nan
                    y_cali[gas]=np.nan
                if sample not in x_cali.index:
                    x_cali=x_cali.append(pd.DataFrame(index=[sample]))
                    y_cali=y_cali.append(pd.DataFrame(index=[sample]))
                x_cali[gas][sample]=step
                y_cali[gas][sample]=integrals.loc[i,gas]
        
        cols=['slope','intercept','r_value','p_value','std_error']
        if method=='iter':
            for gas in gases:
                #slope, intercept, r_value, p_value, std_err=sp.stats.linregress(x_cali[gas].dropna(axis=0).astype(float),y_cali[gas].dropna(axis=0).astype(float))
                x=x_cali[gas].dropna(axis=0).astype(float)
                y=y_cali[gas].dropna(axis=0).astype(float)
                regression=pd.DataFrame([sp.stats.linregress(x,y)],index=[gas],columns=cols)
    
                if gas not in linreg.index:
                    linreg=linreg.append(regression,verify_integrity=True)
                else:
                    linreg.loc[[gas]]=regression
                
            n_iter=10
            X_cali=pd.DataFrame()
            temp_linreg=pd.DataFrame(index=release_steps,columns=cols)
            for i in range(n_iter):
                for step in data.index.levels[1]:
                    gas=release_steps[step]
                    X_cali[gas]=data.loc[(slice(None),step),'mass loss in {}'.format(UNITS['sample_mass'])].droplevel(1)
                    for other in set(release_steps) - set([gas]):
                        corr=(data.loc[(slice(None),step),other].droplevel(1)-linreg['intercept'][other])/linreg['slope'][other]
                        X_cali[gas]=np.subtract(X_cali[gas],corr*(corr>0))
                        
                    x=X_cali[gas].values.astype(float)
                    y=data.loc[(slice(None),step),gas].values.astype(float)
                    temp_linreg.update(pd.DataFrame(np.array([sp.stats.linregress(x,y)]),columns=cols,index=[gas]))
            
                linreg.update(temp_linreg)
                
            for step in data.index.levels[1]:
                    gas=release_steps[step]
                    for other in set(release_steps) - set([gas]):
                        x_cali[gas]=x_cali[gas].subtract((data.loc[(slice(None),step),other].droplevel(1)-linreg['intercept'][other])/linreg['slope'][other])
            
        #regression
        for gas in gases:
            x_cali.update(x_cali[gas]/(MOLAR_MASS.getfloat(gas)))
            x=x_cali[gas].dropna(axis=0).astype(float)
            y=y_cali[gas].dropna(axis=0).astype(float)
            regression=pd.DataFrame([sp.stats.linregress(x,y)],index=[gas],columns=cols)

            if gas not in linreg.index:
                linreg=linreg.append(regression,verify_integrity=True)
            else:
                linreg.loc[[gas]]=regression

        if method=='co_oxi':
            x_cali['CO'].update(x_cali['CO']-((data.loc[(slice(None),1),'CO2']-linreg.loc['CO2','intercept'])/linreg.loc['CO2','slope']).values)
            
            x=x_cali['CO'].dropna(axis=0).astype(float)
            y=y_cali['CO'].dropna(axis=0).astype(float)
            regression=pd.DataFrame([sp.stats.linregress(x,y)],index=['CO'],columns=cols)
            linreg.loc[['CO']]=regression
        
        if method=='mlr':
            Y_cali=data.loc[(slice(None),slice(None)),'mass loss in {}'.format(UNITS['sample_mass'])]
            X_cali=data.drop(['mass loss in {}'.format(UNITS['sample_mass'])],axis=1)
            mlr=linear_model.LinearRegression()#RANSACRegressor()#fit_intercept=0.0)
            mlr.fit(X_cali,Y_cali)
            
            for i, gas in enumerate(X_cali.columns):
                linreg.loc[[gas]]=pd.DataFrame([[1/mlr.coef_[i]*MOLAR_MASS.getfloat(gas),0,np.nan,np.nan,np.nan]],index=[gas],columns=cols)
        
        stats=calibration_stats(x_cali,y_cali,linreg)
        
        #saving of 
        try:
            with pd.ExcelWriter('cali.xlsx') as writer:
                linreg.to_excel(writer,sheet_name='linreg')
                x_cali.to_excel(writer,sheet_name='x in {}'.format(UNITS['molar_amount']))
                y_cali.to_excel(writer,sheet_name='y in {}'.format(UNITS['int_ir']))
                stats.to_excel(writer,sheet_name='stats')
                data.to_excel(writer,sheet_name='data')
        except:
            print('Could not write on cali.xlsx. Close file and rerun command.')
            
    #plotting
    if plot:
        for gas in gases:
            fig=plt.figure()
            x=x_cali[gas]
            y=y_cali[gas]

            plt.scatter(x,y,label='data (N = {})'.format(len(x)))

            x_bounds=np.array((min(x),max(x)))
            plt.plot(x_bounds,x_bounds*linreg['slope'][gas]+linreg['intercept'][gas],label='regression',ls='dashed')
            plt.text(max(x),min(y),'y = {:.3f} $\cdot$ x {:+.3f}, $R^2$ = {:.5}'.format(linreg['slope'][gas],linreg['intercept'][gas],linreg['r_value'][gas]**2),horizontalalignment='right')
            plt.xlabel(UNITS['molar_amount'])
            plt.ylabel(UNITS['int_ir'])
            plt.title(get_label(gas))
            plt.xlim(0,max(x)+abs(min(x)))
            plt.legend(loc=0)
            plt.show()

            Y_cali=x.mul(linreg['slope'][gas]).add(linreg['intercept'][gas])
            plt.scatter(Y_cali,y-Y_cali,label='data (N = {})'.format(len(x)))
            plt.hlines(0,min(Y_cali),max(Y_cali))
            plt.xlabel('$\hat{{y}}_i$ {} {}'.format(SEP,UNITS['int_ir']))
            plt.ylabel('$y_i-\hat{{y}}_i$ {} {}'.format(SEP,UNITS['int_ir']))
            plt.title('Residual plot: {}'.format(get_label(gas)))
            plt.legend(loc=0)
            plt.show()

        plt.figure()
        for gas in gases:
            x=x_cali[gas]
            y=y_cali[gas]

            plt.scatter(x,y,label='data {} (N = {})'.format(get_label(gas),len(x)))

            x=np.array((min(x),max(x)))
            plt.plot(x,x*linreg['slope'][gas]+linreg['intercept'][gas],label='regression {}'.format(get_label(gas)),ls='dashed')
            plt.xlabel(UNITS['molar_amount'])
            plt.ylabel(UNITS['int_ir'])
            plt.xlim(0,max(x)+abs(min(x)))
            plt.legend(loc=0)
        plt.show()
    if mode=='recalibrate':
        print('Calibration completed, data is stored under {}/cali.xlsx'.format(PATHS['dir_calibration']))

    os.chdir(PATHS['dir_home'])
    return linreg,stats
