import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
from ..plotting import get_label
from ..config import UNITS, SEP, DPI, BOUNDS

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

def fitting(TG_IR,centers,labels,func,tol_center=BOUNDS.getfloat('delta_tm'),max_hwhm=BOUNDS.getfloat('hwhm_max'),y_axis='orig',plot=False,save=True,init_offs=[0,0,0]):
    gases=list(centers.columns.values)
    if ('co2' in gases) and ('co' in gases):
        gases.remove('co2')
        gases.insert(0,'co2')
    #thresholds for fit parameters
    ref_mass=TG_IR.info['reference_mass']
    
    #initializing output DataFrame
    peaks=pd.DataFrame(columns=['center','height','hwhm','area','mmol','mmol_per_mg'],index=[group+'_'+gas.upper() for gas in gases for group in labels[gas].dropna()]+[gas.upper() for gas in gases])
    sumsqerr=pd.DataFrame(index=[TG_IR.info['name']],columns=gases)

    #cycling through gases
    FTIR=TG_IR.ir.copy()
    for gas in gases:
        #print(gas)
        params=pd.DataFrame(columns=['min_center','init_center','max_center','min_height','init_height','max_height','min_hwhm','init_hwhm','max_hwhm'],index=[group+'_'+gas.upper() for group in labels[gas].dropna()])

        if gas=='h2o':
            FTIR[gas]-=baseline_als(FTIR[gas])
            
        #molar desorption
        tot_area=np.sum(TG_IR.ir[gas])
        if gas == 'h2o':
            tot_area=np.sum(TG_IR.ir[gas][TG_IR.ir['Ts']>TG_IR.info['dry_temp']])
        tot_mol=(tot_area-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]
        peaks['area'][gas.upper()]=tot_area
        peaks['mmol'][gas.upper()]=tot_mol
        peaks['mmol_per_mg'][gas.upper()]=tot_mol/TG_IR.info[ref_mass]
        

        
        if y_axis=='rel':
            FTIR.update(FTIR[gas]/tot_area*tot_mol)
        
        #initial guesses
        init_centers=np.array(centers[gas].dropna().astype(float))+init_offs[0]
        num_curves=len(init_centers)
        max_heigth=max(FTIR[gas])
        init_heights=(BOUNDS.getfloat('height_0')+init_offs[1])*max_heigth*np.ones(num_curves)
        init_hwhms=(BOUNDS.getfloat('hwhm_0')+init_offs[2])*max_hwhm*np.ones(num_curves)
        
        #lower bounds for fit
        min_heights=0*np.ones(num_curves)
        min_hwhms=0*np.ones(num_curves)
        min_centers=init_centers-tol_center
        
        #upper bounds for fit
        max_heights=max_heigth*np.ones(num_curves)
        max_hwhms=max_hwhm*np.ones(num_curves)
        max_centers=init_centers+tol_center
        
        #print(min_centers,init_centers,max_centers)
        params['init_center']=init_centers
        params['init_height']=init_heights
        params['init_hwhm']=init_hwhms
        params['min_center']=min_centers
        params['min_height']=min_heights
        params['min_hwhm']=min_hwhms
        params['max_center']=max_centers
        params['max_height']=max_heights
        params['max_hwhm']=max_hwhms
        #assumptions by figueiredo
        a=True
        if (gas=='co' and a==True):
            try:
                loc=labels[gas][labels[gas]=='anhydrides'].index[0]

                init_heights[loc]=peaks['height']['anhydrides_CO2']/TG_IR.linreg['slope']['co2']*TG_IR.linreg['slope']['co']
                min_heights[loc]=init_heights[loc]*.9
                max_heights[loc]=init_heights[loc]*1.1

                init_hwhms[loc]=peaks['hwhm']['anhydrides_CO2']
                min_hwhms[loc]=init_hwhms[loc]*.9
                max_hwhms[loc]=init_hwhms[loc]*1.1
            except:
                pass
            
            
        #guesses and bounds
        init_guess=np.concatenate((init_heights,init_centers,init_hwhms))
        min_param=np.concatenate((min_heights,min_centers,min_hwhms))
        max_param=np.concatenate((max_heights,max_centers,max_hwhms))
    
        #actual fitting
        x=FTIR['Ts']
        try:
            popt,pcov=sp.optimize.curve_fit(func,x,FTIR[gas],p0=init_guess,bounds=(min_param,max_param))
        except:
            print('Failed to fit {} signal'.format(gas.upper()))
            break
        #return values
        for i in range(num_curves):
            group=labels[gas][i]+'_'+gas.upper()
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
        profiles['Ts']=x
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
            fitting.xaxis.set_ticks(np.arange(0, 1000, 50))
            
            #plotting of fit
            fitting.plot(x,data,label='data',lw=2,zorder=num_curves+1)#,ls='',marker='x',markevery=2,c='cyan')
            fitting.plot(x,fit,label='fit',lw=2,zorder=num_curves+2)
        for i in range(0,num_curves):
            y=gaussian(x,popt[i],popt[i+num_curves],popt[i+2*num_curves])
            profiles[labels[gas][i]]=y
            if plot:
                fitting.text(popt[num_curves+i],popt[i],labels[gas][i],zorder=num_curves+3+i)
                fitting.plot(x,y,linestyle='dashed',zorder=i)
        if plot:
            fitting.legend()
            fitting.set_xlabel('{} {} ${}$'.format(UNITS['ts'], SEP, UNITS['ts']))
            if y_axis=='orig':
                fitting.set_ylabel('{} {} ${}$'.format(get_label(gas), SEP, UNITS['ir']))
            elif y_axis=='rel':
                fitting.set_ylabel('{} {} ${}\,{}^{{-1}}\,{}^{{-1}}$'.format(get_label(gas), SEP, UNITS['molar_amount'], UNITS['mass'], UNITS['time']))

            #mark center on x-axis
            fitting.scatter(popt[num_curves:2*num_curves],np.zeros(num_curves),marker=7,color='k',s=100,zorder=num_curves+3)

            #plotting of absolute difference
            abs_max=0.05*max(data)
            
            error.text(0,abs_max,'SQERR: {:.2e}'.format(sumsqerr[gas][TG_IR.info['name']]))#,'SQERR: '+'%.2E'% Decimal(sumsqerr[gas][TG_IR.info['name']]),va='bottom')
            error.plot(x,diff)
            error.hlines(0,min(x),max(x),ls='dashed')
            error.set_xlabel('{} {} ${}$'.format(UNITS['ts'], SEP, UNITS['ts']))
            error.set_ylabel('error {} ${}$'.format(SEP, UNITS['ir']))
            error.set_ylim(-abs_max,abs_max)
            plt.show()
            fig.savefig(TG_IR.info['name']+'_'+gas+'.png', bbox_inches='tight', dpi=DPI)
        if save:
            try:
                with pd.ExcelWriter(TG_IR.info['name']+y_axis+'.xlsx',engine='openpyxl', mode='a') as writer:
                    profiles.to_excel(writer,sheet_name=gas)
                    params.to_excel(writer,sheet_name=gas+'_param')
            except:
                with pd.ExcelWriter(TG_IR.info['name']+y_axis+'.xlsx',engine='openpyxl') as writer:
                    profiles.to_excel(writer,sheet_name=gas)
                    params.to_excel(writer,sheet_name=gas+'_param')

    if save:                
        with pd.ExcelWriter(TG_IR.info['name']+y_axis+'.xlsx',engine='openpyxl', mode='a') as writer:
                peaks.astype(float).to_excel(writer,sheet_name='summary')
    return peaks.astype(float),sumsqerr