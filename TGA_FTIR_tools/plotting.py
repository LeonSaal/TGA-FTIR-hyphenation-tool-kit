import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import copy

from .config import PATHS, DPI,LABELS, SEP, UNITS, PARAMS, MOLAR_MASS, SAVGOL
from .input_output.general import time


def get_label(key):
    if key in LABELS:
        return LABELS[key]
    try:
        if int(key)in LABELS:
            return LABELS[int(key)]
    except:
        return str(key)

def plot_TGA(TG_IR,plot,save=False,x_axis='sample_temp',y_axis='orig',ylim=[None,None],xlim=[None,None],legend=True):
    fig, TGA=plt.subplots()
    
    fig.subplots_adjust(right=.8)
    
    DTG=TGA.twinx()
    
    x=copy.deepcopy(TG_IR.tga[x_axis])
    if y_axis=='rel':
        y=(TG_IR.tga[plot]/TG_IR.tga['sample_mass'][0])*100
        yDTG=TG_IR.tga['dtg']*60/TG_IR.tga['sample_mass'][0]
        ylabelDTG=r'{} {} $ \%\,min^{{-1}}$'.format(PARAMS['dtg'],SEP)
        if plot=='sample_mass':
            ylabel='{} {} $\%$'.format(PARAMS['sample_mass'],SEP)
        elif plot=='heat_flow':
            ylabel='{} {} $ {}\,{}^{{-1}}$'.format(PARAMS['heat_flow'],SEP,UNITS['heat_flow'],UNITS['sample_mass'])
        
    elif y_axis=='orig':
        y=TG_IR.tga[plot]
        yDTG=TG_IR.tga['dtg']*60
        ylabelDTG='{} {} ${}\,min^{{-1}}$'.format(PARAMS['dtg'],SEP,UNITS['sample_mass'])
        ylabel='{} {} ${}$'.format(PARAMS[plot],SEP,UNITS[plot])

    TGA.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    if x_axis=='time':
        x=x/60
        temp=TGA.twinx()
        temp.plot(x,TG_IR.tga['sample_temp'],label='{} {} ${}$'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']),ls='dashed',color='black')
        temp.spines['right'].set_position(('axes',1.15))
        temp.set_ylabel('{} {} ${}$'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        if legend:
            temp.legend(loc=1)
        
    gTGA,=TGA.plot(x,y,'r-',label=ylabel)
    gDTG,=DTG.plot(x,yDTG,'b--',label='DTG')
    
    TGA.set_ylabel(ylabel)
    TGA.set_ylim(ylim)
    TGA.set_xlim(xlim)
    DTG.set_ylabel(ylabelDTG)
        
    TGA.yaxis.label.set_color(gTGA.get_color())
    DTG.yaxis.label.set_color(gDTG.get_color())
    
    TGA.set_title('{}, {:.2f} ${}$'.format(TG_IR.info['alias'],TG_IR.info['initial_mass'],UNITS['sample_mass']))
    plt.show()
    
    
    if save:
        path_plots_tga=os.path.join(PATHS['dir_plots'],'TGA')
        if os.path.exists(path_plots_tga)==False:
            os.makedirs(path_plots_tga)
        fig.savefig(os.path.join(path_plots_tga,'{}_TG_{}.png'.format(TG_IR.info['name'],y_axis)), bbox_inches='tight',dpi=DPI)

def plot_FTIR(TG_IR,save=False,gases=[],x_axis='sample_temp',y_axis='orig',xlim=[None,None],legend=True):
    gases=set([gas.upper() for gas in gases])
    colors =plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    x=copy.deepcopy(TG_IR.ir[x_axis])
    try:
        calibrated=set(TG_IR.linreg.index)
    except:
        calibrated=set()
    on_axis=set(TG_IR.info['gases'])
    if len(gases) == 0:
        if y_axis=='rel':
            intersection=calibrated &  on_axis
            if len(on_axis - calibrated)!=0:
                print('{} not calibrated. Proceeding with {}.'.format(' and '.join([gas for gas in list(on_axis - calibrated)]),' and '.join([gas for gas in intersection])))
            gases=intersection
        elif y_axis=='orig':
            gases=on_axis
    else:
        if y_axis=='rel':
            gases=set(gases)
            intersection=calibrated &  on_axis & gases
            if len(gases - calibrated)!=0:
                print('{} not calibrated.'.format(' and '.join([gas for gas in (gases - calibrated)])))
            if len(intersection)!=0:
                print('Proceeding with {}.'.format(' and '.join([gas for gas in intersection])))
                gases=intersection
            else:
                print('None of supplied gases was found in IR data and calibrated.')
                return
        elif y_axis=='orig':
            gases=set(gases)
            intersection = on_axis & gases
            if len(gases - on_axis)!=0:
                print('{} not found in IR data.'.format(' and '.join([gas for gas in (gases - on_axis)])))
            if len(intersection)!=0:
                print('Proceeding with {}.'.format(' and '.join([gas for gas in intersection])))
                gases=intersection
            else:
                print('None of supplied gases was found in IR data.')   
                return
    gases=list(gases)        
        
    #setup figure and plot first gas
    graphs=[]
    fig, g0=plt.subplots()
    fig.subplots_adjust(right=.8)
    graphs.append(g0)
    graphs[0].set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    
    if x_axis=='time':
        x=x/60
    if y_axis=='orig':
        graphs[0].set_ylabel('{} {} ${}$'.format(get_label(gases[0]),SEP,UNITS['ir']))
        graphs[0].yaxis.label.set_color(colors[0])
        #graphs[0].plot(x,TG_IR.ir[gases[0]])
    elif y_axis=='rel':
        graphs[0].set_ylabel('${}\,{}^{{-1}}\,{}^{{-1}}$'.format(UNITS['molar_amount'],UNITS['sample_mass'],UNITS['time']))
        #graphs[0].plot(x,TG_IR.ir[gases[0]],label=get_label(gases[0]))

    
    #append secondary, third... y-axis on right side
     
    for i,gas in enumerate(gases):
        if y_axis=='orig':
            y=TG_IR.ir[gas]
            
            if i>0:
                graphs.append(graphs[0].twinx())
                graphs[i].spines['right'].set_position(('axes',1+(i-1)*.1))
            graphs[i].plot(x,y, color=colors[i])
            graphs[i].set_ylabel('{} {} {}'.format(get_label(gas),SEP,UNITS['ir']))
            graphs[i].yaxis.label.set_color(colors[i])
            
        elif y_axis=='rel':
            tot_area=np.sum(TG_IR.ir[gas])
            tot_mol=(tot_area-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]/TG_IR.info[TG_IR.info['reference_mass']]
            y=TG_IR.ir[gas]/tot_area*tot_mol
            graphs[0].plot(x,y,label=get_label(gas))

    if legend: #and y_axis!='orig':
        fig.legend()
        
    graphs[0].set_title('{}, {:.2f} ${}$'.format(TG_IR.info['alias'],TG_IR.info['initial_mass'],UNITS['sample_mass']))
    graphs[0].set_xlim(xlim)
    plt.show()
    
    if save:
        path_plots_ir=os.path.join(PATHS['dir_plots'],'IR')
        if os.path.exists(path_plots_ir)==False:
            os.makedirs(path_plots_ir)
        fig.savefig(os.path.join(path_plots_ir,'{}_IR_{}.png'.format(TG_IR.info['name'],y_axis)), bbox_inches='tight',dpi=DPI)
        
def FTIR_to_DTG(TG_IR,x_axis='sample_temp',save=False,gases=[],legend=True,y_axis=None,xlim=[None,None]):
    gases=set([gas.upper() for gas in gases])
    try:
        calibrated=set(TG_IR.linreg.index)
    except:
        calibrated=set()
    on_axis=set(TG_IR.info['gases'])
    if len(gases) == 0:
        intersection=calibrated &  on_axis
        if calibrated != on_axis:
            print('{} not calibrated. Proceeding with {}.'.format(' and '.join([gas for gas in list(on_axis - calibrated)]),' and '.join([gas for gas in intersection])))
        gases=intersection
    
    else:
        intersection=calibrated &  on_axis & gases
        if len(gases-calibrated)!=0:
            print('{} not calibrated.'.format(' and '.join([gas for gas in (gases - calibrated)])))
        if len(intersection)!=0:
            print('Proceeding with {}.'.format(' and '.join([gas for gas in intersection])))
            gases=intersection
        else:
            print('None of supplied gases was found in IR data and calibrated.')
            return

    data=pd.merge(TG_IR.tga,TG_IR.ir,how='left',on=['time','sample_temp']).dropna()
    DTG=-sp.signal.savgol_filter(data['sample_mass'], 13, 3, deriv=1)
    
    x=data[x_axis]
    y=np.zeros((len(gases),len(TG_IR.ir)))
    out=pd.DataFrame()
    out['time']=data['time']
    out['sample_temp']=data['sample_temp']
    cumul=np.zeros(len(TG_IR.ir))
    for i,gas in enumerate(gases):
        tot_area=np.sum(TG_IR.ir[gas])
        tot_mass=(tot_area-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]*MOLAR_MASS.getfloat(gas)
        y[i][:]=TG_IR.ir[gas]/tot_area*tot_mass
        out[gas]=y[i][:]
        cumul+=y[i][:]
   
    #setup
    fig=plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(8, 1)
    stack= fig.add_subplot(gs[:-1, 0])
    stack.set_title('{}, {:.2f} {}'.format(TG_IR.info['alias'],TG_IR.info['initial_mass'],UNITS['sample_mass']))
    error = fig.add_subplot(gs[-1,0],sharex=stack)

    #plotting of fit
    stack.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    error.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    
    if x_axis=='time':
        x=x/60
        temp=stack.twinx()
        temp.plot(x,data['sample_temp'],ls='dashed',color='black',label='T')
        temp.set_ylabel('{} {} ${}$'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        if legend:
            temp.legend(loc=1)
        
    stack.stackplot(x,y,labels=[get_label(gas) for gas in gases])
    stack.plot(x,DTG,label=PARAMS['dtg'])
    stack.set_ylabel('{}, {} {} ${}\,{}^{{-1}}$'.format(PARAMS['dtg'],', '.join([get_label(gas) for gas in gases]),SEP,UNITS['sample_mass'],UNITS['time']))
    if legend:
        stack.legend()
    error.plot(x,DTG-cumul)
    error.hlines(0,min(x),max(x),ls='dashed')
    error.set_ylabel('$\Delta$ {}'.format(PARAMS['dtg']))
    stack.set_xlim(xlim)
    plt.show()
    if save:
        path_plot_irdtg=os.path.join(PATHS['dir_plots'],'IRDTG')
        if os.path.exists(path_plot_irdtg)==False:
            os.makedirs(path_plot_irdtg)
        fig.savefig(os.path.join(path_plot_irdtg,'{}_IRDTG.png'.format(TG_IR.info['name'])), bbox_inches='tight',dpi=DPI)
        out['dtg']=DTG
        out.to_excel(os.path.join(PATHS['dir_output'],TG_IR.info['name']+'_IRDTG.xlsx'))
        
def plots(TG_IR,plot,x_axis='sample_temp',y_axis='orig',ylim=[None,None],xlim=[None,None],gas=None,save=False,legend=True):
    #out=pd.DataFrame()
    if plot=='TG':
        ylabel='sample_mass'
    else:
        ylabel=plot.lower()
    if plot=='IR':
        if gas==None:
            print('Supply \'gas = \'')
            return
        else:
            gas=gas.upper()
            if gas not in TG_IR[0].ir.columns:
                print('{} was not found in IR data.'.format(gas))
                return
        if y_axis=='rel':
            if gas not in TG_IR[0].linreg.index:
                print('{} is not calibrated.'.format(gas))
                return
    fig,ax=plt.subplots()
    ax.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    if plot!='IR':
        if y_axis=='orig':
            ax.set_ylabel('{} {} ${}$'.format(PARAMS[ylabel],SEP,UNITS[ylabel]))
        elif y_axis=='rel':
            ax.set_ylabel('{} {} $\%$'.format(PARAMS[ylabel],SEP))
    elif plot=='IR':
        if y_axis=='orig':
            ax.set_ylabel('{} {} ${}$'.format(get_label(gas),SEP,UNITS[ylabel]))
        elif y_axis=='rel':
            ax.set_ylabel('{} {} ${}\,{}^{{-1}}$'.format(get_label(gas),SEP,UNITS['molar_amount'],UNITS['sample_mass']))     
        ax.set_ylabel('{} {} {} ${}^{{-1}}\,{}^{{-1}}$'.format(PARAMS[ylabel],SEP,UNITS[ylabel],UNITS['sample_mass'],UNITS['time']))
            
            
    for obj in TG_IR:
        if plot=='TG':
            x=copy.deepcopy(obj.tga[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y=obj.tga['sample_mass']
            elif y_axis=='rel':
                y=obj.tga['sample_mass']/obj.info[obj.info['reference_mass']]
            ax.plot(x,y,label=obj.info['alias'])
        if plot=='heat flow':
            x=copy.deepcopy(obj.tga[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y=obj.tga['heat_flow']
            elif y_axis=='rel':
                y=obj.tga['heat_flow']/obj.info[obj.info['reference_mass']]
            ax.plot(x,y,label=obj.info['alias'])
        if plot=='IR':
            x=copy.deepcopy(obj.ir[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y=obj.ir[gas]
                ax.plot(x,y,label=obj.info['alias'])#+', dry mass: '+str(round(obj.info['reference_mass'],2))+' mg'
            elif y_axis=='rel':
                y=obj.ir[gas]/obj.linreg['slope'][gas]/obj.info[obj.info['reference_mass']]
                ax.plot(x,y,label='{}, {:.2f} {}'.format(obj.info['alias'],obj.info['initial_mass'],UNITS['sample_mass']))
            
        
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    if legend:
        ax.legend()
    plt.show()
    
    if save:
        path_plots=PATHS['dir_plots']
        if os.path.exists(path_plots)==False:
            os.makedirs(path_plots)
        fig.savefig(os.path.join(path_plots,'_'.join([time(),plot,gas,y_axis]))+'.png', bbox_inches='tight',dpi=DPI)
