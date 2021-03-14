import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import scipy as sp
import copy

from .config import PATHS, DPI, LABELS, SEP, UNITS, PARAMS, MOLAR_MASS, SAVGOL
from .input_output.general import time

WINDOW_LENGTH=SAVGOL.getint('window_length')
POLYORDER=SAVGOL.getint('POLYORDER')

def get_label(key):
    "get labels to put in plots"
    if key in LABELS:
        return LABELS[key]
    try:
        if int(key)in LABELS:
            return LABELS[int(key)]
    except:
        return str(key)
    
def ylim_auto(x, y, xlim):
    "truncate x and y according to xlim"
    x_min = xlim[0]
    x_max = xlim[1]
    if pd.isnull(xlim[0]): x_min = x.min() 
    if pd.isnull(xlim[1]): x_max = x.max() 
    x = x[(x >= x_min) & (x <= x_max)]  
    y = y[x.index]
    ylim = [None, None]   # reset ylim
    
    return x, y, ylim

def plot_TGA(TG_IR, plot, save=False, x_axis='sample_temp', y_axis='orig', ylim = 'auto', xlim=[None,None], legend=True):
    "plot TG data"
    
    #setting up plot
    fig, TGA=plt.subplots()
    
    fig.subplots_adjust(right=.8)
    
    DTG=TGA.twinx()
    
    x = copy.deepcopy(TG_IR.tga[x_axis])
    
    #if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis

    
    # adjusting y data and setting axis labels according to y_axis
    if y_axis=='rel':
        y=(TG_IR.tga[plot]/TG_IR.info[TG_IR.info['reference_mass']])*100
        yDTG=TG_IR.tga['dtg']*60/TG_IR.info[TG_IR.info['reference_mass']]*100
        ylabelDTG=r'{} {} $ \%\,min^{{-1}}$'.format(PARAMS['dtg'],SEP)
        if plot=='sample_mass':
            ylabel='{} {} $\%$'.format(PARAMS['sample_mass'],SEP)
        elif plot=='heat_flow':
            ylabel='{} {} $ {}\,{}^{{-1}}$'.format(PARAMS['heat_flow'],SEP,UNITS['heat_flow'],UNITS['sample_mass'])
        
    elif y_axis=='orig':
        y=TG_IR.tga[plot]
        yDTG=TG_IR.tga['dtg']*60   # turning dtg from mg/s in mg/min
        ylabelDTG='{} {} ${}\,min^{{-1}}$'.format(PARAMS['dtg'],SEP,UNITS['sample_mass'])
        ylabel='{} {} ${}$'.format(PARAMS[plot],SEP,UNITS[plot])

    TGA.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    
    # adjusting x data if x_axis == time and constructing y-axis for temperature
    if x_axis=='time':
        x=x/60
        temp=TGA.twinx()
        temp.plot(x,TG_IR.tga['sample_temp'],label='{} {} ${}$'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']),ls='dashed',color='black')
        temp.spines['right'].set_position(('axes',1.15))
        temp.set_ylabel('{} {} ${}$'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        if legend:
            temp.legend(loc=1)
    
    if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis
        x, y, ylim = ylim_auto(x, y, xlim)    
        x, yDTG, ylim = ylim_auto(x, yDTG, xlim)    
    
    # actual plotting    
    gTGA,=TGA.plot(x,y,'r-',label=ylabel)
    gDTG,=DTG.plot(x,yDTG,'b--',label='DTG')
    
    TGA.set_ylabel(ylabel)
    TGA.set_ylim(ylim)
    TGA.set_xlim(xlim)
    DTG.set_ylabel(ylabelDTG)
        
    TGA.yaxis.label.set_color(gTGA.get_color())
    DTG.yaxis.label.set_color(gDTG.get_color())
    
    TGA.xaxis.set_minor_locator(ticker.AutoMinorLocator())  # switch on minor ticks on each axis
    TGA.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    DTG.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    
    TGA.set_title('{}, {} = {:.2f} ${}$'.format(TG_IR.info['alias'], TG_IR.info['reference_mass'], TG_IR.info[TG_IR.info['reference_mass']],UNITS['sample_mass']))
    plt.show()
    
    
    if save:
        path_plots_tga=os.path.join(PATHS['dir_plots'],'TGA')
        if os.path.exists(path_plots_tga)==False:
            os.makedirs(path_plots_tga)
        fig.savefig(os.path.join(path_plots_tga,'{}_TG_{}.png'.format(TG_IR.info['name'],y_axis)), bbox_inches='tight',dpi=DPI)

def plot_FTIR(TG_IR, save=False, gases=[], x_axis='sample_temp', y_axis='orig', xlim=[None,None], legend=True):
    "plot IR data"
    gases=set([gas.upper() for gas in gases])
    colors =plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    x=copy.deepcopy(TG_IR.ir[x_axis])
    
    # catching possible input errors
    try:
        calibrated = set(TG_IR.linreg.index)
    except:
        calibrated = set()
    on_axis = set(TG_IR.info['gases'])
    if len(gases) == 0:
        if y_axis=='rel':
            intersection = calibrated & on_axis
            if len(on_axis - calibrated)!=0:
                print('{} not calibrated. Proceeding with {}.'.format(' and '.join([gas for gas in list(on_axis - calibrated)]),
                                                                      ' and '.join([gas for gas in intersection])))
            gases = intersection
        elif y_axis=='orig':
            gases=on_axis
    else:
        if y_axis=='rel':
            gases=set(gases)
            intersection = calibrated & on_axis & gases
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
        
    #setup figure 
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
    elif y_axis=='rel':
        graphs[0].set_ylabel('${}\,{}^{{-1}}\,{}^{{-1}}$'.format(UNITS['molar_amount'],UNITS['sample_mass'],UNITS['time']))

    # plot data and append secondary, third... y-axis on right side if necessary
    for i,gas in enumerate(gases):
        if y_axis=='orig':
            y=TG_IR.ir[gas]
            
            if i>0:
                graphs.append(graphs[0].twinx())
                graphs[i].spines['right'].set_position(('axes',1+(i-1)*.1))
            graphs[i].plot(x,y, color=colors[i])
            graphs[i].set_ylabel('{} {} {}'.format(get_label(gas),SEP,UNITS['ir']))
            graphs[i].yaxis.label.set_color(colors[i])
            graphs[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())      # switch on minor ticks on each axis
            
        elif y_axis=='rel':
            tot_area=np.sum(TG_IR.ir[gas])
            tot_mol=(tot_area-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]/TG_IR.info[TG_IR.info['reference_mass']]
            y=TG_IR.ir[gas]/tot_area*tot_mol
            graphs[0].plot(x,y,label=get_label(gas))

    if legend and y_axis == 'rel':
        fig.legend()
    
    graphs[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())  # switch on minor ticks on each axis
    graphs[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        
    graphs[0].set_title('{}, {}: {:.2f}$\,{}$'.format(TG_IR.info['alias'], TG_IR.info['reference_mass'], TG_IR.info[TG_IR.info['reference_mass']],UNITS['sample_mass']))
    graphs[0].set_xlim(xlim)
    plt.show()
    
    if save:
        path_plots_ir=os.path.join(PATHS['dir_plots'],'IR')
        if os.path.exists(path_plots_ir)==False:
            os.makedirs(path_plots_ir)
        fig.savefig(os.path.join(path_plots_ir,'{}_IR_{}.png'.format(TG_IR.info['name'],y_axis)), bbox_inches='tight',dpi=DPI)
        
def FTIR_to_DTG(TG_IR, x_axis='sample_temp', save=False, gases=[], legend=True, y_axis=None, xlim=[None,None]):
    "reconstructing DTG from calibrated IR data"
    gases_temp=set([gas.upper() for gas in gases])
    
    # catching possible input errors
    try:
        calibrated = set(TG_IR.linreg.index)
    except:
        calibrated = set()
    on_axis=set(TG_IR.info['gases'])
    if len(gases) == 0:
        intersection=calibrated &  on_axis
        if calibrated != on_axis:
            print('{} not calibrated. Proceeding with {}.'.format(' and '.join([gas for gas in list(on_axis - calibrated)]),' and '.join([gas for gas in intersection])))
        gases_temp=list(intersection)
    
    else:
        intersection=calibrated &  on_axis & gases_temp
        if len(gases_temp-calibrated)!=0:
            print('{} not calibrated.'.format(' and '.join([gas for gas in (gases_temp - calibrated)])))
        if len(intersection)!=0:
            print('Proceeding with {}.'.format(' and '.join([gas for gas in intersection]))) 
            gases_temp=list(intersection)
            gases_temp.sort(key = lambda i: gases.index(i) if i in gases else len(gases))
        else:
            print('None of supplied gases was found in IR data and calibrated.')
            return
    
    # setup IR data and calculating DTG 
    gases=gases_temp
    data=pd.merge(TG_IR.tga,TG_IR.ir,how='left',on=['time','sample_temp']).dropna()
    DTG=-sp.signal.savgol_filter(data['sample_mass'], 13, 3, deriv=1)
    
    x=data[x_axis]
    y=np.zeros((len(gases),len(TG_IR.ir)))
    out=pd.DataFrame()
    out['time']=data['time']
    out['sample_temp']=data['sample_temp']
    cumul=np.zeros(len(TG_IR.ir))
    
    # calibrating IR data
    for i,gas in enumerate(gases):
        tot_area=np.sum(TG_IR.ir[gas])
        tot_mass=(tot_area-TG_IR.linreg['intercept'][gas])/TG_IR.linreg['slope'][gas]*MOLAR_MASS.getfloat(gas)
        y[i][:]=TG_IR.ir[gas]/tot_area*tot_mass
        out[gas]=y[i][:]
        cumul+=y[i][:]
   
    # setup figure
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(8, 1)
    stack = fig.add_subplot(gs[:-1, 0])
    stack.set_title('{}, {} = {:.2f} ${}$'.format(TG_IR.info['alias'], TG_IR.info['reference_mass'], TG_IR.info[TG_IR.info['reference_mass']],UNITS['sample_mass']))
    error = fig.add_subplot(gs[-1,0],sharex=stack)

    stack.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    error.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    
    # actual plotting
    if x_axis=='time':
        x=x/60
        temp=stack.twinx()
        temp.plot(x,data['sample_temp'],ls='dashed',color='black',label='T')
        temp.set_ylabel('{} {} ${}$'.format(PARAMS['sample_temp'],SEP,UNITS['sample_temp']))
        if legend:
            temp.legend()
        
    stack.stackplot(x,y,labels=[get_label(gas) for gas in gases])
    stack.plot(x,DTG,label=PARAMS['dtg'])
    stack.set_ylabel('{}, {} {} ${}\,{}^{{-1}}$'.format(PARAMS['dtg'],', '.join([get_label(gas) for gas in gases]),SEP,UNITS['sample_mass'],UNITS['time']))
    
    stack.xaxis.set_minor_locator(ticker.AutoMinorLocator())  # switch on minor ticks on each axis
    stack.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    
    if legend:
        stack.legend()
    
    # plot error of reconstruction    
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
        
def plots(TG_IR_objs, plot, x_axis='sample_temp', y_axis='orig', ylim = 'auto', xlim=[None,None], gas=None, save=False, legend=True, reference_mass='reference_mass'):
    "overlay plots from different objects"
    
    # setting up axis-labels and catching possible input errors
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
            if gas not in TG_IR_objs[0].ir.columns:
                print('{} was not found in IR data.'.format(gas))
                return

        # just to see if supplied gas is calibrated or not
        calibrated = set()
        for TG_IR in TG_IR_objs:
            try:
                calibrated.update(set(TG_IR.linreg.index))
            except:
                if y_axis=='rel':
                    print('{} is not calibrated for {}'.format(gas, TG_IR.info["name"]))
          
        if y_axis=='rel':
            if (calibrated == set()):
                print('{} is not calibrated. Change y_axis to \'orig\' and rerun command.'.format(gas))
                return
    fig,ax=plt.subplots()
    ax.set_xlabel('{} {} ${}$'.format(PARAMS[x_axis.lower()],SEP,UNITS[x_axis.lower()]))
    if plot!='IR':
        if y_axis=='orig':
            ax.set_ylabel('{} {} ${}$'.format(PARAMS[ylabel],SEP,UNITS[ylabel]))
        elif y_axis=='rel':
            if plot == 'DTG':
                ax.set_ylabel('{} {} $\%\,min^{{-1}}$'.format(PARAMS[ylabel],SEP))
            else:
                ax.set_ylabel('{} {} $\%$'.format(PARAMS[ylabel],SEP))
    elif plot=='IR':
        if y_axis=='orig':
            ax.set_ylabel('{} {} ${}$'.format(get_label(gas),SEP,UNITS[ylabel]))
        elif y_axis=='rel':
            ax.set_ylabel('{} {} ${}\,{}^{{-1}}$'.format(get_label(gas),SEP,UNITS['molar_amount'],UNITS['sample_mass']))     
            
    # actual plotting       
    for obj in TG_IR_objs:
        if reference_mass=='reference_mass':
            ref_mass=obj.info[obj.info[reference_mass]]
        else:
            ref_mass=obj.info[reference_mass]
            
        if plot=='TG':
            x=copy.deepcopy(obj.tga[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y=obj.tga['sample_mass']
            elif y_axis=='rel':
                y=100*obj.tga['sample_mass']/ref_mass
            if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(x,y,label = '{}, {}: {:.2f}$\,{}$'.format(obj.info['alias'], obj.info['reference_mass'], obj.info[obj.info['reference_mass']],UNITS['sample_mass']))
        if plot=='DTG':
            x=copy.deepcopy(obj.tga[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y = obj.tga['dtg']*60
            elif y_axis=='rel':
                y = obj.tga['dtg']*60 / ref_mass * 100
            if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(x,y,label = '{}, {}: {:.2f}$\,{}$'.format(obj.info['alias'], obj.info['reference_mass'], obj.info[obj.info['reference_mass']],UNITS['sample_mass']))
        if plot=='heat_flow':
            x=copy.deepcopy(obj.tga[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y=obj.tga['heat_flow']
            elif y_axis=='rel':
                y=obj.tga['heat_flow']/ref_mass
            if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis
                x, y, ylim_temp = ylim_auto(x, y, xlim)
            ax.plot(x,y,label = '{}, {}: {:.2f}$\,{}$'.format(obj.info['alias'], obj.info['reference_mass'], obj.info[obj.info['reference_mass']],UNITS['sample_mass']))
        if plot=='IR':
            x=copy.deepcopy(obj.ir[x_axis])
            if x_axis=='time':
                x/=60
            if y_axis=='orig':
                y = obj.ir[gas]
                if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis
                    x, y, ylim_temp = ylim_auto(x, y, xlim)
                ax.plot(x,y,label = '{}, {}: {:.2f}$\,{}$'.format(obj.info['alias'], obj.info['reference_mass'], obj.info[obj.info['reference_mass']],UNITS['sample_mass']))
            elif y_axis=='rel':
                y = obj.ir[gas] / obj.linreg['slope'][gas] / ref_mass
                if (ylim == 'auto'):   # only select relevant range of x data, to auto-scale the y axis
                    x, y, ylim_temp = ylim_auto(x, y, xlim)
                ax.plot(x,y,label = '{}, {}: {:.2f}$\,{}$'.format(obj.info['alias'], obj.info['reference_mass'], obj.info[obj.info['reference_mass']],UNITS['sample_mass']))
    
    if (ylim == 'auto'):   # reset ylim to [None,None]
        ylim = ylim_temp
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    
    if legend:
        ax.legend()
        
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())  # switch on minor ticks on each axis
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    plt.show()
    
    if save:
        path_plots=PATHS['dir_plots']
        if os.path.exists(path_plots)==False:
            os.makedirs(path_plots)
        if (gas == None):
            gas = ''
        fig.savefig(os.path.join(path_plots,'_'.join([time(),plot,gas,y_axis]))+'.png', bbox_inches='tight',dpi=DPI)
