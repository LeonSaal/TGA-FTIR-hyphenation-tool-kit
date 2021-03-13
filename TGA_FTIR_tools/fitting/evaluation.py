import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from ..config import PATHS, DPI
from ..input_output.general import time


def summarize(path, select_groups = [], condense_results = True):
# import robustness results from robustness_in_mmol_per_mg.xlsx in 'path'
    robustness_summary = pd.read_excel(path + '/robustness_in_mmol_per_mg.xlsx', 
                                       sheet_name = 'summary',
                                       header = 0,
                                       index_col = 1)
    # recreate  multiindex
    for i in range(len(robustness_summary['samples'])-1):
        if pd.isnull(robustness_summary.iloc[i+1, 0]):
            robustness_summary.iloc[i+1, 0] = robustness_summary.iloc[i, 0]
    robustness_summary = robustness_summary.set_index(['samples', robustness_summary.index])
        
# Transpose the df, select relevant columns and rename them
    try:
        df = robustness_summary.T.loc[:, (slice(None), ['mean', 'stddev', 'limits'])].rename(columns={'mean': 'mmol_per_mg', 'limits': 'label'})
    except:
        df = robustness_summary.T.loc[:, (slice(None), ['mean', 'stddev'])].rename(columns={'mean': 'mmol_per_mg'})
        
# condense results based on groups, independent from gases
    if condense_results == True:
        # extract gases and groups from index
        gases = list(set([re.split('_| ',group)[-1] for group in df.index]))
        groups = list(set([re.split('_| ',group)[0] for group in df.index if group not in gases]))
        # select samples from df, keeping their sequence
        samples = df.columns.get_level_values('samples')
        seen = set()
        samples = [x for x in samples if not (x in seen or seen.add(x))]
        # initialize DataFrame for results
        columns = df.columns
        summarized = pd.DataFrame(index = groups, columns = columns)
        for group in groups:
            group_set = df.loc[df.index.map(lambda x: x.startswith(group))].copy()
            
            # set values < LOD or < LOQ to 0.0
            for sample in samples:
                try:
                    mask = pd.notnull(group_set.loc[:,(sample, 'label')])
                    group_set.loc[mask,(sample, ['mmol_per_mg', 'stddev'])] = 0.0   
                except:
                    pass
            
            if (group == 'anhydrides'):    
                summarized.loc[group,(slice(None), 'mmol_per_mg')] = group_set.loc[:,(slice(None), 'mmol_per_mg')].mean(axis=0)
                # error propagation for the mean of stddev
                group_set.loc[:,(slice(None), 'stddev')] = group_set.loc[:,(slice(None), 'stddev')]**2
                summarized.loc[group,(slice(None), 'stddev')] = group_set.loc[:,(slice(None), 'stddev')].sum(axis=0)**0.5 / len(group_set.loc[:,(slice(None), 'stddev')])
            else:
                summarized.loc[group,(slice(None), 'mmol_per_mg')] = group_set.loc[:,(slice(None), 'mmol_per_mg')].sum(axis=0)
                # error propagation for the sum of stddev
                group_set.loc[:,(slice(None), 'stddev')] = group_set.loc[:,(slice(None), 'stddev')]**2
                summarized.loc[group,(slice(None), 'stddev')] = group_set.loc[:,(slice(None), 'stddev')].sum(axis=0)**0.5
           
            # create label for condensed group
            for sample in samples:
                if (summarized.loc[group,(sample, 'mmol_per_mg')] == 0.0):
                    try:
                        mask = pd.notnull(group_set.loc[:,(sample, 'label')])   # select rows with labels in group_set
                        label_set = set(group_set.loc[mask,(sample, 'label')])  # create a set of these labels
                        if (len(label_set) == 2):   # if there are two items in the set
                            label_text = '< LOQ'
                        elif (len(label_set) == 1): # if there is only one item
                            label_text = str(label_set)
                        summarized.loc[group,(sample, 'label')] = label_text   # set the label according to label_text
                    except:
                        pass
    else:
        # Values are not manipulated to preserve details including the labels
        summarized = df
    
    if (len(select_groups) > 0):
        summarized = summarized.loc[select_groups,:]
        
    return summarized


def concatenate(*dfs):
    df = pd.concat(list(dfs), axis=1, sort=False)
    
    # extract samples from column level
    samples = df.columns.get_level_values('samples')
    seen = set()
    samples = [x for x in samples if not (x in seen or seen.add(x))]
    # set NaN in 'mmol_per_mg' to n.d. in 'label'
    for sample in samples:
        try:
            df.loc[df.loc[:,(sample, 'mmol_per_mg')].isnull(), (sample, 'label')] = 'n.d.'
            df.loc[df.loc[:,(sample, 'label')].isnull(), (sample, 'label')] = ''
        except:
            pass
    
    return df


def bar_plot_results(df, show_groups = [], group_by = 'samples', save = False):
    if (len(show_groups) > 0):
        df = df.reindex(index = show_groups)   # select and order df index/groups
    
    # select samples from df, keeping their sequence
    samples = df.columns.get_level_values('samples')
    seen = set()
    samples = [x for x in samples if not (x in seen or seen.add(x))]
    
    # set labeled values to 0.0
    for sample in samples:
        df.loc[df.loc[:, (sample, 'label')] != '', (sample, 'mmol_per_mg')] = 0.0
        df.loc[df.loc[:, (sample, 'label')] != '', (sample, 'stddev')] = 0.0
    
    fig, ax = plt.subplots()
    
    if group_by == 'samples':
        x = samples   # samples as x-ticks
        groups = df.index   # groups for each sample
    elif group_by == 'groups':
        x = df.index   # groups as x-ticks
        groups = samples   # samples for each group
    
    x_num = np.arange(len(x))   # define x-ticks nummerically
    width = 1 / (len(groups)+1)   # calculate width of a bar; +1 is the gap to bars of the next sample
    
    # find maximum y value in the data, for adjusting of labels on bars
    y_max = 0.0
    for group in groups:
        if group_by == 'samples':
            max_group = df.loc[group, (slice(None), 'mmol_per_mg')].max()
        elif group_by == 'groups':
            max_group = df.loc[x, (group, 'mmol_per_mg')].max()
        
        if (max_group > y_max): 
            y_max = max_group
    
    # rotate x-axis labels and labels on bars dependent on length of data to plot
    y_lab_rotation = 'horizontal'
    y_lab_space = y_max / 100.0
    if ((2 * len(x) + len(groups)) > 15):
        y_lab_rotation = 'vertical'
        y_lab_space = y_max / 50.0
        
    if (len(x) > 5):
        plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')   # rotate x-axis labels
        
    # loop through the bars (defines as groups)
    for i,group in enumerate(groups):        
        # define y values according to 
        if group_by == 'samples':
            y = df.loc[group, (slice(None), 'mmol_per_mg')]
            y_err = df.loc[group, (slice(None), 'stddev')]
            try:
                y_lab = df.loc[group, (slice(None), 'label')]
            except:
                pass
        elif group_by == 'groups':
            y = df.loc[x, (group, 'mmol_per_mg')]
            y_err = df.loc[x, (group, 'stddev')]
            try:
                y_lab = df.loc[x, (group, 'label')]
            except:
                pass
        
        x_i = (x_num - width*(len(groups)-1)/2 + i*width)
                
        ax.bar(x_i, y, yerr = y_err, capsize=5,
               width = width*0.9, label = group, zorder=10)   # color = c,
        
        # add labels to bars
        for j in range(len(x_i)):
            x_j = x_i[j]
            y_j = y.iloc[j] + y_lab_space
            text = y_lab.iloc[j]
            ax.text(x_j, y_j, text, horizontalalignment = 'center', rotation = y_lab_rotation)
    
    ax.legend()
    ax.yaxis.grid(True)
    
    ax.set(title = 'summary plot with errors from robustness testing')
    ax.set_ylabel('surface oxygen groups in $mmol\,mg^{-1}$')
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x)
    
    fig.tight_layout()
    plt.show()
    
    if save:
        path_plots_eval = os.path.join(PATHS['dir_plots'],'EVALUATION')
        if os.path.exists(path_plots_eval)==False:
            os.makedirs(path_plots_eval)
        fig.savefig(os.path.join(path_plots_eval,'{}_{}_by_{}.png'.format(time(), samples, group_by)), bbox_inches='tight',dpi=DPI)
    
    return