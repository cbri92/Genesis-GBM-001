# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:17:46 2020

@author: Caterina Brighi
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.patches as ptch
import matplotlib.pyplot as plt
import os
from pylab import plot, xlim, legend, setp

#%%Functions

# function for setting the colors of the box plots pairs
def setBoxColors(bp):
    
    '''This function sets the colors of the box plot pairs as blue and red.'''
    
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['medians'][1], color='red')

def box_plot(dataframe, column_list, ylabel, color, ymin='Auto', ymax='Auto', tick_min='Auto', tick_max='Auto'):
    
    '''This function generates a box plot of the column_list of the dataframe specified.
    You need to provide ylabel as a string, y_min, y_max, tick_min, tick_max as floats and color as a list of strings.
    The default ymin, ymax, tick_min and tick_max is 'Auto'. '''

    ax.set_ylabel(ylabel, fontsize=15, wrap=True)
    if ymin == 'Auto' and ymax == 'Auto':
        ax.set_ylim(auto=True)
    else:
        ax.set_ylim(bottom=ymin, top=ymax)
        
    if tick_min == 'Auto' and tick_max == 'Auto':
        bp = dataframe.boxplot(column=column_list, fontsize=15, return_type='both', patch_artist=True, grid=False)
        for i, f in zip(range(len(column_list)), color):     
            bp[0].findobj(ptch.Patch)[i].set_facecolor(f)
            bp[0].findobj(ptch.Patch)[i].set_edgecolor("black")
            bp[0].findobj(ptch.Patch)[i].set_linewidth(1.0)
    else:
        ax.yaxis.set_major_locator(plt.MultipleLocator(tick_max))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(tick_min))
        bp = dataframe.boxplot(column=column_list, fontsize=15, return_type='both', patch_artist=True, grid=False)
        for i, f in zip(range(len(column_list)), color):     
            bp[0].findobj(ptch.Patch)[i].set_facecolor(f)
            bp[0].findobj(ptch.Patch)[i].set_edgecolor("black")
            bp[0].findobj(ptch.Patch)[i].set_linewidth(1.0)

#%%Statistical analysis

data_supradir = 'Path to data directory' #Set working directory            
fileDir = data_supradir

Trial_stats_results = pd.ExcelWriter(fileDir +'Trial_statAnalysis_results.xlsx')

#Calculate volumes ratio statistics
volumes_ratio = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='Volumes ratio', index_col='Subject_ID')

#Calculate FET roi on FETPET in FETCT stats
FETroionFETPET = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='FETroi_onFETPET', index_col='Subject_ID')

FETroiSUVmean = FETroionFETPET["Mean intensity [SUV]"]
FETroiSUVmean = pd.DataFrame(FETroiSUVmean)
FETroiSUVmean.rename({'Mean intensity [SUV]': 'FET roi on FETPET SUVmean'}, axis=1, inplace=True)

FETroiSUVmax = FETroionFETPET["Max intensity [SUV]"]
FETroiSUVmax = pd.DataFrame(FETroiSUVmax)
FETroiSUVmax.rename({'Max intensity [SUV]': 'FET roi on FETPET SUVmax'}, axis=1, inplace=True)

#Calculate TER roi on TERPET in FETCT stats
TERroionTERPET = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='TERroi_onTERPET', index_col='Subject_ID')

TERroiSUVmean = TERroionTERPET["Mean intensity [SUV]"]
TERroiSUVmean = pd.DataFrame(TERroiSUVmean)
TERroiSUVmean.rename({'Mean intensity [SUV]': 'TER roi on TERPET SUVmean'}, axis=1, inplace=True)

TERroiSUVmax = TERroionTERPET["Max intensity [SUV]"]
TERroiSUVmax = pd.DataFrame(TERroiSUVmax)
TERroiSUVmax.rename({'Max intensity [SUV]': 'TER roi on TERPET SUVmax'}, axis=1, inplace=True)


#Calculate TER roi on TERPET in TERCT stats
TERroionTERPETorig = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='TERroi_onTERPETorig', index_col='Subject_ID')

TERroiSUVmeanOrig = TERroionTERPETorig["Mean intensity [SUV]"]
TERroiSUVmeanOrig = pd.DataFrame(TERroiSUVmeanOrig)
TERroiSUVmeanOrig.rename({'Mean intensity [SUV]': 'TER roi on TERPETorig SUVmean'}, axis=1, inplace=True)

TERroiSUVmaxOrig = TERroionTERPETorig["Max intensity [SUV]"]
TERroiSUVmaxOrig = pd.DataFrame(TERroiSUVmaxOrig)
TERroiSUVmaxOrig.rename({'Max intensity [SUV]': 'TER roi on TERPETorig SUVmax'}, axis=1, inplace=True)



#Calculate rois TBR stats
TBRval = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='TBR', index_col='Subject_ID')
# TBRval = TBRval.drop(columns=['FET stdTBRmean MI', 'FET stdTBRmean CS', 'TER stdTBRmean MI', 'TER stdTBRmean CS'])

#Group patients in PsP group (FET<2) and Recurrent tumour group (FET>2)

def PsP_vs_RT_classification(name_FET_TBRmean_col, name_FET_TBRmeanSTD_col, method):
    
    '''This function classifies subjects into two groups based on their FET TBRmean value: PsP with FET<2.0 and RT with FET>2.0.
   Generates two dataframes, PsP_df and RT_df, including Subject_ID as index and values of PSMA SUVmean and PSMA SUVmax.'''
   
    FETTBRmean = TBRval.loc[:,name_FET_TBRmean_col]
    FETTBRmeanSTD = TBRval.loc[:,name_FET_TBRmeanSTD_col]
    x=FETTBRmean+FETTBRmeanSTD
    maskMin = x<2.0
    PsP_group = FETTBRmean[maskMin]
    PsP_array = []
    for i in PsP_group.index:
        PsP_array.append([i, TERroiSUVmean.loc[i,'TER roi on TERPET SUVmean'], TERroiSUVmax.loc[i, 'TER roi on TERPET SUVmax']])
    PsP_array = np.array(PsP_array)
    PsP_df = pd.DataFrame(PsP_array, columns=['Subject_ID', 'PsP PSMA SUVmean '+method, 'PsP PSMA SUVmax '+method])  
    PsP_df = PsP_df.set_index('Subject_ID')
    
    maskMax = x>2.0
    RT_group = FETTBRmean[maskMax]
    RT_array = []
    for i in RT_group.index:
        RT_array.append([i, TERroiSUVmean.loc[i,'TER roi on TERPET SUVmean'], TERroiSUVmax.loc[i, 'TER roi on TERPET SUVmax']])
    RT_array = np.array(RT_array)
    RT_df = pd.DataFrame(RT_array, columns=['Subject_ID', 'RT PSMA SUVmean '+method, 'RT PSMA SUVmax '+method])
    RT_df = RT_df.set_index('Subject_ID')
    
    return PsP_df, RT_df

PsP_df_MI, RT_df_MI = PsP_vs_RT_classification('FET TBRmean MI', 'FET stdTBRmean MI', 'MI')
PsP_df_CS, RT_df_CS = PsP_vs_RT_classification('FET TBRmean CS', 'FET stdTBRmean CS', 'CS')

PsPRT_dfs = pd.ExcelWriter(fileDir +'PsPvsRT_results.xlsx') 
PsP_df_MI.to_excel(PsPRT_dfs, sheet_name='PsP_df_MI', index=True)
RT_df_MI.to_excel(PsPRT_dfs, sheet_name='RT_df_MI', index=True)
PsP_df_CS.to_excel(PsPRT_dfs, sheet_name='PsP_df_CS', index=True)
RT_df_CS.to_excel(PsPRT_dfs, sheet_name='RT_df_CS', index=True)
PsPRT_dfs.save()

PsP_df_MI=pd.read_excel(fileDir + 'PsPvsRT_results.xlsx', sheet_name='PsP_df_MI', index_col='Subject_ID')
RT_df_MI=pd.read_excel(fileDir + 'PsPvsRT_results.xlsx', sheet_name='RT_df_MI', index_col='Subject_ID')
PsP_df_CS=pd.read_excel(fileDir + 'PsPvsRT_results.xlsx', sheet_name='PsP_df_CS', index_col='Subject_ID')
RT_df_CS=pd.read_excel(fileDir + 'PsPvsRT_results.xlsx', sheet_name='RT_df_CS', index_col='Subject_ID')

#Calculate FETminCE roi on TERPET in FETCT stats
FETminCEonTERPET = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='FETminCE_onTERPET', index_col='Subject_ID')

FETminCESUVmean = FETminCEonTERPET["Mean intensity [SUV]"]
FETminCESUVmean = pd.DataFrame(FETminCESUVmean)
FETminCESUVmean.rename({'Mean intensity [SUV]': 'FETminCE on TERPET SUVmean'}, axis=1, inplace=True)

FETminCESUVmax = FETminCEonTERPET["Max intensity [SUV]"]
FETminCESUVmax = pd.DataFrame(FETminCESUVmax)
FETminCESUVmax.rename({'Max intensity [SUV]': 'FETminCE on TERPET SUVmax'}, axis=1, inplace=True)

#Calculate FETmaxCE roi on TERPET in FETCT stats
FETmaxCEonTERPET = pd.read_excel(fileDir + 'Group_stats_results.xlsx', sheet_name='FETmaxCE_onTERPET', index_col='Subject_ID')

FETmaxCESUVmean = FETmaxCEonTERPET["Mean intensity [SUV]"]
FETmaxCESUVmean = pd.DataFrame(FETmaxCESUVmean)
FETmaxCESUVmean.rename({'Mean intensity [SUV]': 'FETmaxCE on TERPET SUVmean'}, axis=1, inplace=True)

FETmaxCESUVmax = FETmaxCEonTERPET["Max intensity [SUV]"]
FETmaxCESUVmax = pd.DataFrame(FETmaxCESUVmax)
FETmaxCESUVmax.rename({'Max intensity [SUV]': 'FETmaxCE on TERPET SUVmax'}, axis=1, inplace=True)


#Put all rois stats together into a final df
final = pd.concat([volumes_ratio, FETroiSUVmean, FETroiSUVmax, TERroiSUVmean, TERroiSUVmax, TERroiSUVmeanOrig, TERroiSUVmaxOrig, TBRval,
                   PsP_df_MI, RT_df_MI, PsP_df_CS, RT_df_CS, FETminCESUVmean, FETminCESUVmax, FETmaxCESUVmean, FETmaxCESUVmax], axis=1)
final_stats = final.describe()
final_stats.to_excel(Trial_stats_results, sheet_name='Trial stats', index=True)


#Perform statistical tests


storage = []
listnames = [['FET TBRmean MI', 'TER TBRmean MI'], ['FET TBRmax MI', 'TER TBRmax MI'],['FET TBRmean CS', 'TER TBRmean CS'],
             ['FET TBRmax CS', 'TER TBRmax CS'],['FET TBRmean MI', 'FET TBRmean CS'],['FET TBRmax MI', 'FET TBRmax CS'],
             ['TER TBRmean MI', 'TER TBRmean CS'],['TER TBRmax MI', 'TER TBRmax CS'],
             ['PsP PSMA SUVmean MI', 'RT PSMA SUVmean MI'],
             ['PsP PSMA SUVmax MI', 'RT PSMA SUVmax MI'],
             ['PsP PSMA SUVmean CS', 'RT PSMA SUVmean CS'],
             ['PsP PSMA SUVmax CS', 'RT PSMA SUVmax CS'],
             ['FETminCE on TERPET SUVmean', 'FETmaxCE on TERPET SUVmean'],
             ['FETminCE on TERPET SUVmax', 'FETmaxCE on TERPET SUVmax']]

for i in listnames:
    if i in listnames[8:]:
        t, p = stats.mannwhitneyu(final[i[0]].dropna(), final[i[1]].dropna())
    elif i in listnames[:8]:
        t, p = stats.wilcoxon(final[i[0]], final[i[1]])
    
    if p < 0.0001:
        ast = '****'
    elif p < 0.001:
        ast = '***'
    elif p < 0.01:
        ast = '**'
    elif p < 0.05:
        ast = '*'
    else: 
        ast = 'ns' 

    storage.append([t,p,ast])
    
storage = np.array(storage)


statistical_tests = pd.DataFrame(storage, index=['FET vs TER TBRmean MI', 'FET vs TER TBRmax MI', 'FET vs TER TBRmean CS', 'FET vs TER TBRmax CS',
                                                 'FET TBRmean MI vs CS', 'FET TBRmax MI vs CS', 'TER TBRmean MI vs CS', 'TER TBRmax MI vs CS',
                                                 'PSMA SUVmean PsP vs RT MI', 'PSMA SUVmax PsP vs RT MI',
                                                 'PSMA SUVmean PsP vs RT CS', 'PSMA SUVmax PsP vs RT CS',
                                                 'FETmaxCE vs FETminCE TER SUVmean', 'FETmaxCE vs FETminCE TER SUVmax'], 
                                 columns=['test stats', 'p value', 'significance'])

statistical_tests.to_excel(Trial_stats_results, sheet_name='Statistical tests', index=True)

Trial_stats_results.save() 



#%%Plotting graphs analysis

#Make a directory for the plots
os.mkdir(fileDir +'/Analysis Results Plots')

#Set paths to plot dir
plot_dir = fileDir +'/Analysis Results Plots'


#Box plot tumour volume ratios
fig, ax = plt.subplots()
box_plot(final, ['TER tum volume/FET tum volume'], 'PSMA/FET tumour volumes [a.u.]', 'cyan', ymin=0.0, ymax=1.5, tick_min=0.1, tick_max=0.5)
plt.savefig(plot_dir + '/VolumesRatio.png')     
plt.close()

#Box plot SUVmean of FET and TER rois
xtick_labels = ['FET', 'PSMA']
fig, ax = plt.subplots()
box_plot(final, ['FET roi on FETPET SUVmean', 'TER roi on TERPET SUVmean'], 'SUVmean [g mL\u207b\u00b9]', ['magenta','magenta'], ymin=0.0, ymax=6.0, tick_min=0.5, tick_max=1.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmeanFETTERrois.png')     
plt.close()

#Box plot SUVmax of FET and TER rois
fig, ax = plt.subplots()
box_plot(final, ['FET roi on FETPET SUVmax', 'TER roi on TERPET SUVmax'], 'SUVmax [g mL\u207b\u00b9]', ['lightgreen','lightgreen'], ymin=0.0, ymax=12.0, tick_min=0.5, tick_max=1.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmaxFETTERrois.png')     
plt.close()

#Box plot TBRmean MI 
fig, ax = plt.subplots()
box_plot(final, ['FET TBRmean MI', 'TER TBRmean MI'], 'TBRmean MI [a.u.]', ['orange','orange'], tick_min=1.0, tick_max=5.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/TBRmeanMI.png')     
plt.close()

#Box plot TBRmax MI 
fig, ax = plt.subplots()
box_plot(final, ['FET TBRmax MI', 'TER TBRmax MI'], 'TBRmax MI [a.u.]', ['blue','blue'], tick_min=1.0, tick_max=5.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/TBRmaxMI.png')     
plt.close()

#Box plot TBRmean CS 
fig, ax = plt.subplots()
box_plot(final, ['FET TBRmean CS', 'TER TBRmean CS'], 'TBRmean CS [a.u.]', ['orange','orange'], tick_min=5.0, tick_max=10.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/TBRmeanCS.png')     
plt.close()

#Box plot TBRmax CS 
fig, ax = plt.subplots()
box_plot(final, ['FET TBRmax CS', 'TER TBRmax CS'], 'TBRmax CS [a.u.]', ['blue','blue'], tick_min=1.0, tick_max=5.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/TBRmaxCS.png')     
plt.close()

#Box plot PSMA SUVmean PsP vs RT MI
xtick_labels = ['PsP', 'RT']
fig, ax = plt.subplots()
box_plot(final, ['PsP PSMA SUVmean MI', 'RT PSMA SUVmean MI'], 'PSMA SUVmean [g mL\u207b\u00b9]',['yellow','yellow'], tick_min=0.1, tick_max=0.5)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmean_PsPvsRT_MI.png') 
plt.close()

#Box plot PSMA SUVmax PsP vs RT MI
xtick_labels = ['PsP', 'RT']
fig, ax = plt.subplots()
box_plot(final, ['PsP PSMA SUVmax MI', 'RT PSMA SUVmax MI'], 'PSMA SUVmax [g mL\u207b\u00b9]',['green','green'], tick_min=0.5, tick_max=1.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmax_PsPvsRT_MI.png') 
plt.close()

#Box plot PSMA SUVmean PsP vs RT CS
xtick_labels = ['PsP', 'RT']
fig, ax = plt.subplots()
box_plot(final, ['PsP PSMA SUVmean CS', 'RT PSMA SUVmean CS'], 'PSMA SUVmean [g mL\u207b\u00b9]',['yellow','yellow'], tick_min=0.1, tick_max=0.5)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmean_PsPvsRT_CS.png') 
plt.close()

#Box plot PSMA SUVmax PsP vs RT CS
xtick_labels = ['PsP', 'RT']
fig, ax = plt.subplots()
box_plot(final, ['PsP PSMA SUVmax CS', 'RT PSMA SUVmax CS'], 'PSMA SUVmax [g mL\u207b\u00b9]',['green','green'], tick_min=0.5, tick_max=1.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmax_PsPvsRT_CS.png') 
plt.close()


#Box plot PSMA SUVmean CE recurrent tum vs CE CR changes
xtick_labels = ['CE recurrent tumour', 'CE CR changes']
fig, ax = plt.subplots()
box_plot(final, ['FETmaxCE on TERPET SUVmean', 'FETminCE on TERPET SUVmean'], 'PSMA SUVmean [g mL\u207b\u00b9]', ['red','red'], tick_min=0.1, tick_max=0.5)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmean_CEtest.png')     
plt.close()

#Box plot PSMA SUVmax CE recurrent tum vs CE CR changes
xtick_labels = ['CE recurrent tumour', 'CE CR changes']
fig, ax = plt.subplots()
box_plot(final, ['FETmaxCE on TERPET SUVmax', 'FETminCE on TERPET SUVmax'], 'PSMA SUVmax [g mL\u207b\u00b9]', ['violet','violet'], tick_min=0.5, tick_max=1.0)
ax.set_xticklabels(xtick_labels, rotation=0, fontsize=15)
plt.savefig(plot_dir + '/SUVmax_CEtest.png')     
plt.close()



#TBRmean MI vs CS comparison  
FET = [final['FET TBRmean MI'].tolist(), final['FET TBRmean CS'].tolist()]
PSMA = [final['TER TBRmean MI'].tolist(), final['TER TBRmean CS'].tolist()]

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True)
#ax.set_ylabel('TBRmean [a.u.]', fontsize=12, wrap=True)
# first boxplot pair
bp = ax2.boxplot(FET, positions = [1, 2], widths = 0.6)
setBoxColors(bp)
# second boxplot pair
bp1 = ax.boxplot(PSMA, positions = [4, 5], widths = 0.6)
setBoxColors(bp1)

# zoom-in / limit the view to different portions of the data
ax2.set_ylim(1.0, 4.0)  # outliers only
ax.set_ylim(5, 85)  # most of the data


# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

ax2.set_ylabel('TBRmean [a.u.]', fontsize=12)
ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)

## set axes limits and labels
xlim(0,6)
#ylim(0,80)
ax2.set_xticklabels(['FET', 'PSMA'], fontsize=12)
ax2.set_xticks([1.5, 4.5])
#ax.yaxis.set_major_locator(plt.MultipleLocator(10.0))
#ax.yaxis.set_minor_locator(plt.MultipleLocator(5.0))

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # bottom-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-right diagonal


kwargs.update(transform=ax2.transAxes)  # switch to the top axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # top-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # top-rig


# draw temporary red and blue lines and use them to create a legend
hB, = plot([1,1],'b-')
hR, = plot([1,1],'r-')
legend((hB, hR),('MI', 'CS'), loc='lower right')
hB.set_visible(False)
hR.set_visible(False)

plt.savefig(plot_dir + '/TBRmeanMIvsCS.png')     
plt.close()


#TBRmax MI vs CS comparison  
FET = [final['FET TBRmax MI'].tolist(), final['FET TBRmax CS'].tolist()]
PSMA = [final['TER TBRmax MI'].tolist(), final['TER TBRmax CS'].tolist()]

fig, (ax, ax2) = plt.subplots(2, 1, sharex=True)
#ax.set_ylabel('TBRmean [a.u.]', fontsize=12, wrap=True)
# first boxplot pair
bp = ax2.boxplot(FET, positions = [1, 2], widths = 0.6)
setBoxColors(bp)
# second boxplot pair
bp1 = ax2.boxplot(PSMA, positions = [4, 5], widths = 0.6)
setBoxColors(bp1)

bp1 = ax.boxplot(PSMA, positions = [4, 5], widths = 0.6)
setBoxColors(bp1)

# zoom-in / limit the view to different portions of the data
ax2.set_ylim(1.0, 5.0)  # outliers only
ax.set_ylim(5, 50)  # most of the data


# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

ax2.set_ylabel('TBRmax [a.u.]', fontsize=12)
ax2.yaxis.set_label_coords(0.05, 0.5, transform=fig.transFigure)

## set axes limits and labels
xlim(0,6)
#ylim(0,80)
ax2.set_xticklabels(['FET', 'PSMA'], fontsize=12)
ax2.set_xticks([1.5, 4.5])
#ax.yaxis.set_major_locator(plt.MultipleLocator(10.0))
#ax.yaxis.set_minor_locator(plt.MultipleLocator(5.0))

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-d, +d), (-d, +d), **kwargs)        # bottom-left diagonal
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-right diagonal


kwargs.update(transform=ax2.transAxes)  # switch to the top axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # top-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # top-rig


# draw temporary red and blue lines and use them to create a legend
hB, = plot([1,1],'b-')
hR, = plot([1,1],'r-')
legend((hB, hR),('MI', 'CS'), loc='lower right')
hB.set_visible(False)
hR.set_visible(False)

plt.savefig(plot_dir + '/TBRmaxMIvsCS.png')     
plt.close()


