# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 17:15:56 2021

@author: Kosuke
"""

#importing libraries for data analysis
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import linregress

###formatting for plotting
plt.rcParams['font.sans-serif']="Arial"
plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['ytick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams["figure.figsize"] = [4.25,4.25]

#importing data and tidying it
rfudf = pd.read_csv('Fig3a_Formatted_Raw_Data.csv')
##col_info is a csv that specifies the reaction in each well
##i.e. A1 - this enzyme, with this concentration of metal, this rep#
col_info = pd.read_csv('Fig3a_Well_Info.csv')
#Arranging into tidydata format
melt = rfudf.melt(id_vars = 'Time', 
                  var_name = 'Well',
                  value_name = 'rfu')
full_data = pd.merge(melt, 
                col_info, 
                how = 'left', 
                on = 'Well')

#defining variables and conditions
conditions = [100, 10, 1, 0]
species = ['25', '27', '32', '35',
            '40', '41', '42', '49',
            '50', '26', '28', '29',
            '30', '31', '33', '34',
            '36', '37', '38', '39',
            '43', '44', '45', '46',
            '47', '48', '51', 'WT',
            'GFP', '-']
mutations = ['Blank','216/516','214/442','159/513','167/517','158/512',
             '169/516','237/439','238/516','156/513','161/516','215/520',
             '167/513','169/512','156/516','238/442','170/513','158/514',
             '167/516','191/513','159/512','191/517','215/524','169/510',
             '159/517','191/520','169/514','214/524','2TAG-sfGFP', 'WT']
time = full_data.iloc[0:41,0]
replicate = [1,2,3]
species = sorted(species)
conditions = sorted(conditions)
##setting up sliding window length of five points
window = 5

#setting up plots and subplots
fig, ax = plt.subplots(nrows = 6,
                       ncols = 5, 
                       sharex = True, 
                       sharey = True)
colors = plt.cm.Greens(np.linspace(0.3,0.9,len(conditions)))

##setting up dictionaries to store values
dict = {}
slope_rep_dict = {}
slope_dict = {}

#setting up loops to calculate rates
for enz in species:
    for conc in conditions:
        
        ###only look at data of one enz @ one ni conc
        mini_data = full_data.loc[(full_data['species'] == enz)&
            (full_data['Ni (uM)'] == conc)]
        
        ###initialize an empty slope variable for all slopes 
        ### over a window of time
        allslope = []
        allslopestd = []
        slope_reps = []
        
        ###set up sliding window
        for t in range(0,len(time[0:-(window-1)])):
            timepoints = time[t:t+window]
            timepoints = timepoints.tolist()
            
            temp_slope = []
            temp_rsquared = []
            
            for rep in replicate:

            ##initialize empty array for rfu values
                rfu = []

            ###for each point in the sliding window, find all rfu
            ###that correspondong to that time point
                for x in timepoints:
                    rfu_reps = mini_data['rfu'].loc[(mini_data['Time'] == x)&(mini_data['replicate'] == rep)]
                    rfu_reps = rfu_reps.tolist()
                    rfu.append(rfu_reps)
                
                ## change rfu from list of tuples to just a list
                rfu = [y[0] for y in rfu]
                    
                ###after populating rfu for a window of time
                ###(timepoints), do regression on those points
                var = linregress(timepoints, rfu)
                slope = var[0]
                rsquared = var[2]*var[2]
            
                ###append those values to temp_slope,temp_rsquared
                temp_slope.append(slope)
                temp_rsquared.append(rsquared)

            ###keep replicates for plotting
            slope_reps.append(temp_slope)
            
            ###after the loop over replicates, we have replicates of slopes 
            ###for each window of time. now avg them together. 
            avgslope = np.mean(temp_slope)
            stdevslope = np.std(temp_slope)
            
            ###store the avg slope in allslope, which captures slope over time
            allslope.append(avgslope)
            allslopestd.append(stdevslope)
            
        ##once all slopes over all time have been calculated, append that 
        ##to dictionary to store the value
        dict[enz+str(conc)] = (allslope, allslopestd)
        slope_rep_dict[enz+str(conc)] = slope_reps
    
        ###find the max slope and its stdev over time period, keep replicates
        maxslope = max(dict[enz+str(conc)][0])
        maxindex = dict[enz+str(conc)][0].index(maxslope)
        maxslopereps = slope_rep_dict[enz+str(conc)][maxindex]
        maxslopestd = dict[enz+str(conc)][1][maxindex]
    
        ###store maxslope, stdev at maxslope in new dict
        slope_dict[enz+str(conc)] = (maxslope,maxslopereps,maxslopestd, enz, conc)

###plotting bar plots
count1 = 0
count2 = 0
count3 = 0

test = pd.DataFrame.from_dict(slope_dict, orient = 'index')
test2 = test.explode(1)
    
for enz in species:
    
    ###plotting individual replicates
    swarmdata = test2.loc[test2[3] == enz]
    sns.swarmplot(x = swarmdata[4], y = swarmdata[1], size = 2, ax = ax[count1,count2], color = 'black', alpha = 0.5)
    ax[count1,count2].set(xlabel = None, ylabel = None)
    
    ###plotting bar charts for rates
    for conc in conditions:
        ax[count1,count2].bar(str(conc),slope_dict[enz+str(conc)][0], yerr = slope_dict[enz+str(conc)][2],color = colors[conditions.index(conc)],
                              error_kw = {'capsize': 2, 'lw': 0.5, 'capthick': 0.5})
        ax[count1,count2].set_title(mutations[count3], fontsize = 7)
    count2 = count2 + 1
    if count2 > 4:
        count1 = count1 + 1
        count2 = 0
    count3 = count3 + 1

for ax in fig.axes:
    plt.sca(ax)
    plt.xticks(rotation=45)


###adjustments on plots
fig.text(0.5,0.05, 'Ni (μM)', ha='center', va='center') 
fig.text(0.05, 0.5, 'Rate (AU/min) at 30°C', va = 'center', rotation = 'vertical')
fig.subplots_adjust(left = 0.15, bottom = 0.13, wspace = 0.4, hspace = 0.8)
fig.savefig('Figure 3a Graph', bbox_inches='tight')

