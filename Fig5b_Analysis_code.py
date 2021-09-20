# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:43:50 2021

@author: Jewett Lab
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#plotting settings
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
plt.rcParams["figure.figsize"] = [5,2]

#copy data from N1-129/210524_combined_data.xlsx
data = pd.read_csv('Fig5b_Normalized_Data.csv', index_col = 'Ni (μM)')

#plotting
sns.heatmap(data, vmin = 0, vmax = 0.35, cmap = 'mako_r', square = True, linewidths = 0.5,cbar_kws={"shrink": 0.5})
#plt.savefig('210812_n1-129_heatmap.pdf', bbox_inches = 'tight')