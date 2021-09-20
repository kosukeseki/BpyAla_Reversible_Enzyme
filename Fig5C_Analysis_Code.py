# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 13:27:16 2021

@author: Kosuke
"""
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np

#fitting to michaelis-menten with substrate inhibition
def func(x,vm,km,ki):
    return vm * x / (km + x * (1 + x/ki))

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
plt.rcParams["figure.figsize"] = [5.5,2]
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)

#make square axis for nice graphs
def make_square_axis(ax):
    ax.set_aspect(1 / ax.get_data_ratio())

metal_df = pd.read_csv('Fig5c_Metal_Data.csv')
nometal_df = pd.read_csv('Fig5c_No_Metal_Data.csv')

#calculated numbers from michaelis-menten kinetics experiment
x_metal = metal_df['Luciferin (uM)']
y_metal = metal_df['156_rfu']
yerror_metal = metal_df['156_stdev']

x_nometal = nometal_df['Luciferin (uM)']
y_nometal = nometal_df['156_rfu']
yerror_nometal = nometal_df['156_stdev']

#intiial guess
p0_m = [10000,5,600]
p0 = [10000,5,600]

#curve fitting
metal_popt, pcov = curve_fit(func,x_metal,y_metal,p0_m)
popt,pcov = curve_fit(func,x_nometal,y_nometal,p0)

#getting fits
xfit = np.linspace(0,750,1000)
yfit_metal = func(xfit,*metal_popt)
yfit_nometal = func(xfit,*popt)

#plotting
fig = plt.figure()
ax = fig.add_subplot(111)

#plotting the fit
plt.plot(xfit,yfit_metal, color = 'black')
plt.plot(xfit,yfit_nometal, color = 'black')

#plotting true nubmers
plt.errorbar(x_nometal,y_nometal,yerr = yerror_nometal, fmt = 'o',color = 'gray', ecolor = 'black', capsize = 3, label = '- Ni\n$K_m$ = {:.3f}\n$k_{{cat}}$ = {:.0f}'.format(popt[1], popt[0]))
plt.errorbar(x_metal,y_metal,yerr = yerror_metal, fmt = 'o',color = 'green', ecolor = 'black', capsize = 3, label = '+ Ni\n$K_m$ = {:.3f}\n$k_{{cat}}$ = {:.0f}'.format(metal_popt[1], metal_popt[0]))

#formatting graph
plt.xlabel('Luciferin (μM)')
plt.ylabel('RLU')
make_square_axis(ax)
plt.title('202/532')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.legend(loc=(1.04,0), frameon = False)
#plt.savefig('210603_n1-130-2_156_combined.pdf', bbox_inches = 'tight')