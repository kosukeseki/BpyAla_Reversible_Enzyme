# BpyAla_Reversible_Enzyme
Data analysis codes for analyzing enzyme kinetic data in 384-well plates. 

This code was run on Widows 10 Home Edition, 64-bit operating system. Code was written on Spyder 3.3.2 running IPython 7.2.0 using Python 3.7.1. Code was also tested on Jupyter Notebook by a seperate user and was able to replicate the results. Python was installed using Anaconda. 

To run data, place all files included in this repository to an arbitrary directory. Open any of the .py files on Spyder 3.3.2 and run script. 

For figures made in Figure 3, the scripts outputs a bar graph of enzyme rates of mutants tested in this manuscript. Output slope data is outputted to the "slope_dict" dictionary, which contains each enzyme variant and condition as keys and the average slope, replicates of slope, and the standard deviation as its values in tuple form. Each enzyme variant is labelled in numbers between 25-51, which can be linked to the mutations made by cross-referencing the "species" and "mutations" variables. This script generally takes about a minute to run. 

For figures made in Figure 5, the 5b scripts output a heatmap of Pluc luciferase activity in response to a variety of Nickel (II) concentrations. The data file included in this set is represents the average RLU value from the experiments, whose raw data can be found in the supplemental data file. Raw data from this Michaelis-Menten analysis is also availabe in the supplemental data file. 156 repreents the 202/532 variant, while 186 represents the 108/508 variant.
