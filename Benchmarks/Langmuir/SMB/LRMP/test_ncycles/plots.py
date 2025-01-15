# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 08:24:33 2024

@author: jespfra
"""


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd


data_CADETDG = pd.read_csv('CADETDGConvergence.csv')
data_CADETJulia = pd.read_csv('CADETJuliaConvergence.csv')


# Makers
markers = ['*', '^','h']
markerss = []

# Creating custom legend handles
patch1 = mpatches.Patch(color='purple', label='DG-Julia')
patch3 = mpatches.Patch(color='orangered', label='DG-C++')


# Get the unique settings from the 'settings' column
unique_settings = data_CADETDG['Setting'].unique()

i = 0
plt.figure()
for setting in unique_settings:
    # Find the index for where the settings match the current unique setting
    idxx = data_CADETDG['Setting'] == setting
    
    plt.plot(data_CADETDG['nCycles'][idxx],data_CADETDG['rtime'][idxx], '--', marker = markers[i], color='orangered',label= f"DG-C++, {setting}",markersize=12)
    plt.plot(data_CADETJulia['nCycles'][idxx],data_CADETJulia['rtime'][idxx], '--', marker = markers[i], color='purple',label=f'DG-Julia, {setting}',markersize=12)
    
    markerss.append(plt.Line2D([], [], color='black', marker=f'{markers[i]}', linestyle='None', label=f"$N_d^b$ = {data_CADETDG['polyDeg'][idxx].iloc[0]}, $N_e$ = {data_CADETDG['nCells'][idxx].iloc[0]}",markersize=12))

    i += 1

handles = [patch1, patch3] + markerss[:]
plt.xlabel('Number of cycles', fontsize=12)
plt.ylabel('Compute time / s', fontsize=12)
plt.tick_params(axis='both', which='major', labelsize=12)
plt.legend(ncol=1, handles=handles, fontsize=11)
plt.savefig('plot_runtime_cycles.svg',format = 'svg',dpi = 1200)

