# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 14:14:43 2025

@author: pbzit
"""


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pandas as pd
import numpy as np



def plot_convergence(CADETFVdata,CADETJuliadata,CADETDGdata,c_analytical,saveLocation="",model="LRM"):
    colorFV = 'deeppink'
    colorDG_exact = plt.cm.Blues(range(230, 90, -30))   # Nuances of red
    colorDG_inexact =  plt.cm.Reds(range(255, 90, -40))  # Nuances of red
    colorDGJulia_exact = plt.cm.Greens(range(230, 90, -30))   # Nuances of red
    colorDGJulia_inexact =  plt.cm.Purples(range(255, 90, -30))  # Nuances of red
    colors_profile1 = ['g','r','b']
    colors_profile2 = ['m', 'c', 'y']
    markersize = 14


    fig,ax = plt.subplots(figsize=(9, 9)) #figsize=(15, 13)
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['MAE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)

    
    # Plot CADET-DG 
    ax.loglog(CADETDGdata['DOF'][:],CADETDGdata["MAE"][:],'-', label = f'DG-C++' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[0])
    
    
    # Julia plot
    ax.loglog(CADETJuliadata['DOF'][:],CADETJuliadata["MAE"][:],'.--', label = f'DG-Julia' ,markersize=markersize, linewidth=2,color=colorDGJulia_inexact[0])


    ax.set_xlabel('Degrees of freedom', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, framealpha=1.0)
    fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_convergence_'+saveLocation+'.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
    # Plot profiles
    fig,ax = plt.subplots(figsize=(9, 7))
    ax.plot(c_analytical['C1'], linewidth=3)
    ax.set_xlabel('Time / s', fontsize=25)
    ax.set_ylabel('Concentration / mM', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=22)
    fig.subplots_adjust(bottom=0.1)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_profiles_'+saveLocation+'.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
#%%
# Plotting the acyclic case 
CADETFVdata = pd.read_csv('acyclic/CADETFVConvergence.csv')
CADETJuliadata = pd.read_csv('acyclic/CADETJuliaConvergence.csv')
CADETDGdata = pd.read_csv('acyclic/CADETDGConvergence.csv')
c_analytical = pd.read_csv('acyclic/Semi-analytical_acyclic.csv')
plot_convergence(CADETFVdata,CADETJuliadata[:],CADETDGdata[:], c_analytical, 'acyclic')



# Plotting the cyclic case 
CADETFVdata = pd.read_csv('cyclic/CADETFVConvergence.csv')
CADETJuliadata = pd.read_csv('cyclic/CADETJuliaConvergence.csv')
CADETDGdata = pd.read_csv('cyclic/CADETDGConvergence.csv')
c_analytical = pd.read_csv('cyclic/Semi-analytical_cyclic.csv')
plot_convergence(CADETFVdata,CADETJuliadata[:],CADETDGdata[:], c_analytical, 'cyclic')