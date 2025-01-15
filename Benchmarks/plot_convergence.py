# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 12:07:39 2023

@author: jespfra
A module that contains a plot function, a plot initiator and all the function calls to generate plots to the batch benchmarks

"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import pandas as pd
import numpy as np


def plot_convergence(CADETFVdata,CADETJuliadata,CADETDGdata=[],saveLocation="",model="LRM"):
    
    colorFV = 'deeppink'
    colorDG_exact = plt.cm.Blues(range(230, 90, -30))   # Nuances of red
    colorDG_inexact =  plt.cm.Reds(range(255, 90, -40))  # Nuances of red
    colorDGJulia_exact = plt.cm.Greens(range(230, 90, -30))   # Nuances of red
    colorDGJulia_inexact =  plt.cm.Purples(range(255, 90, -30))  # Nuances of red
    colors_profile1 = ['g','r','b']
    colors_profile2 = ['m', 'c', 'y']
    markersize = 14

    if model == "GRM":
        tag = "polyDegPoreu"
        plottag = "$N_d^p$"
    else:
        tag = "polyDegu"
        plottag = "$N_d^b$"
        
        
    
    fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['maxE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_e"],'.--', label = 'DG-Julia, Exact')
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_i"],'.--', label = 'DG-Julia, Collocation')
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'DG-Julia, Collocation, {plottag}={CADETJuliadata[tag][idx].min()}',markersize=markersize, linewidth=2,color=colorDGJulia_inexact[i] )

    
    # Plot CADET-DG 
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_e"][idx],'-', label = f'DG-C++, Exact, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_exact[i])
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_i"][idx],'-', label = f'DG-C++, Collocation, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[i])


    ax.set_xlabel('Degrees of freedom', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_convergence.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
    # DOF Rutime plot
    fig,ax = plt.subplots(figsize=(11.5, 12.5))
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['runtime'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_e"],'.--', label = 'DG-Julia, Exact')
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_i"],'.--', label = 'DG-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'DG-Julia, Collocation, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_inexact[i])
        

        
    # Plot DG-C++
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_e"][idx],'-', label = f'DG-C++, Exact, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_exact[i])
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_i"][idx],'-', label = f'DG-C++, Collocation, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[i])
            
    ax.set_xlabel('Degrees of freedom', fontsize=25)
    ax.set_ylabel('Compute time / s', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.45)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_runtime.svg'),format = 'svg',dpi = 1200)
    
    
    # Runtime Error plot
    fig,ax = plt.subplots(figsize=(11.5, 10))
    ax.loglog(CADETFVdata['runtime'],CADETFVdata['maxE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    # plt.loglog(CADETJuliadata["runtime_e"],CADETJuliadata["maxError_e"],'.--', label = 'DG-Julia, Exact')
    # plt.loglog(CADETJuliadata["runtime_i"],CADETJuliadata["maxError_i"],'.--', label = 'DG-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'DG-Julia, Collocation, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_inexact[i])

        
    # Plot CADET-DG 
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            ax.loglog(CADETDGdata['runtime_e'][idx],CADETDGdata["maxE_e"][idx],'-', label = f'DG-C++, Exact, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_exact[i])
            ax.loglog(CADETDGdata['runtime_i'][idx],CADETDGdata["maxE_i"][idx],'-', label = f'DG-C++, Collocation, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[i])
            
   
    ax.set_xlabel('Compute time / s', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.4)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_err_runtime.svg'),format = 'svg',dpi = 1200)
    
    
    
    
    
    
    ############ Plotting only Collocation ############
    fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['maxE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_e"],'.--', label = 'DG-Julia, Exact')
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_i"],'.--', label = 'DG-Julia, Collocation')
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        # ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'DG-Julia, {plottag}={CADETJuliadata[tag][idx].min()}',markersize=markersize+4, linewidth=2,color=colorDGJulia_inexact[i] )
        # ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_fbdf"][idx],'.--', label = f'DG-Julia, FBDF, {plottag}={CADETJuliadata[tag][idx].min()}',markersize=markersize, linewidth=2,color=colorDGJulia_exact[i] )

    
    # Plot CADET-DG 
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            # ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_e"][idx],'-', label = f'DG-C++, Exact, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_exact[i])
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["maxE_i"][idx],'-', label = f'DG-C++, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[i])


    ax.set_xlabel('Degrees of freedom', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.38)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_convergence_i.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
    # DOF Rutime plot
    fig,ax = plt.subplots(figsize=(11.5, 10)) #
    ax.loglog(CADETFVdata['DOF'],CADETFVdata['runtime'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_e"],'.--', label = 'DG-Julia, Exact')
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_i"],'.--', label = 'DG-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        # ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'DG-Julia, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize+4, linewidth=2,color=colorDGJulia_inexact[i])
        # ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_fbdf"][idx],'.--', label = f'DG-Julia, FBDF, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        

        
    # Plot CADET-DG
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            # ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_e"][idx],'-', label = f'DG-C++, Exact, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_exact[i])
            ax.loglog(CADETDGdata['DOF'][idx],CADETDGdata["runtime_i"][idx],'-', label = f'DG-C++, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[i])
            
    ax.set_xlabel('Degrees of freedom', fontsize=25)
    ax.set_ylabel('Compute time / s', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.38)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_runtime_i.svg'),format = 'svg',dpi = 1200)
    
    
    # Runtime Error plot
    fig,ax = plt.subplots(figsize=(11.5, 11)) #
    ax.loglog(CADETFVdata['runtime'],CADETFVdata['maxE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    # plt.loglog(CADETJuliadata["runtime_e"],CADETJuliadata["maxError_e"],'.--', label = 'DG-Julia, Exact')
    # plt.loglog(CADETJuliadata["runtime_i"],CADETJuliadata["maxError_i"],'.--', label = 'DG-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        # ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'DG-Julia, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize+4, linewidth=2,color=colorDGJulia_inexact[i])
        # ax.loglog(CADETJuliadata['runtime_fbdf'][idx],CADETJuliadata["maxE_fbdf"][idx],'.--', label = f'DG-Julia, FBDF, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])

        
    # Plot CADET-DG 
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        for i in range(CADETDGdata[tag].nunique()):
            idx = slice(i * CADETDGdata["nCellu"].nunique(), (i + 1) * CADETDGdata["nCellu"].nunique())
            # ax.loglog(CADETDGdata['runtime_e'][idx],CADETDGdata["maxE_e"][idx],'-', label = f'DG-C++, Exact, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_exact[i])
            ax.loglog(CADETDGdata['runtime_i'][idx],CADETDGdata["maxE_i"][idx],'-', label = f'DG-C++, {plottag}={CADETDGdata[tag][idx].min()}' ,markersize=markersize, marker = '*', linewidth=2,color=colorDG_inexact[i])
            
   
    ax.set_xlabel('Compute time / s', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.38)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_err_runtime_i.svg'),format = 'svg',dpi = 1200)
    
    
    
    
    ############ Plotting only Julia exact vs. collocation ############
    fig,ax = plt.subplots(figsize=(11.5, 10)) #figsize=(15, 13)
    # ax.loglog(CADETFVdata['DOF'],CADETFVdata['maxE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_e"],'.--', label = 'DG-Julia, Exact')
    # plt.loglog(CADETJuliadata['DOF'],CADETJuliadata["maxError_i"],'.--', label = 'DG-Julia, Collocation')
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'DG-Julia, Collocation, {plottag}={CADETJuliadata[tag][idx].min()}',markersize=markersize, linewidth=2,color=colorDGJulia_inexact[i] )



    ax.set_xlabel('Degrees of freedom', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend()
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.38)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_convergence_julia.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
    # DOF Rutime plot
    fig,ax = plt.subplots(figsize=(11.5, 10))
    # ax.loglog(CADETFVdata['DOF'],CADETFVdata['runtime'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_e"],'.--', label = 'DG-Julia, Exact')
    # plt.plot(CADETJuliadata['DOF'],CADETJuliadata["runtime_i"],'.--', label = 'DG-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['DOF'][idx],CADETJuliadata["runtime_i"][idx],'.--', label = f'DG-Julia, Collocation, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_inexact[i])
        
        
       
    ax.set_xlabel('Degrees of freedom', fontsize=25)
    ax.set_ylabel('Compute time / s', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.38)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_runtime_julia.svg'),format = 'svg',dpi = 1200)
    
    
    # Runtime Error plot
    fig,ax = plt.subplots(figsize=(11.5, 11))
    # ax.loglog(CADETFVdata['runtime'],CADETFVdata['maxE'],':', label = 'FV-C++',markersize=markersize, marker = '^', linewidth=2, color = colorFV)
    # plt.loglog(CADETJuliadata["runtime_e"],CADETJuliadata["maxError_e"],'.--', label = 'DG-Julia, Exact')
    # plt.loglog(CADETJuliadata["runtime_i"],CADETJuliadata["maxError_i"],'.--', label = 'DG-Julia, Collocation')
    
    for i in range(CADETJuliadata[tag].nunique()):
        idx = slice(i * CADETJuliadata["nCellu"].nunique(), (i + 1) * CADETJuliadata["nCellu"].nunique())
        ax.loglog(CADETJuliadata['runtime_e'][idx],CADETJuliadata["maxE_e"][idx],'.--', label = f'DG-Julia, Exact, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_exact[i])
        ax.loglog(CADETJuliadata['runtime_i'][idx],CADETJuliadata["maxE_i"][idx],'.--', label = f'DG-Julia, Collocation, {plottag}={CADETJuliadata[tag][idx].min()}' ,markersize=markersize, linewidth=2,color=colorDGJulia_inexact[i])

        

   
    ax.set_xlabel('Compute time / s', fontsize=25)
    # ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.set_ylabel('Max abs error / mM', fontsize=25)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.tick_params(axis='both', which='major', labelsize=22)
    # plt.title('LRM Langmuir')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
    ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
    fig.subplots_adjust(bottom=0.38)  # Adjust this value as needed
    plt.savefig(os.path.join(saveLocation,'Plot_err_runtime_julia.svg'),format = 'svg',dpi = 1200)
    
    
    
    
    
    
    ############ CSS evaluations ############
    if os.path.exists(path + "CSS/CADETJuliaConvergenceCSS.csv"):
        
        # Split the data into coupled and OCA-FPI data 
        CSSData_coupled = pd.read_csv(path + "CSS/CADETJuliaConvergenceCSS.csv", 
                                      usecols=[
                                          'DOF',
                                          'nCellu',
                                          'polyDegu',
                                          'polyDegPoreu',
                                          'maxE_coupled_i',
                                          'L2_coupled_i',
                                          'L2_coupled_analytical_i',
                                          'runtime_iter_coupled_i',
                                          ])
        
        CSSData_oca_fpi = pd.read_csv(path + "CSS/CADETJuliaConvergenceCSS.csv", 
                                      usecols=[
                                          'DOF',
                                          'nCellu',
                                          'polyDegu',
                                          'polyDegPoreu',
                                          'maxE_oca_fpi_i',
                                          'L2_oca_fpi_i',
                                          'L2_oca_fpi_analytical_i',
                                          'runtime_iter_oca_fpi_i',
                                          ])
        # Remove zeros
        CSSData_coupled = CSSData_coupled[CSSData_coupled['L2_coupled_i'] != 0].reset_index(drop=True)
        CSSData_oca_fpi = CSSData_oca_fpi[CSSData_oca_fpi['L2_oca_fpi_i'] != 0].reset_index(drop=True)
        
        
        
        if os.path.exists(path + "CSS/CADETDGConvergenceCSS.csv"):
            CSSData_coupled_DG = pd.read_csv(path + "CSS/CADETDGConvergenceCSS.csv",
                                             usecols=[
                                            'DOF',
                                            'nCellu',
                                            'polyDegu',
                                            'polyDegPoreu',
                                            'maxE_coupled_i',
                                            'L2_coupled_i',
                                            'L2_coupled_analytical_i',
                                            'runtime_iter_coupled_i',
                                            'runtime_iter_coupled_i_tot',])
            CSSData_coupled_DG = CSSData_coupled_DG[CSSData_coupled_DG['L2_coupled_i'] != 0].reset_index(drop=True)
            
            
            if 'maxE_oca_fpi_i' in pd.read_csv(path + "CSS/CADETDGConvergenceCSS.csv"):
                CSSData_oca_fpi_DG = pd.read_csv(path + "CSS/CADETDGConvergenceCSS.csv", 
                                              usecols=[
                                                  'DOF',
                                                  'nCellu',
                                                  'polyDegu',
                                                  'polyDegPoreu',
                                                  'maxE_oca_fpi_i',
                                                  'L2_oca_fpi_i',
                                                  'L2_oca_fpi_analytical_i',
                                                  'runtime_iter_oca_fpi_i',
                                                  'runtime_iter_oca_fpi_i_tot',
                                                  ])
                CSSData_oca_fpi_DG = CSSData_oca_fpi_DG[CSSData_oca_fpi_DG['L2_oca_fpi_i'] != 0].reset_index(drop=True)
            
        
        if os.path.exists(path + "CSS/CADETFVConvergenceCSS.csv"):
            CSSData_coupled_FV = pd.read_csv(path + "CSS/CADETFVConvergenceCSS.csv",
                                             usecols=[
                                            'DOF',
                                            'nCellu',
                                            'polyDegu',
                                            'polyDegPoreu',
                                            'maxE_coupled_i',
                                            'L2_coupled_i',
                                            'L2_coupled_analytical_i',
                                            'runtime_iter_coupled_i',
                                            'runtime_iter_coupled_i_tot',])
            CSSData_coupled_FV = CSSData_coupled_FV[CSSData_coupled_FV['L2_coupled_i'] != 0].reset_index(drop=True)        
        
        
        # Here we want to plot simulation time vs. L2, MaxE, DOF vs. RT
        lines = ['--', '-', '--',':']
        markers = ['.','*', '^','h']
        colors3 = ['tab:orange', 'tab:red', 'tab:cyan', 'tab:gray', 'tab:blue', 
                   'tab:green', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:olive', 
                   'gold', 'lightblue', 'lightgreen', 'lightgray', 
                   'lightpink', 'lightcoral', 'lime', 'darkblue', 
                   'darkgreen', 'darkred']
        
        nCells = CSSData_coupled['nCellu'].unique()
        nPoly = CSSData_coupled[tag].unique()
        
        # number of iterations vs. L2 
        # Unique cells, unique polynomials, eval indices and then plot 
        fig,ax = plt.subplots(figsize=(9, 9)) #
        j = 0
        for i in range(CSSData_coupled[tag].nunique()):
            for l in range(CSSData_coupled['nCellu'].nunique()):
                
                # Index coupled approach 
                idx_coupled = CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.min()
                idx_coupled = slice(idx_coupled, CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.max(), 1)
                ax.semilogy(range(1, idx_coupled.stop - idx_coupled.start +1), CSSData_coupled['L2_coupled_i'][idx_coupled], lines[0], label = f'DG-Julia, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                
                # Index OCA-FPI 
                idx_oca_fpi = CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.min()
                idx_oca_fpi = slice(idx_oca_fpi, CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.max(), 1)
                ax.semilogy(range(1, idx_oca_fpi.stop - idx_oca_fpi.start +1), CSSData_oca_fpi['L2_oca_fpi_i'][idx_oca_fpi], lines[1], label = f'OCA-FPI, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                j += 1
            
        ax.set_xlabel('Number of iterations', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('$L_2$ Norm', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        # plt.savefig(os.path.join(saveLocation,'CSS/Plot_CSS_L2.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        
        
        # Simulation time vs. L2 
        # Unique cells, unique polynomials, eval indices and then plot 
        j = 0
        fig,ax = plt.subplots(figsize=(9, 9)) #
        for i in range(CSSData_coupled[tag].nunique()):
            for l in range(CSSData_coupled['nCellu'].nunique()):
                
                # Index coupled approach 
                idx_coupled = CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.min()
                idx_coupled = slice(idx_coupled, CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.max(), 1)
                ax.semilogy(CSSData_coupled['runtime_iter_coupled_i'][idx_coupled], CSSData_coupled['L2_coupled_i'][idx_coupled],lines[0], label = f'DG-Julia, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                
                # Index OCA-FPI 
                idx_oca_fpi = CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.min()
                idx_oca_fpi = slice(idx_oca_fpi, CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.max(), 1)
                ax.semilogy(CSSData_oca_fpi['runtime_iter_oca_fpi_i'][idx_oca_fpi], CSSData_oca_fpi['L2_oca_fpi_i'][idx_oca_fpi],lines[1], label = f'OCA-FPI, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                j += 1

            
        ax.set_xlabel('Compute time / s', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('$L_2$ Norm', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.yaxis.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        plt.minorticks_on()
        plt.savefig(os.path.join(saveLocation,'CSS/Plot_CSS_runtime_L2.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        
        # Simulation time vs. maxE
        # Unique cells, unique polynomials, eval indices and then plot 
        fig,ax = plt.subplots(figsize=(9, 9)) #
        j = 0
        for i in range(CSSData_coupled[tag].nunique()):
            for l in range(CSSData_coupled['nCellu'].nunique()):
                
                # Index coupled approach 
                idx_coupled = CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.min()
                idx_coupled = slice(idx_coupled, CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.max(), 1)
                ax.loglog(CSSData_coupled['runtime_iter_coupled_i'][idx_coupled], CSSData_coupled['maxE_coupled_i'][idx_coupled],lines[0], label = f'DG-Julia, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                
                # Index OCA-FPI 
                idx_oca_fpi = CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.min()
                idx_oca_fpi = slice(idx_oca_fpi, CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.max(), 1)
                ax.loglog(CSSData_oca_fpi['runtime_iter_oca_fpi_i'][idx_oca_fpi], CSSData_oca_fpi['maxE_oca_fpi_i'][idx_oca_fpi],lines[1], label = f'OCA-FPI, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                j += 1
            
        ax.set_xlabel('Compute time / s', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'CSS/Plot_CSS_runtime_err.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        
        # iteration vs. Simulation time per iteration
        fig,ax = plt.subplots(figsize=(9, 9)) #
        j = 0
        for i in range(CSSData_coupled[tag].nunique()):
            for l in range(CSSData_coupled['nCellu'].nunique()):
                
                # Index coupled approach 
                idx_coupled = CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.min()
                idx_coupled = slice(idx_coupled, CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.max(), 1)
                ax.plot(range(1, idx_coupled.stop - idx_coupled.start +1), CSSData_coupled['runtime_iter_coupled_i'][idx_coupled], lines[0], label = f'DG-Julia, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                
                # Index OCA-FPI 
                idx_oca_fpi = CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.min()
                idx_oca_fpi = slice(idx_oca_fpi, CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.max(), 1)
                ax.plot(range(1, idx_oca_fpi.stop - idx_oca_fpi.start +1), CSSData_oca_fpi['runtime_iter_oca_fpi_i'][idx_oca_fpi], lines[1], label = f'OCA-FPI, $N_e^b$ = {int(nCells[l])}, {plottag}={int(nPoly[i])}' ,markersize=2, linewidth=2, color = colors3[j % len(colors3)])
                j += 1
            
        ax.set_xlabel('Iteration', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('Cumulated runtime / s', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        ax.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'CSS/Plot_CSS_runtime.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')

        
        
        ### Plots with only the resulting stats 
        DOF = []
        runtime_coupled = []
        maxE_coupled = []
        runtime_oca_fpi = []
        maxE_oca_fpi = []
        markerss = []
        markerss.append(plt.Line2D([], [], color='deeppink', marker='^', linestyle='None', label='FV-C++',markersize=markersize))
        
        
        DOF_DG = []
        runtime_coupled_DG = []
        maxE_coupled_DG = []
        runtime_oca_fpi_DG = []
        maxE_oca_fpi_DG = []
        
        
        DOF_FV = []
        runtime_coupled_FV = []
        maxE_coupled_FV = []
        runtime_oca_fpi_FV = []
        maxE_oca_fpi_FV = []
        
        
        for i in range(CSSData_coupled[tag].nunique()):
            for l in range(CSSData_coupled['nCellu'].nunique()):
                
                # Index coupled approach 
                idx_coupled = CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.min()
                idx_coupled = slice(idx_coupled, CSSData_coupled[(CSSData_coupled['nCellu'] == nCells[l]) & (CSSData_coupled[tag] == nPoly[i])].index.max()+1, 1)
                DOF.append(CSSData_coupled['DOF'][idx_coupled].min())
                maxE_coupled.append(CSSData_coupled['maxE_coupled_i'][idx_coupled].iloc[-1])
                runtime_coupled.append(CSSData_coupled['runtime_iter_coupled_i'][idx_coupled].iloc[-1])
    
                
                # Index OCA-FPI 
                idx_oca_fpi = CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.min()
                idx_oca_fpi = slice(idx_oca_fpi, CSSData_oca_fpi[(CSSData_oca_fpi['nCellu'] == nCells[l]) & (CSSData_oca_fpi[tag] == nPoly[i])].index.max()+1, 1)
                maxE_oca_fpi.append(CSSData_oca_fpi['maxE_oca_fpi_i'][idx_oca_fpi].iloc[-1])
                runtime_oca_fpi.append(CSSData_oca_fpi['runtime_iter_oca_fpi_i'][idx_oca_fpi].iloc[-1])
                
                
            markerss.append(plt.Line2D([], [], color='black', marker=f'{markers[i]}', linestyle='None', label=f'{plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize))
        
        
        if os.path.exists(path + "CSS/CADETDGConvergenceCSS.csv"):
            nPoly_DG = CSSData_coupled_DG[tag].unique()
            nCells_DG = CSSData_coupled_DG['nCellu'].unique()
            for i in range(CSSData_coupled_DG[tag].nunique()):
                for l in range(CSSData_coupled_DG['nCellu'].nunique()):
                    
                    # Index coupled approach 
                    idx_coupled = CSSData_coupled_DG[(CSSData_coupled_DG['nCellu'] == nCells_DG[l]) & (CSSData_coupled_DG[tag] == nPoly_DG[i])].index.min()
                    idx_coupled = slice(idx_coupled, CSSData_coupled_DG[(CSSData_coupled_DG['nCellu'] == nCells_DG[l]) & (CSSData_coupled_DG[tag] == nPoly_DG[i])].index.max()+1, 1)
                    DOF_DG.append(CSSData_coupled_DG['DOF'][idx_coupled].min())
                    maxE_coupled_DG.append(CSSData_coupled_DG['maxE_coupled_i'][idx_coupled].iloc[-1])
                    runtime_coupled_DG.append(CSSData_coupled_DG['runtime_iter_coupled_i'][idx_coupled].max())
                    
                    
                    # Index OCA-FPI 
                    try:
                        idx_oca_fpi = CSSData_oca_fpi_DG[(CSSData_oca_fpi_DG['nCellu'] == nCells_DG[l]) & (CSSData_oca_fpi_DG[tag] == nPoly_DG[i])].index.min()
                        idx_oca_fpi = slice(idx_oca_fpi, CSSData_oca_fpi_DG[(CSSData_oca_fpi_DG['nCellu'] == nCells_DG[l]) & (CSSData_oca_fpi_DG[tag] == nPoly_DG[i])].index.max()+1, 1)
                        maxE_oca_fpi_DG.append(CSSData_oca_fpi_DG['maxE_oca_fpi_i'][idx_oca_fpi].iloc[-1])
                        runtime_oca_fpi_DG.append(CSSData_oca_fpi_DG['runtime_iter_oca_fpi_i'][idx_oca_fpi].max())
                    except:
                        0
        
            comparison = pd.DataFrame({'SpeedUp_C_Jl' : [i / j for i, j in zip(runtime_oca_fpi, runtime_coupled)],
                                       'SpeedUp_C_DG' : [i / j for i, j in zip(runtime_coupled_DG, runtime_coupled)],
                                      'SpeedUp_C_Jl_avg' : np.mean([i / j for i, j in zip(runtime_oca_fpi, runtime_coupled)])*np.ones(len(runtime_oca_fpi)),
                                      'SpeedUp_C_DG_avg' : np.mean([i / j for i, j in zip(runtime_coupled_DG, runtime_coupled)])*np.ones(len(runtime_oca_fpi))})  
            comparison.to_csv(path + 'CSS/comparison.csv')
        
        
        if os.path.exists(path + "CSS/CADETFVConvergenceCSS.csv"):
            nCells_FV = CSSData_coupled_FV['nCellu'].unique()
            for l in range(CSSData_coupled_FV['nCellu'].nunique()):
                
                # Index coupled approach 
                idx_coupled = CSSData_coupled_FV[(CSSData_coupled_FV['nCellu'] == nCells_FV[l])].index.min()
                idx_coupled = slice(idx_coupled, CSSData_coupled_FV[(CSSData_coupled_FV['nCellu'] == nCells_FV[l])].index.max()+1, 1)
                DOF_FV.append(CSSData_coupled_FV['DOF'][idx_coupled].min())
                maxE_coupled_FV.append(CSSData_coupled_FV['maxE_coupled_i'][idx_coupled].iloc[-1])
                runtime_coupled_FV.append(CSSData_coupled_FV['runtime_iter_coupled_i'][idx_coupled].max())
        
        
        # Creating custom legend handles
        patch1 = mpatches.Patch(color='purple', label='DG-Julia')
        patch2 = mpatches.Patch(color='green', label='OCA-FPI')
        patch3 = mpatches.Patch(color='orangered', label='DG-C++')
    
        
        handles = [markerss[0]] + [patch1, patch2, patch3] + markerss[1:]
        
        
        # DOF vs. Simulation time
        fig,ax = plt.subplots(figsize=(9, 9))
        for i in range(CSSData_coupled[tag].nunique()):
            ran = slice(i*CSSData_coupled['nCellu'].nunique(), (i+1)*CSSData_coupled['nCellu'].nunique())
            ax.loglog(DOF[ran],runtime_coupled[ran],markers[i],color="purple", label = f'DG-Julia, {plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize)
            ax.loglog(DOF[ran],runtime_oca_fpi[ran],markers[i],color="green", label = f'OCA-FPI, {plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize)
            
            
        if os.path.exists(path + "CSS/CADETDGConvergenceCSS.csv"):
            for i in range(CSSData_coupled_DG[tag].nunique()):
                ran = slice(i*CSSData_coupled_DG['nCellu'].nunique(), (i+1)*CSSData_coupled_DG['nCellu'].nunique())
                ax.loglog(DOF_DG[ran],runtime_coupled_DG[ran],markers[i],color="orangered", label = f'DG-C++, {plottag} = {int(CSSData_coupled_DG[tag].unique()[i])}',markersize=markersize)
                if len(runtime_oca_fpi_DG)>2:
                    ax.loglog(DOF_DG[ran],runtime_oca_fpi_DG[ran],markers[i],color="red", label = f'OCA-FPI DG, {plottag} = {int(CSSData_coupled_DG[tag].unique()[i])}',markersize=markersize)
        
        
        if os.path.exists(path + "CSS/CADETFVConvergenceCSS.csv"):
            ax.loglog(DOF_FV,runtime_coupled_FV, '^', color="deeppink", label = 'FV-C++',markersize=markersize)
                
          
         
        ax.set_xlabel('Degrees of freedom', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('Compute time / s', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        #Set x-ticks 
        if path == "Linear/SMB/LRM/" or path == "Langmuir/SMB/LRMP/":
            ax.set_xticks([400, 600, 1000, 3000])
            ax.set_xticklabels(['400', '600', '1000', '3000'])
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        plt.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0, handles=handles)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'CSS/Plot_CSS_DOF_runtime.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        
        # DOF vs. MAE
        fig,ax = plt.subplots(figsize=(9, 9))
        for i in range(CSSData_coupled[tag].nunique()):
            ran = slice(i*CSSData_coupled['nCellu'].nunique(), (i+1)*CSSData_coupled['nCellu'].nunique())
            ax.loglog(DOF[ran],maxE_coupled[ran],markers[i],color="purple", label = f'DG-Julia, {plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize)
            ax.loglog(DOF[ran],maxE_oca_fpi[ran],markers[i],color="green", label = f'OCA-FPI, {plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize)
            
            
        if os.path.exists(path + "CSS/CADETDGConvergenceCSS.csv"):
            for i in range(CSSData_coupled_DG[tag].nunique()):
                ran = slice(i*CSSData_coupled_DG['nCellu'].nunique(), (i+1)*CSSData_coupled_DG['nCellu'].nunique())
                ax.loglog(DOF_DG[ran],maxE_coupled_DG[ran],markers[i],color="orangered", label = f'DG-C++, {plottag} = {int(CSSData_coupled_DG[tag].unique()[i])}',markersize=markersize)
                if len(runtime_oca_fpi_DG)>2:
                    ax.loglog(DOF_DG[ran],maxE_oca_fpi_DG[ran],markers[i],color="red", label = f'OCA-FPI, {plottag} = {int(CSSData_coupled_DG[tag].unique()[i])}',markersize=markersize)
        
        
        if os.path.exists(path + "CSS/CADETFVConvergenceCSS.csv"):
            ax.loglog(DOF_FV,maxE_coupled_FV, '^', color="deeppink", label = 'FV-C++',markersize=markersize)
                
          
        ax.set_xlabel('Degrees of freedom', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        plt.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0, handles=handles)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'CSS/Plot_CSS_convergence.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
        
        # MaxE vs. Simulation time
        fig,ax = plt.subplots(figsize=(9, 9))
        for i in range(CSSData_coupled[tag].nunique()):
            ran = slice(i*CSSData_coupled['nCellu'].nunique(), (i+1)*CSSData_coupled['nCellu'].nunique())
            ax.loglog(runtime_coupled[ran],maxE_coupled[ran],markers[i],color="purple", label = f'DG-Julia, {plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize)
            ax.loglog(runtime_oca_fpi[ran],maxE_oca_fpi[ran],markers[i],color="green", label = f'OCA-FPI, {plottag} = {int(CSSData_coupled[tag].unique()[i])}',markersize=markersize)
            
            
        if os.path.exists(path + "CSS/CADETDGConvergenceCSS.csv"):
            for i in range(CSSData_coupled_DG[tag].nunique()):
                ran = slice(i*CSSData_coupled_DG['nCellu'].nunique(), (i+1)*CSSData_coupled_DG['nCellu'].nunique())
                ax.loglog(runtime_coupled_DG[ran],maxE_coupled_DG[ran],markers[i],color="orangered", label = f'DG-C++, {plottag} = {int(CSSData_coupled_DG[tag].unique()[i])}',markersize=markersize)
                if len(runtime_oca_fpi_DG)>2:
                    ax.loglog(runtime_oca_fpi_DG[ran],maxE_oca_fpi_DG[ran],markers[i],color="red", label = f'OCA-FPI, {plottag} = {int(CSSData_coupled_DG[tag].unique()[i])}',markersize=markersize)
        
        if os.path.exists(path + "CSS/CADETFVConvergenceCSS.csv"):
            ax.loglog(runtime_coupled_FV,maxE_coupled_FV, '^', color="deeppink", label = 'FV-C++',markersize=markersize)
           
        
        ax.set_xlabel('Compute time / s', fontsize=25)
        # ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.set_ylabel('Max abs error / mM', fontsize=25)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.tick_params(axis='both', which='major', labelsize=22)
        if path == "SMA/SMB/GRM/":
            ax.set_xticks([40, 60, 100, 300, 800])
            ax.set_xticklabels(['40', '60', '100', '300', '800'])
        # plt.title('LRM Langmuir')
        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
        plt.legend(fontsize=20, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, framealpha=1.0, handles=handles)
        fig.subplots_adjust(bottom=0.25)  # Adjust this value as needed
        plt.savefig(os.path.join(saveLocation,'CSS/Plot_runtime_err.svg'),format = 'svg',dpi = 1200, bbox_inches='tight')
    
    
    
    


# A function to find data and initiate the plot function 
def plot_initiator(path,no_bind = False):

        
    CADETFVdata = pd.read_csv(path + 'CADETFVConvergence.csv')
    CADETJuliadata = pd.read_csv(path + 'CADETJuliaConvergence.csv')
    if os.path.exists(path + "CADETDGConvergence.csv"):
        CADETDGdata = pd.read_csv(path + 'CADETDGConvergence.csv')
    else:
        CADETDGdata = []
        
    
    
    # Compare the speed up and convergence between DG-C++ and DG-Julia
    if type(CADETDGdata) == pd.core.frame.DataFrame:
        comparison = pd.DataFrame({'SpeedUp' : CADETDGdata['runtime_i'] / CADETJuliadata['runtime_i'], 
                                  'AvgSpeedUp' : np.mean(CADETDGdata['runtime_i'] / CADETJuliadata['runtime_i'])*np.ones(len(CADETDGdata['runtime_i'])),
                                  'DiffMae' : CADETDGdata['maxE_i'] - CADETJuliadata['maxE_i']})  
        comparison.to_csv(path + 'comparison.csv')
    
   
        
    # Plots
    if path[-4] == "G": # If using the GRM
 
        plot_convergence(CADETFVdata,CADETJuliadata[:],CADETDGdata[:],path, "GRM")
    else:
        plot_convergence(CADETFVdata,CADETJuliadata[:],CADETDGdata[:],path)


    
# plot profiles 
# Load semi-analytical solution to plot profiles
def plot_profiles(path):
    
    if path == "Linear/SMB/LRM/":
        profileData = pd.read_csv(path + 'Semi-analytical_LRM_Linear.csv', delimiter=",")
        profileCSSData = pd.read_csv(path + 'CSS/Semi-analytical_LRM_Linear_CSS.csv', delimiter=",")
    elif path == "Langmuir/SMB/LRMP/":
        profileData = pd.read_csv(path + 'Semi-analytical_LRMP_Langmuir.csv', delimiter=",")
        profileCSSData = pd.read_csv(path + 'CSS/Semi-analytical_LRMP_Langmuir_CSS.csv', delimiter=",")
    elif path == "SMA/SMB/LRMP/":
        profileData = pd.read_csv(path + 'Semi-analytical_LRMP_SMA.csv', delimiter=",")
        profileCSSData = pd.read_csv(path + 'CSS/Semi-analytical_LRMP_SMA_CSS.csv', delimiter=",")
    elif path == "SMA/SMB/GRM/":
        profileData = pd.read_csv(path + 'Semi-analytical_GRM_SMA.csv', delimiter=",")
        profileCSSData = pd.read_csv(path + 'CSS/Semi-analytical_GRM_SMA_CSS.csv', delimiter=",")
    else:
        print('get out')
    profileDataTime = pd.read_csv(path + 'Profiles_data.csv', delimiter=",")['time']
    
    
    profileData['diff'] = np.zeros(len(profileData))
    if path[0] != 'L':
        idxx = 'C2_R'
    else:
        idxx = 'C0_E'
    for i in range(1,len(profileData['diff'])):
        profileData.loc[i,'diff'] = abs(profileData[idxx][i] - profileData[idxx][i-1])
    
    # last two cycles plot
    if 'C3_E' not in profileData.columns:
        diff_1 = profileData.index[profileData['diff'] > 1]
    else:
        diff_1 = profileData.index[profileData['diff'] > 0.1]
         
         
    ############ plot profiles ############
    # Complete Extract plot
    plt.figure()
    if 'C3_E' not in profileData.columns:
        plt.plot(profileDataTime,profileData['C0_E'], label = 'A')
        plt.plot(profileDataTime,profileData['C1_E'], label = 'B')
    else:
        plt.plot(profileDataTime,profileData['C1_E'], label = 'A')
        plt.plot(profileDataTime,profileData['C2_E'], label = 'B')
        plt.plot(profileDataTime,profileData['C3_E'], label = 'C')
    plt.xlabel('Time / s')
    plt.ylabel('Concentration / mM')
    plt.legend()
    plt.savefig(os.path.join(path,'Plot_profiles_complete_E.svg'),format = 'svg',dpi = 1200)
    
    plt.figure()
    if 'C3_E' not in profileData.columns:
        plt.plot(profileDataTime,profileData['C0_R'], label = 'A')
        plt.plot(profileDataTime,profileData['C1_R'], label = 'B')
    else:
        plt.plot(profileDataTime,profileData['C1_R'], label = 'A')
        plt.plot(profileDataTime,profileData['C2_R'], label = 'B')
        plt.plot(profileDataTime,profileData['C3_R'], label = 'C')
    plt.xlabel('Time / s')
    plt.ylabel('Concentration / mM')
    plt.legend()
    plt.savefig(os.path.join(path,'Plot_profiles_complete_R.svg'),format = 'svg',dpi = 1200)
    
   
    
    plt.figure()
    if 'C3_E' not in profileData.columns:
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C0_E'][diff_1[-2]:], label = 'A')
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C1_E'][diff_1[-2]:], label = 'B')
    else:
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C1_E'][diff_1[-2]:], label = 'A')
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C2_E'][diff_1[-2]:], label = 'B')
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C3_E'][diff_1[-2]:], label = 'C')
    plt.xlabel('Time / s')
    plt.ylabel('Concentration / mM')
    plt.legend()
    plt.savefig(os.path.join(path,'Plot_profiles_last2_E.svg'),format = 'svg',dpi = 1200)
    
    
    
    plt.figure()
    if 'C3_E' not in profileData.columns:
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C0_R'][diff_1[-2]:], label = 'A')
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C1_R'][diff_1[-2]:], label = 'B')
    else:
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C1_R'][diff_1[-2]:], label = 'A')
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C2_R'][diff_1[-2]:], label = 'B')
        plt.plot(profileDataTime[diff_1[-2]:],profileData['C3_R'][diff_1[-2]:], label = 'C')
    plt.xlabel('Time / s')
    plt.ylabel('Concentration / mM')
    plt.legend()
    plt.savefig(os.path.join(path,'Plot_profiles_last2_R.svg'),format = 'svg',dpi = 1200)
    
    
    # Combined plot
    fig,ax = plt.subplots(2,2,figsize=(24, 10)) #figsize=(11.5, 10)
    #Counting number of components as number of x-1
    if 'C3_E' not in profileData.columns:
        ax[0,0].plot(profileDataTime, profileData['C0_E'], label='A')
        ax[0,0].plot(profileDataTime,profileData['C1_E'], label = 'B')
        
        ax[0,1].plot(profileDataTime,profileData['C0_R'], label = 'A')
        ax[0,1].plot(profileDataTime,profileData['C1_R'], label = 'B')
        
        ax[1,0].plot(profileDataTime[diff_1[-2]:],profileData['C0_E'][diff_1[-2]:], label = 'A')
        ax[1,0].plot(profileDataTime[diff_1[-2]:],profileData['C1_E'][diff_1[-2]:], label = 'B')
        
        ax[1,1].plot(profileDataTime[diff_1[-2]:],profileData['C0_R'][diff_1[-2]:], label = 'A')
        ax[1,1].plot(profileDataTime[diff_1[-2]:],profileData['C1_R'][diff_1[-2]:], label = 'B')
        
        
    else:
        ax[0,0].plot(profileDataTime,profileData['C1_E'], label = 'A')
        ax[0,0].plot(profileDataTime,profileData['C2_E'], label = 'B')
        ax[0,0].plot(profileDataTime,profileData['C3_E'], label = 'C')
        
        ax[0,1].plot(profileDataTime,profileData['C1_R'], label = 'A')
        ax[0,1].plot(profileDataTime,profileData['C2_R'], label = 'B')
        ax[0,1].plot(profileDataTime,profileData['C3_R'], label = 'C')
        
        ax[1,0].plot(profileDataTime[diff_1[-2]:],profileData['C1_E'][diff_1[-2]:], label = 'A')
        ax[1,0].plot(profileDataTime[diff_1[-2]:],profileData['C2_E'][diff_1[-2]:], label = 'B')
        ax[1,0].plot(profileDataTime[diff_1[-2]:],profileData['C3_E'][diff_1[-2]:], label = 'C')
        
        ax[1,1].plot(profileDataTime[diff_1[-2]:],profileData['C1_R'][diff_1[-2]:], label = 'A')
        ax[1,1].plot(profileDataTime[diff_1[-2]:],profileData['C2_R'][diff_1[-2]:], label = 'B')
        ax[1,1].plot(profileDataTime[diff_1[-2]:],profileData['C3_R'][diff_1[-2]:], label = 'C')
        
        
    fig.text(0.075, 0.5, "Concentration / mM", ha="center", va="center", rotation=90, fontsize=30)
    ax[1,0].set_xlabel('Time / s', fontsize=30)
    ax[1,1].set_xlabel('Time / s', fontsize=30)
    ax[0,0].legend(fontsize=25)
    ax[0,0].tick_params(axis='both', labelsize=22)
    ax[0,1].tick_params(axis='both', labelsize=22)
    ax[1,0].tick_params(axis='both', labelsize=22)
    ax[1,1].tick_params(axis='both', labelsize=22)
    
    
    ax[0,0].set_title('Extract', fontsize=30)
    ax[0,1].set_title('Raffinate', fontsize=30)

    plt.savefig(os.path.join(path,'Plot_profiles_combined.svg'),format = 'svg',dpi = 1200)
    
    
    
    
    
    
    # CSS combined profiles    
    
    fig,ax = plt.subplots(2,1,figsize=(10, 7)) #figsize=(11.5, 10)
    #Counting number of components as number of x-1
    if 'C3_E' not in profileCSSData.columns:        
        ax[0].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C0_E'][:], label = 'A')
        ax[0].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C1_E'][:], label = 'B')
        
        ax[1].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C0_R'][:], label = 'A')
        ax[1].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C1_R'][:], label = 'B')
        
    else:
        ax[0].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C1_E'], label = 'A')
        ax[0].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C2_E'], label = 'B')
        ax[0].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C3_E'], label = 'C')
        
        ax[1].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C1_R'], label = 'A')
        ax[1].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C2_R'], label = 'B')
        ax[1].plot(profileDataTime[:len(profileCSSData)], profileCSSData['C3_R'], label = 'C')
        
        
    fig.text(0.075, 0.5, "Concentration / mM", ha="center", va="center", rotation=90, fontsize=20)
    ax[1].set_xlabel('Time / s', fontsize=20)
    ax[0].legend(fontsize=20)
    
    ax[0].tick_params(axis='both', labelsize=20)
    ax[1].tick_params(axis='both', labelsize=20)
    
    ax[0].set_title('Extract', fontsize=20)
    ax[1].set_title('Raffinate', fontsize=20)
    fig.subplots_adjust(left=0.15, bottom=None, right=None, top=None, wspace=None, hspace=.4)
    # fig.tight_layout(h_pad=2)
    plt.savefig(os.path.join(path,'CSS/Plot_profiles_combined.svg'),format = 'svg',dpi = 1200)


#%%
# Run the convergence function to plot the Benchmarks 
path = "Linear/SMB/LRM/"
plot_initiator(path)
# plot_profiles(path)


path = "Langmuir/SMB/LRMP/"
plot_initiator(path)
# plot_profiles(path)


path = "SMA/SMB/LRMP/"
plot_initiator(path)
# plot_profiles(path)


path = "SMA/SMB/GRM/"
plot_initiator(path)
# plot_profiles(path)




