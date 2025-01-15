# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:12:37 2023

@author: jespfra
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import timeit
import pandas as pd

from cadet import Cadet
# Cadet.cadet_path = r'C:\Users\pbzit\source\Test\out\install\aRELEASE\bin\cadet-cli'
Cadet.cadet_path = r'C:\Users\pbzit\source\Test\out\install\aRELEASE\bin\cadet-cli'


#%% General model options


def model(ncol, polydeg, polydegpore, exactInt=True, state_y=[0], state_ydot=[0]):
    #Setting up the model
    smb_model = Cadet()


    #Speciy number of unit operations: input, column and output, 3
    smb_model.root.input.model.nunits = 8
    
    #Specify # of components (salt,proteins)
    n_comp  = 2
    
    #First unit operation: inlet
    ## Feed
    smb_model.root.input.model.unit_000.unit_type = 'INLET'
    smb_model.root.input.model.unit_000.ncomp = n_comp
    smb_model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'
    
    
    ## Eluent
    smb_model.root.input.model.unit_001.unit_type = 'INLET'
    smb_model.root.input.model.unit_001.ncomp = n_comp
    smb_model.root.input.model.unit_001.inlet_type = 'PIECEWISE_CUBIC_POLY'
    
    ## Extract 
    smb_model.root.input.model.unit_002.ncomp = n_comp
    smb_model.root.input.model.unit_002.unit_type = 'OUTLET'
    
    ## Raffinate
    smb_model.root.input.model.unit_003.ncomp = n_comp
    smb_model.root.input.model.unit_003.unit_type = 'OUTLET'
    
    
    ## Columns
    smb_model.root.input.model.unit_004.unit_type = 'LUMPED_RATE_MODEL_WITH_PORES'
    smb_model.root.input.model.unit_004.ncomp = n_comp 

    ## Geometry
    smb_model.root.input.model.unit_004.col_porosity = 0.6
    smb_model.root.input.model.unit_004.par_porosity = 0.2
    smb_model.root.input.model.unit_004.col_dispersion = 1e-5
    smb_model.root.input.model.unit_004.col_length = 0.25
    smb_model.root.input.model.unit_004.cross_section_area = 1 #From Lubke2007, is not important
    smb_model.root.input.model.unit_004.film_diffusion = [3.3e-3, 3.3e-3]
    smb_model.root.input.model.unit_004.par_radius = 1.0e-4


    #Isotherm specification
    smb_model.root.input.model.unit_004.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
    smb_model.root.input.model.unit_004.adsorption.is_kinetic = False    # Kinetic binding
    smb_model.root.input.model.unit_004.adsorption.MCL_KA = [0.1, 0.05]      # m^3 / (mol * s)   (mobile phase)
    smb_model.root.input.model.unit_004.adsorption.MCL_KD = [1,1]      # 1 / s (desorption)
    smb_model.root.input.model.unit_004.adsorption.MCL_QMAX = [10,10]      # 1 / s (desorption)
    
    #To write out last output to check for steady state
    smb_model.root.input['return'].WRITE_SOLUTION_LAST = True
    
    #Initial conditions
    if len(state_y) == 1:
        smb_model.root.input.model.unit_004.init_c = [0,0]
        smb_model.root.input.model.unit_004.init_q = [0,0] 
    else:
        smb_model.root.input.model.init_state_y = state_y
        smb_model.root.input.model.init_state_ydot = state_ydot
    
    
    
    #Discretization
    smb_model.root.input.model.unit_004.discretization.ncol = ncol
    smb_model.root.input.model.unit_004.discretization.npar = ncol
    
    smb_model.root.input.model.unit_004.discretization.nbound = np.ones(n_comp,dtype=int)
    
    
    smb_model.root.input.model.unit_004.discretization.par_disc_type = 'EQUIDISTANT_PAR'    
    smb_model.root.input.model.unit_004.discretization.use_analytic_jacobian = 1
    smb_model.root.input.model.unit_004.discretization.reconstruction = 'WENO'
    smb_model.root.input.model.unit_004.discretization.gs_type = 1
    smb_model.root.input.model.unit_004.discretization.max_krylov = 0
    smb_model.root.input.model.unit_004.discretization.max_restarts = 10
    smb_model.root.input.model.unit_004.discretization.schur_safety = 1.0e-8

    smb_model.root.input.model.unit_004.discretization.weno.boundary_model = 0
    smb_model.root.input.model.unit_004.discretization.weno.weno_eps = 1e-10
    smb_model.root.input.model.unit_004.discretization.weno.weno_order = 3
    
    ### Copy column models
    smb_model.root.input.model.unit_005 = smb_model.root.input.model.unit_004
    smb_model.root.input.model.unit_006 = smb_model.root.input.model.unit_004
    smb_model.root.input.model.unit_007 = smb_model.root.input.model.unit_004
    
    
    #% Input and connections
    n_cycles = 1
    switch_time = ts #s
    
    #Sections
    smb_model.root.input.solver.sections.nsec = 4*n_cycles
    smb_model.root.input.solver.sections.section_times = [0]
    for i in range(n_cycles):
        smb_model.root.input.solver.sections.section_times.append((4*i+1)*switch_time)
        smb_model.root.input.solver.sections.section_times.append((4*i+2)*switch_time)
        smb_model.root.input.solver.sections.section_times.append((4*i+3)*switch_time)
        smb_model.root.input.solver.sections.section_times.append((4*i+4)*switch_time)    
    
    ## Feed and Eluent concentration
    smb_model.root.input.model.unit_000.sec_000.const_coeff = [CfA, CfB] #Inlet flowrate concentration
    smb_model.root.input.model.unit_001.sec_000.const_coeff = [0, 0] #Desorbent stream
    
    
    #Connections
    smb_model.root.input.model.connections.nswitches = 4
    
    smb_model.root.input.model.connections.switch_000.section = 0
    smb_model.root.input.model.connections.switch_000.connections = [
        4, 5, -1, -1, Q4,#flowrates, Q, m3/s
        5, 6, -1, -1, Q4,
        6, 7, -1, -1, Q2,
        7, 4, -1, -1, Q2,
        0, 4, -1, -1, QF,
        1, 6, -1, -1, QD,
        4, 3, -1, -1, QR,
        6, 2, -1, -1, QE
    ]
    
    smb_model.root.input.model.connections.switch_001.section = 1
    smb_model.root.input.model.connections.switch_001.connections = [
        4,	5,	-1,	-1,	Q2,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q4,
        6,	7,	-1,	-1,	Q4,
        7,	4,	-1,	-1,	Q2,
        0,	5,	-1,	-1,	QF,
        1,	7,	-1,	-1,	QD,
        5,	3,	-1,	-1,	QR,
        7,	2,	-1,	-1,	QE
    ]
    
    smb_model.root.input.model.connections.switch_002.section = 2
    smb_model.root.input.model.connections.switch_002.connections = [
        4,	5,	-1,	-1,	Q2,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q2,
        6,	7,	-1,	-1,	Q4,
        7,	4,	-1,	-1,	Q4,
        0,	6,	-1,	-1,	QF,
        1,	4,	-1,	-1,	QD,
        6,	3,	-1,	-1,	QR,
        4,	2,	-1,	-1,	QE
    ]
    
    smb_model.root.input.model.connections.switch_003.section = 3
    smb_model.root.input.model.connections.switch_003.connections = [
        4,	5,	-1,	-1,	Q4,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q2,
        6,	7,	-1,	-1,	Q2,
        7,	4,	-1,	-1,	Q4,
        0,	7,	-1,	-1,	QF,
        1,	5,	-1,	-1,	QD,
        7,	3,	-1,	-1,	QR,
        5,	2,	-1,	-1,	QE
    ]
    
    #solution times
    tt = [0.0]
    for i in range(1, n_cycles*4 + 1):
        ttt = [tt[-1] + val for val in [1e-10, 5e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]]
        tt.extend(ttt + list(np.linspace(0.1 + (i-1)*ts, i*ts, 10*ts)))

    smb_model.root.input.solver.user_solution_times = tt
    
    
    #Time 
    # Tolerances for the time integrator
    smb_model.root.input.solver.time_integrator.abstol = 1e-12 #absolute tolerance
    smb_model.root.input.solver.time_integrator.algtol = 1e-10
    smb_model.root.input.solver.time_integrator.reltol = 1e-10 #Relative tolerance
    smb_model.root.input.solver.time_integrator.init_step_size = 1e-10
    smb_model.root.input.solver.time_integrator.max_steps = 1000000
    
    
    
    #Solver options in general (not only for column although the same)
    smb_model.root.input.model.solver.gs_type = 1
    smb_model.root.input.model.solver.max_krylov = 0
    smb_model.root.input.model.solver.max_restarts = 10
    smb_model.root.input.model.solver.schur_safety = 1e-8
    
    # Number of cores for parallel simulation
    smb_model.root.input.solver.nthreads = 1
    
    
    #Specify which results we want to return
    # Return data
    smb_model.root.input['return'].split_components_data = 0
    smb_model.root.input['return'].split_ports_data = 0
    smb_model.root.input['return'].unit_000.write_solution_bulk = 0
    smb_model.root.input['return'].unit_000.write_solution_inlet = 0
    smb_model.root.input['return'].unit_000.write_solution_outlet = 0
    smb_model.root.input['return'].unit_002.write_solution_bulk = 0
    smb_model.root.input['return'].unit_002.write_solution_inlet = 0
    smb_model.root.input['return'].unit_002.write_solution_outlet = 1
    
    
    # Copy settings to the other unit operations
    smb_model.root.input['return'].unit_001 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_003 = smb_model.root.input['return'].unit_002
    smb_model.root.input['return'].unit_004 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_005 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_006 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_007 = smb_model.root.input['return'].unit_000
    
    
    #Saving data
    smb_model.filename = 'ModelFV.h5'
    smb_model.save()
    
    return smb_model
    
                       
    

#Import analytical solution from Jan
c_analytical = pd.read_csv('Semi-analytical_LRMP_Langmuir_CSS.csv')


data = []
runtime = []
maxE_e = []
runtime = []
DOF = []

ncols = [8, 16, 32, 64]


# Determining flows 
CfA = 10                             #mol/m3
CfB = 10                             #mol/m3

ka = [0.1, 0.05]
kd = [1.0,1.0]
qmax = [10.0,10.0]


m1 = 1.2*ka[0]*qmax[0]
m2 = 1.2*ka[1]*qmax[1]
m3 = 0.9*ka[0]*qmax[0]
m4 = 0.8*ka[1]*qmax[1]
ColLength = 0.25
ts = 180 #s

eps_c = 0.6
eps_p = 0.2
eps_t = eps_c + (1-eps_c) * eps_p

Q1 = -(m1*eps_t - m1 - eps_t)*ColLength/ts
Q2 = -(m2*eps_t - m2 - eps_t)*ColLength/ts
Q3 = -(m3*eps_t - m3 - eps_t)*ColLength/ts
Q4 = -(m4*eps_t - m4 - eps_t)*ColLength/ts  
QD = -ColLength*(m1*eps_t - m4*eps_t - m1 + m4)/ts
QE = -ColLength*(m1*eps_t - m2*eps_t - m1 + m2)/ts
QF = ColLength*(m2*eps_t - m3*eps_t - m2 + m3)/ts
QR = -ColLength*(m3*eps_t - m4*eps_t - m3 + m4)/ts



exec(open('../../../../evaluate_convergence_css.py').read())
Evaluate_CADETDG_css(model=model, c_analytical=c_analytical, nCells=ncols, polyDeg=[0], polyDegPore=0, nCellsPar=1, max_iter=150*8, css_tol=1e-6)

import os
os.remove('ModelFV.h5')