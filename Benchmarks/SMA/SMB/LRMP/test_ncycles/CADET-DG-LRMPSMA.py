# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:12:37 2023

@author: jespfra
"""


import numpy as np
import matplotlib.pyplot as plt
import timeit
import time
import pandas as pd

from cadet import Cadet
Cadet.cadet_path = r'C:\Users\pbzit\source\Test\out\install\aRELEASE\bin\cadet-cli'

#%% General model options

def model(ncol,polydeg,ncycles):
    #Setting up the model
    smb_model = Cadet()


    #Speciy number of unit operations: input, column and output, 3
    smb_model.root.input.model.nunits = 8
    
    #Specify # of components (salt,proteins)
    n_comp  = 4
    
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
    smb_model.root.input.model.unit_004.col_porosity = 0.37
    smb_model.root.input.model.unit_004.par_porosity = 0.75
    smb_model.root.input.model.unit_004.col_dispersion = 1e-5
    smb_model.root.input.model.unit_004.col_length = 0.014
    #smb_model.root.input.model.unit_004.total_porosity = 0.37+0.75*(1-0.37)
    smb_model.root.input.model.unit_004.cross_section_area = 1
    smb_model.root.input.model.unit_004.film_diffusion = [3.3e-3,3.3e-3,3.3e-3,3.3e-3]
    smb_model.root.input.model.unit_004.par_radius = 4.5e-5


    #Isotherm specification
    smb_model.root.input.model.unit_004.adsorption_model = 'STERIC_MASS_ACTION'
    smb_model.root.input.model.unit_004.adsorption.is_kinetic = 1    # Kinetic binding
    smb_model.root.input.model.unit_004.adsorption.sma_ka = [0, 35.5e-3, 1.59e-3, 7.70e-3]      # m^3 / (mol * s)   (mobile phase)
    smb_model.root.input.model.unit_004.adsorption.sma_kd = [0, 1, 1,1]      # 1 / s (desorption)
    smb_model.root.input.model.unit_004.adsorption.sma_lambda = 1200    #max binding capacity
    smb_model.root.input.model.unit_004.adsorption.sma_nu = [0.0, 4.7, 5.29, 3.7]    #characteristic charge
    smb_model.root.input.model.unit_004.adsorption.sma_sigma = [0.0, 11.83, 10.6, 10.00] #shielding factor

    #Initial conditions
    smb_model.root.input.model.unit_004.init_c = [282,0,0,0]
    smb_model.root.input.model.unit_004.init_q = [1200,0,0,0] #salt starts at max capacity
    
    
    ### Grid cells in column and particle: the most important ones - ensure grid-independent solutions
    smb_model.root.input.model.unit_004.discretization.SPATIAL_METHOD = "DG"
    smb_model.root.input.model.unit_004.discretization.nelem = ncol 
    
    #Polynomial order 
    smb_model.root.input.model.unit_004.discretization.polydeg = polydeg
    smb_model.root.input.model.unit_004.discretization.exact_integration = 0

    ### Bound states - for zero the compound does not bind, >1 = multiple binding sites
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
    
    #To write out last output to check for steady state
    smb_model.root.input['return'].WRITE_SOLUTION_LAST = True



    #% Input and connections
    n_cycles = ncycles
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
    smb_model.root.input.model.unit_000.sec_000.const_coeff = [Cs0, CfA, CfB, CfC] #Inlet flowrate concentration
    smb_model.root.input.model.unit_001.sec_000.const_coeff = [Cs0, 0, 0, 0] #Desorbent stream
    
    
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
    smb_model.root.input.solver.user_solution_times = np.linspace(0, n_cycles*4*switch_time, int(n_cycles*4*switch_time*10)+1)
    
    
    #Time 
    # Tolerances for the time integrator
    smb_model.root.input.solver.time_integrator.abstol = 1e-12 #absolute tolerance
    smb_model.root.input.solver.time_integrator.algtol = 1e-10
    smb_model.root.input.solver.time_integrator.reltol = 1e-10 #Relative tolerance
    smb_model.root.input.solver.time_integrator.init_step_size = 1e-10
    smb_model.root.input.solver.time_integrator.max_steps = 1000000
    smb_model.root.input.solver.time_integrator.USE_MODIFIED_NEWTON = 1
    
    
    
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


    #Saving data - do not touch
    smb_model.filename = 'ModelDG.h5'
    smb_model.save()

    #run model
    start = timeit.default_timer()
    data = smb_model.run()
    stop = timeit.default_timer()
    
    if data.returncode == 0:
        print("Simulation completed successfully")
        smb_model.load()   
        df = pd.DataFrame({'Time': smb_model.root.output.solution.solution_times,
                       'C0_E': smb_model.root.output.solution.unit_002.solution_outlet[:,0],
                       'C1_E': smb_model.root.output.solution.unit_002.solution_outlet[:,1],
                       'C2_E': smb_model.root.output.solution.unit_002.solution_outlet[:,2],
                       'C3_E': smb_model.root.output.solution.unit_002.solution_outlet[:,3],
                       'C0_R': smb_model.root.output.solution.unit_003.solution_outlet[:,0],
                       'C1_R': smb_model.root.output.solution.unit_003.solution_outlet[:,1],
                       'C2_R': smb_model.root.output.solution.unit_003.solution_outlet[:,2],
                       'C3_R': smb_model.root.output.solution.unit_003.solution_outlet[:,3],})
        return df, smb_model.root.meta.time_sim, stop - start
    else:
        print(data)
        return [],[],[]
    
    


ncolumns = 4
# setting 1: nCells = 4, polyDeg = 4
# setting 2: nCells = 16, polyDeg = 5
# setting 1: nCells = 32, polyDeg = 6

# Determining flows 
Cs0, CfA, CfB, CfC = 282, 1, 1, 1    #mol/m3
Lambda = 1200                    #mol/m3
Keq = [35.5e-3, 1.59e-3, 7.70e-3,]           #Equilibrium constant, A, B, C
v = [4.7, 5.29, 3.7]                      #Characteristic charge


m1 = Keq[1]*(Lambda/Cs0)**v[1]*1.5 
m2 = 5.0 #
m3 = 10.0 # 
m4 = 0.7256 #
ColLength = 1.4e-2
ts = 180 # Switch time, s


# Geometry
eps_p = 0.75 #Column porosity
eps_c = 0.37 #particle porosity
eps_t = eps_c + (1-eps_c) * eps_p

Q1 = -(m1*eps_t - m1 - eps_t)*ColLength/ts
Q2 = -(m2*eps_t - m2 - eps_t)*ColLength/ts
Q3 = -(m3*eps_t - m3 - eps_t)*ColLength/ts
Q4 = -(m4*eps_t - m4 - eps_t)*ColLength/ts
QD = -ColLength*(m1*eps_t - m4*eps_t - m1 + m4)/ts
QE = -ColLength*(m1*eps_t - m2*eps_t - m1 + m2)/ts
QF = ColLength*(m2*eps_t - m3*eps_t - m2 + m3)/ts
QR = -ColLength*(m3*eps_t - m4*eps_t - m3 + m4)/ts

# Settings configurations
settings = [
    {'nCells': 2, 'polyDeg': 4, 'setting_name': 'Setting 1'},
    {'nCells': 4, 'polyDeg': 5, 'setting_name': 'Setting 2'},
    {'nCells': 8, 'polyDeg': 6, 'setting_name': 'Setting 3'}
]

# Number of cycles to evaluate
n_cycles = [5, 10, 15, 20, 30, 50]

# List to store results
results = []

# Loop through each setting and each cycle count
for setting in settings:
    for ncycle in n_cycles:
        # Evaluate the model for each setting and cycle
        df, rtime, rtime_tot = model(setting['nCells'], setting['polyDeg'], ncycle)
        
        # Append the result as a dictionary
        results.append({
            'Setting': setting['setting_name'],
            'nCells': setting['nCells'],
            'polyDeg': setting['polyDeg'],
            'nCycles': ncycle,
            'rtime': rtime,
            'rtime_tot':rtime_tot
        })

# Convert results list into a pandas DataFrame
df_results = pd.DataFrame(results)
df_results.to_csv('CADETDGConvergence.csv')


import os
os.remove('ModelDG.h5')
