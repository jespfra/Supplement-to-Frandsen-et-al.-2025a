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

def model(nelem,polydeg,exactInt):
    #Setting up the model
    smb_model = Cadet()


    #Speciy number of unit operations: input, column and output, 3
    smb_model.root.input.model.nunits = 12
    
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
    smb_model.root.input.model.unit_004.unit_type = 'LUMPED_RATE_MODEL_WITHOUT_PORES'
    smb_model.root.input.model.unit_004.ncomp = n_comp 

    ## Geometry
    smb_model.root.input.model.unit_004.total_porosity = 0.38
    smb_model.root.input.model.unit_004.col_dispersion = 3.81e-6
    smb_model.root.input.model.unit_004.col_length = 5.36e-1
    smb_model.root.input.model.unit_004.cross_section_area = 5.31e-4 #From Lubke2007, is not important


    #Isotherm specification
    smb_model.root.input.model.unit_004.adsorption_model = 'LINEAR'
    smb_model.root.input.model.unit_004.adsorption.is_kinetic = False    # Kinetic binding
    smb_model.root.input.model.unit_004.adsorption.LIN_KA = [0.54, 0.28]      # m^3 / (mol * s)   (mobile phase)
    smb_model.root.input.model.unit_004.adsorption.LIN_KD = [1,1]      # 1 / s (desorption)
    
    #Initial conditions
    smb_model.root.input.model.unit_004.init_c = [0,0]
    smb_model.root.input.model.unit_004.init_q = [0,0] #salt starts at max capacity
    
    
    ### Grid cells in column and particle: the most important ones - ensure grid-independent solutions
    smb_model.root.input.model.unit_004.discretization.SPATIAL_METHOD = "DG"
    smb_model.root.input.model.unit_004.discretization.nelem = nelem 
    
    #Polynomial order 
    smb_model.root.input.model.unit_004.discretization.polydeg = polydeg
    smb_model.root.input.model.unit_004.discretization.exact_integration = exactInt

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
    smb_model.root.input.model.unit_008 = smb_model.root.input.model.unit_004
    smb_model.root.input.model.unit_009 = smb_model.root.input.model.unit_004
    smb_model.root.input.model.unit_010 = smb_model.root.input.model.unit_004
    smb_model.root.input.model.unit_011 = smb_model.root.input.model.unit_004
    
    #To write out last output to check for steady state
    smb_model.root.input['return'].WRITE_SOLUTION_LAST = True



    #% Input and connections
    n_cycles = 10
    switch_time = ts #s
    
    #Sections
    smb_model.root.input.solver.sections.nsec = 8*n_cycles
    smb_model.root.input.solver.sections.section_times = [0]
    for i in range(n_cycles):
        smb_model.root.input.solver.sections.section_times.append((8*i+1)*switch_time)
        smb_model.root.input.solver.sections.section_times.append((8*i+2)*switch_time)
        smb_model.root.input.solver.sections.section_times.append((8*i+3)*switch_time)
        smb_model.root.input.solver.sections.section_times.append((8*i+4)*switch_time)    
        smb_model.root.input.solver.sections.section_times.append((8*i+5)*switch_time)    
        smb_model.root.input.solver.sections.section_times.append((8*i+6)*switch_time)    
        smb_model.root.input.solver.sections.section_times.append((8*i+7)*switch_time)    
        smb_model.root.input.solver.sections.section_times.append((8*i+8)*switch_time)    
    
    ## Feed and Eluent concentration
    smb_model.root.input.model.unit_000.sec_000.const_coeff = [2.78,2.78] #Inlet flowrate concentration
    smb_model.root.input.model.unit_001.sec_000.const_coeff = [0, 0] #Desorbent stream
    
    
    #Connections
    smb_model.root.input.model.connections.nswitches = 8
    
    smb_model.root.input.model.connections.switch_000.section = 0
    smb_model.root.input.model.connections.switch_000.connections =[
        4, 5, -1, -1, Q3,#flowrates, Q, m3/s
        5, 6, -1, -1, Q4,
        6, 7, -1, -1, Q4,
        7, 8, -1, -1, Q4,
        8, 9, -1, -1, Q1,
        9, 10, -1, -1, Q2,
        10, 11, -1, -1, Q2,
        11, 4, -1, -1, Q2,
        0, 4, -1, -1, QF,
        1, 8, -1, -1, QD,
        5, 3, -1, -1, QR,
        9, 2, -1, -1, QE
    ]

    smb_model.root.input.model.connections.switch_001.section = 1
    smb_model.root.input.model.connections.switch_001.connections =[
        4, 5, -1, -1, Q2,#flowrates, Q, m3/s
        5, 6, -1, -1, Q3,
        6, 7, -1, -1, Q4,
        7, 8, -1, -1, Q4,
        8, 9, -1, -1, Q4,
        9, 10, -1, -1, Q1,
        10, 11, -1, -1, Q2,
        11, 4, -1, -1, Q2,
        0, 5, -1, -1, QF,
        1, 9, -1, -1, QD,
        6, 3, -1, -1, QR,
        10, 2, -1, -1, QE
    ]

    smb_model.root.input.model.connections.switch_002.section = 2
    smb_model.root.input.model.connections.switch_002.connections =[
        4, 5, -1, -1, Q2,#flowrates, Q, m3/s
        5, 6, -1, -1, Q2,
        6, 7, -1, -1, Q3,
        7, 8, -1, -1, Q4,
        8, 9, -1, -1, Q4,
        9, 10, -1, -1, Q4,
        10, 11, -1, -1, Q1,
        11, 4, -1, -1, Q2,
        0, 6, -1, -1, QF,
        1, 10, -1, -1, QD,
        7, 3, -1, -1, QR,
        11, 2, -1, -1, QE
    ]

    smb_model.root.input.model.connections.switch_003.section = 3
    smb_model.root.input.model.connections.switch_003.connections =[
        4, 5, -1, -1, Q2,#flowrates, Q, m3/s
        5, 6, -1, -1, Q2,
        6, 7, -1, -1, Q2,
        7, 8, -1, -1, Q3,
        8, 9, -1, -1, Q4,
        9, 10, -1, -1, Q4,
        10, 11, -1, -1, Q4,
        11, 4, -1, -1, Q1,
        0, 7, -1, -1, QF,
        1, 11, -1, -1, QD,
        8, 3, -1, -1, QR,
        4, 2, -1, -1, QE
    ]


    smb_model.root.input.model.connections.switch_004.section = 4
    smb_model.root.input.model.connections.switch_004.connections =[
        4, 5, -1, -1, Q1,#flowrates, Q, m3/s
        5, 6, -1, -1, Q2,
        6, 7, -1, -1, Q2,
        7, 8, -1, -1, Q2,
        8, 9, -1, -1, Q3,
        9, 10, -1, -1, Q4,
        10, 11, -1, -1, Q4,
        11, 4, -1, -1, Q4,
        0, 8, -1, -1, QF,
        1, 4, -1, -1, QD,
        9, 3, -1, -1, QR,
        5, 2, -1, -1, QE
    ]


    smb_model.root.input.model.connections.switch_005.section = 5
    smb_model.root.input.model.connections.switch_005.connections =[
        4, 5, -1, -1, Q4,#flowrates, Q, m3/s
        5, 6, -1, -1, Q1,
        6, 7, -1, -1, Q2,
        7, 8, -1, -1, Q2,
        8, 9, -1, -1, Q2,
        9, 10, -1, -1, Q3,
        10, 11, -1, -1, Q4,
        11, 4, -1, -1, Q4,
        0, 9, -1, -1, QF,
        1, 5, -1, -1, QD,
        10, 3, -1, -1, QR,
        6, 2, -1, -1, QE
    ]


    smb_model.root.input.model.connections.switch_006.section = 6
    smb_model.root.input.model.connections.switch_006.connections =[
        4, 5, -1, -1, Q4,#flowrates, Q, m3/s
        5, 6, -1, -1, Q4,
        6, 7, -1, -1, Q1,
        7, 8, -1, -1, Q2,
        8, 9, -1, -1, Q2,
        9, 10, -1, -1, Q2,
        10, 11, -1, -1, Q3,
        11, 4, -1, -1, Q4,
        0, 10, -1, -1, QF,
        1, 6, -1, -1, QD,
        11, 3, -1, -1, QR,
        7, 2, -1, -1, QE
    ]


    smb_model.root.input.model.connections.switch_007.section = 7
    smb_model.root.input.model.connections.switch_007.connections =[
        4, 5, -1, -1, Q4,#flowrates, Q, m3/s
        5, 6, -1, -1, Q4,
        6, 7, -1, -1, Q4,
        7, 8, -1, -1, Q1,
        8, 9, -1, -1, Q2,
        9, 10, -1, -1, Q2,
        10, 11, -1, -1, Q2,
        11, 4, -1, -1, Q3,
        0, 11, -1, -1, QF,
        1, 7, -1, -1, QD,
        4, 3, -1, -1, QR,
        8, 2, -1, -1, QE
    ]

    
    #solution times
    smb_model.root.input.solver.user_solution_times = np.linspace(0, n_cycles*8*switch_time, int(n_cycles*8*switch_time)+1)
    
    
    
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
    smb_model.root.input.solver.consistent_init_mode = 5 #necessary specifically for this sim
    smb_model.root.input.solver.time_integrator.USE_MODIFIED_NEWTON = 1
    
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
    smb_model.root.input['return'].unit_008 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_009 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_010 = smb_model.root.input['return'].unit_000
    smb_model.root.input['return'].unit_011 = smb_model.root.input['return'].unit_000


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
                       'C0_R': smb_model.root.output.solution.unit_003.solution_outlet[:,0],
                       'C1_R': smb_model.root.output.solution.unit_003.solution_outlet[:,1],})
        return df, smb_model.root.meta.time_sim, stop - start
    else:
        print(data)
        return [], [], []
    
    



#Import analytical solution from Jan
c_analytical = pd.read_csv('Semi-analytical_LRM_Linear.csv')


nCells = [2,4,8,16]
polyDeg = [4,5,6]
ncolumns = 8

ts = 1552  
QD = 4.14e-8
QE = 3.48e-8
QF = 2.00e-8
QR = 2.66e-8
Q2 = 1.05e-7
Q3 = Q2 + QF
Q4 = Q3 - QR
Q1 = Q4 + QD

# df, rt, rt1 = model(4,4,1)

exec(open('../../../benchmark_runner.py').read())
runCadetDG("LRM", c_analytical, polyDeg, nCells, ncolumns)

import os
os.remove("ModelDG.h5")



