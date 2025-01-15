# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 15:12:37 2023

@author: jespfra
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy


from cadet import Cadet
Cadet.cadet_path = r'C:\Users\pbzit\source\Test\out\install\aRELEASE\bin\cadet-cli'
cadet_path = Cadet.cadet_path
#%% General model options

def Cyclic_model1(nelem,polydeg,exactInt):
    
    
    # Setting up the model
    Cyclic_model = Cadet()

    # Speciy number of unit operations: input, column and output, 3
    Cyclic_model.root.input.model.nunits = 4

    # Specify # of components (salt,proteins)
    n_comp = 1

    # First unit operation: inlet
    # Source 1
    Cyclic_model.root.input.model.unit_000.unit_type = 'INLET'
    Cyclic_model.root.input.model.unit_000.ncomp = n_comp
    Cyclic_model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

    # Sink
    Cyclic_model.root.input.model.unit_003.ncomp = n_comp
    Cyclic_model.root.input.model.unit_003.unit_type = 'OUTLET'

    # Unit LRMP2
    Cyclic_model.root.input.model.unit_001.unit_type = 'LUMPED_RATE_MODEL_WITH_PORES'
    Cyclic_model.root.input.model.unit_001.ncomp = n_comp

    # Geometry
    Cyclic_model.root.input.model.unit_001.col_porosity = 0.37
    Cyclic_model.root.input.model.unit_001.par_porosity = 0.75
    Cyclic_model.root.input.model.unit_001.col_dispersion = 2e-7
    Cyclic_model.root.input.model.unit_001.col_length = 1.4e-2
    # From Lubke2007, is not important
    Cyclic_model.root.input.model.unit_001.cross_section_area = 1
    Cyclic_model.root.input.model.unit_001.film_diffusion = 6.9e-6
    Cyclic_model.root.input.model.unit_001.par_radius = 45e-6
    LRMP_Q3 = 3.45*1e-2 / 60 * 0.37

    # Isotherm specification
    Cyclic_model.root.input.model.unit_001.adsorption_model = 'LINEAR'
    Cyclic_model.root.input.model.unit_001.adsorption.is_kinetic = True    # Kinetic binding
    Cyclic_model.root.input.model.unit_001.adsorption.LIN_KA = [3.55]  # m^3 / (mol * s)   (mobile phase)
    Cyclic_model.root.input.model.unit_001.adsorption.LIN_KD = [0.1]      # 1 / s (desorption)
    # Initial conditions
    Cyclic_model.root.input.model.unit_001.init_c = [0]
    Cyclic_model.root.input.model.unit_001.init_q = [0]  # salt starts at max capacity

    # Grid cells in column and particle: the most important ones - ensure grid-independent solutions
    Cyclic_model.root.input.model.unit_001.discretization.SPATIAL_METHOD = "DG"
    Cyclic_model.root.input.model.unit_001.discretization.NELEM = nelem

    # Polynomial order
    Cyclic_model.root.input.model.unit_001.discretization.POLYDEG = polydeg
    Cyclic_model.root.input.model.unit_001.discretization.EXACT_INTEGRATION = exactInt

    # Bound states - for zero the compound does not bind, >1 = multiple binding sites
    Cyclic_model.root.input.model.unit_001.discretization.nbound = np.ones(n_comp, dtype=int)

    Cyclic_model.root.input.model.unit_001.discretization.PAR_DISC_TYPE = 'EQUIDISTANT_PAR'
    Cyclic_model.root.input.model.unit_001.discretization.USE_ANALYTIC_JACOBIAN = 1
    Cyclic_model.root.input.model.unit_001.discretization.RECONSTRUCTION = 'WENO'
    Cyclic_model.root.input.model.unit_001.discretization.GS_TYPE = 1
    Cyclic_model.root.input.model.unit_001.discretization.MAX_KRYLOV = 0
    Cyclic_model.root.input.model.unit_001.discretization.MAX_RESTARTS = 10
    Cyclic_model.root.input.model.unit_001.discretization.SCHUR_SAFETY = 1.0e-8

    Cyclic_model.root.input.model.unit_001.discretization.weno.BOUNDARY_MODEL = 0
    Cyclic_model.root.input.model.unit_001.discretization.weno.WENO_EPS = 1e-10
    Cyclic_model.root.input.model.unit_001.discretization.weno.WENO_ORDER = 3

    # Copy column models
    Cyclic_model.root.input.model.unit_002 = copy.deepcopy(Cyclic_model.root.input.model.unit_001)

    # Unit LRMP2 
    Cyclic_model.root.input.model.unit_002.adsorption.is_kinetic = False    # Kinetic binding
    Cyclic_model.root.input.model.unit_002.adsorption.LIN_KA = [35.5]  # m^3 / (mol * s)   (mobile phase)
    Cyclic_model.root.input.model.unit_002.adsorption.LIN_KD = [1]      # 1 / s (desorption)

    # To write out last output to check for steady state
    Cyclic_model.root.input['return'].WRITE_SOLUTION_LAST = True

    # % Input and connections
    # Sections
    Cyclic_model.root.input.solver.sections.nsec = 2
    Cyclic_model.root.input.solver.sections.section_times = [0, 100, 6000]

    # Feed and Eluent concentration
    Cyclic_model.root.input.model.unit_000.sec_000.const_coeff = [1]  # Inlet flowrate concentration

    Cyclic_model.root.input.model.unit_000.sec_001.const_coeff = [0]  # Inlet flowrate concentration

    # Connections
    Cyclic_model.root.input.model.connections.nswitches = 1

    Cyclic_model.root.input.model.connections.switch_000.section = 0
    Cyclic_model.root.input.model.connections.switch_000.connections = [
        0, 1, -1, -1, LRMP_Q3/2.0,  # flowrates, Q, m3/s
        1, 2, -1, -1, LRMP_Q3,
        2, 1, -1, -1, LRMP_Q3/2.0,
        2, 3, -1, -1, LRMP_Q3/2.0,
    ]

    # solution times
    Cyclic_model.root.input.solver.user_solution_times = np.linspace(0.0, 6000.0, 6001)

    # Time
    # Tolerances for the time integrator
    Cyclic_model.root.input.solver.time_integrator.ABSTOL = 1e-12  # absolute tolerance
    Cyclic_model.root.input.solver.time_integrator.ALGTOL = 1e-10
    Cyclic_model.root.input.solver.time_integrator.RELTOL = 1e-10  # Relative tolerance
    Cyclic_model.root.input.solver.time_integrator.INIT_STEP_SIZE = 1e-10
    Cyclic_model.root.input.solver.time_integrator.MAX_STEPS = 1000000

    # Solver options in general (not only for column although the same)
    Cyclic_model.root.input.model.solver.gs_type = 1
    Cyclic_model.root.input.model.solver.max_krylov = 0
    Cyclic_model.root.input.model.solver.max_restarts = 10
    Cyclic_model.root.input.model.solver.schur_safety = 1e-8
    
    # necessary specifically for this sim
    # Cyclic_model.root.input.solver.consistent_init_mode = 5
    Cyclic_model.root.input.solver.time_integrator.USE_MODIFIED_NEWTON = 1

    # Number of cores for parallel simulation
    Cyclic_model.root.input.solver.nthreads = 1

    # Specify which results we want to return
    # Return data
    Cyclic_model.root.input['return'].split_components_data = 0
    Cyclic_model.root.input['return'].split_ports_data = 0
    Cyclic_model.root.input['return'].unit_000.write_solution_bulk = 0
    Cyclic_model.root.input['return'].unit_000.write_solution_inlet = 0
    Cyclic_model.root.input['return'].unit_000.write_solution_outlet = 0
    Cyclic_model.root.input['return'].unit_001.write_solution_bulk = 0
    Cyclic_model.root.input['return'].unit_001.write_solution_inlet = 0
    Cyclic_model.root.input['return'].unit_001.write_solution_outlet = 1

    # Copy settings to the other unit operations
    Cyclic_model.root.input['return'].unit_002 = Cyclic_model.root.input['return'].unit_001
    Cyclic_model.root.input['return'].unit_003 = Cyclic_model.root.input['return'].unit_001

    
     #Saving data - do not touch
    Cyclic_model.filename = 'ModelDG.h5'
    Cyclic_model.save()

    #run model
    data = Cyclic_model.run()
    Cyclic_model.load()
    
    
    return Cyclic_model
    

#%%
#Import analytical solution from Jan
c_analytical = pd.read_csv('Semi-analytical_cyclic.csv')
c_analytical_1 = pd.read_csv('data.csv')



nCells = [1,2,4,8,16]
polyDeg = 2
runtime = []
MAE = []
DOF = []
polyDegu = []


for NElem in nCells:
    out = Cyclic_model1(NElem,polyDeg,1)
    err = max(abs(out.root.output.solution.unit_003.solution_outlet[1:,0] - c_analytical["C1"]))
    
    MAE.append(err)
    runtime.append(out.root.meta.time_sim)
    DOF.append(1 * NElem * (polyDeg + 1) * 3 * 2) # 3 phases, 2 columns
    polyDegu.append(polyDeg)


df = pd.DataFrame({'DOF': DOF, 'nCellu': nCells,'polyDegu': polyDegu,'runtime': runtime, 'MAE': MAE,})
df.to_csv('CADETDGConvergence.csv')

import os
os.remove("ModelDG.h5")
