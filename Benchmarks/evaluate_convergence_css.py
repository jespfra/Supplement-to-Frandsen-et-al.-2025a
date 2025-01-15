
# A collection of functions used for CSS evaluations
import timeit
import numpy as np
import pandas as pd
import os

from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Inlet
from CADETProcess.processModel import FlowSheet, Process
from CADETProcess.simulator import Cadet as CadetSimulator

from cadet import Cadet
# Cadet.cadet_path = r'C:\Users\jespfra\Anaconda3\bin\cadet-cli'
# Cadet.cadet_path = r'C:\Users\Jespfra\AppData\Local\miniconda3\envs\CADETenv\bin\cadet-cli'
Cadet.cadet_path = r'C:\Users\pbzit\source\Test\out\install\aRELEASE\bin\cadet-cli'

os.path.dirname(os.path.dirname(Cadet.cadet_path))

#%%

def run_model_css(smb_model, model, max_iter = 1000, css_tol = 1e-4, c_analytical = []):
    # This runs the simulation till CSS takes place 
    # It evaluates using L2 norm whether CSS has taken place
    # If not, it loads the states and state derivatives from H5 file to rerun simulations 
    # It runs one complete cycle for each evaluation 
    
    # Number of components
    ncomp = smb_model.root.input.model.unit_004.ncomp
    
    # length of each cycle to evaluate CSS
    cycle_len = len(smb_model.root.input.solver.user_solution_times)
    
    
    # Load initial storage vectors 
    l2norm = []
    l2norm_analytical = []
    maxE = []
    rtime = []
    rtime_tot = []
    
    
    # Run initial simulation 
    data = smb_model.run()
    
    if data.returncode == 0:
        smb_model.load()
    else:
        print(data)
        print("Simulation failed")
        return 0, 0, 0, 0, 0, 0, 0, 0
    
    outletE = smb_model.root.output.solution.unit_002.solution_outlet
    outletR = smb_model.root.output.solution.unit_003.solution_outlet
    
    if smb_model.root.input.model.unit_004.discretization.spatial_method != b'DG':
        nelem = smb_model.root.input.model.unit_004.discretization.ncol
        npoly = 0
        npolypore = 0
        exactint = 0
        if smb_model.root.input.model.unit_004.unit_type == b'GENERAL_RATE_MODEL':
            nelem = np.log(nelem)/np.log(2)
    else:
        nelem = smb_model.root.input.model.unit_004.discretization.nelem
        npoly = smb_model.root.input.model.unit_004.discretization.polydeg
        npolypore = smb_model.root.input.model.unit_004.discretization.par_polydeg
        exactint = smb_model.root.input.model.unit_004.discretization.exact_integration
    
    
    
    # start while loop here 
    j = 0
    rtime0 = 0
    start = timeit.default_timer()
    while j < max_iter:
        
        # update state and state dot in smb_model
        # update state and state dot 
        state_y = smb_model.root.output.last_state_y
        state_ydot = smb_model.root.output.last_state_ydot
        
        # create new instance
        smb_model = model(nelem,npoly,npolypore,exactint,state_y,state_ydot)

        # Run simulation 
        data = smb_model.run()
        if data.returncode == 0:
            smb_model.load()
        else:
            print(data)
            print("Simulation failed")
            return int(j), l2norm, l2norm_analytical, maxE, 0, 0, outletE, outletR
        
        # append outlets 
        outletE = np.append(outletE, smb_model.root.output.solution.unit_002.solution_outlet, axis = 0)
        outletR = np.append(outletR, smb_model.root.output.solution.unit_003.solution_outlet, axis = 0)
        
        # Evaluate CSS
        l2norm_status, l2norm_analytical_status, maxE_status = eval_css(outletE, outletR, ncomp, cycle_len, c_analytical)
        
        # store data
        l2norm.append(l2norm_status)
        l2norm_analytical.append(l2norm_analytical_status)
        maxE.append(maxE_status)
        stop = timeit.default_timer()
        rtime.append(rtime0 + smb_model.root.meta.time_sim)
        rtime_tot.append(stop-start)
        # print(l2norm_status)
        if l2norm_status < css_tol:
            return int(j), l2norm, l2norm_analytical, maxE, rtime, rtime_tot, outletE, outletR
        
        
        # increment 
        rtime0 += smb_model.root.meta.time_sim
        j += 1
        
    
    # If CSS has not been achieved by now, not enought iterations has been make, throw error 
    print(f"Max iterations of {max_iter} has been simulated without achieving CSS. Try increasing the max_iter") 
    return int(j), l2norm, l2norm_analytical, maxE, rtime, rtime_tot, outletE, outletR
    
    
    
def eval_css(outletE, outletR, ncomp, cycle_len, c_analytical = []):
    # this function evaluates the norms and max error compared to an analytical css solution. 
    
    # Eval CSS using L2 norm 
    l2norm_status = 0
    l2norm_analytical_status = 0
    maxE_status = 0
    
    for k in range(ncomp):
        l2norm_status += np.linalg.norm(abs(outletE[-cycle_len:,k] - outletE[-2*cycle_len:-cycle_len,k]), ord=2)
        l2norm_status += np.linalg.norm(abs(outletR[-cycle_len:,k] - outletR[-2*cycle_len:-cycle_len,k]), ord=2)
        if type(c_analytical) == pd.core.frame.DataFrame:
            idxx_E = 'C' + str(k) + '_E' # Extract 
            idxx_R = 'C' + str(k) + '_R' # Raffinate
            l2norm_analytical_status += np.linalg.norm(abs(outletE[-cycle_len:,k] - c_analytical[-cycle_len:][idxx_E]), ord=2)
            l2norm_analytical_status += np.linalg.norm(abs(outletR[-cycle_len:,k] - c_analytical[-cycle_len:][idxx_R]), ord=2)
            maxE_status = max(maxE_status, max(abs(outletE[-cycle_len+2:,k] - c_analytical[-cycle_len+2:][idxx_E])))
            maxE_status = max(maxE_status, max(abs(outletR[-cycle_len+2:,k] - c_analytical[-cycle_len+2:][idxx_R])))
    return l2norm_status, l2norm_analytical_status, maxE_status
    
    
        
    
def Evaluate_CADETDG_css(model, c_analytical=[], nCells=[8], polyDeg=[4], polyDegPore=[0], nCellsPar=[1], max_iter=1000, css_tol=1e-4):
    # A function to benchmark evaluation of convergence to CSS.
	# It makes use of the either the coupled approach or using the OCA-FPI approach. 
    # the SMB setup should be model(nCells, polyDeg, polyDegPore) or model(nCells, polyDeg, exactInt)
    
    maxE_coupled_i = []
    L2_coupled_i = []
    L2_coupled_analytical_i = []
    runtime_iter_coupled_i = []
    runtime_iter_coupled_i_tot = []
    iter_coupled_i = []

    maxE_oca_fpi_i = []
    L2_oca_fpi_i = []
    L2_oca_fpi_analytical_i = []
    runtime_iter_oca_fpi_i = []
    runtime_iter_oca_fpi_i_tot = []
    iter_oca_fpi_i = []
	

    DOF = []
    nCellu = []
    polyDegu = []
    polyDegPoreu = []
    
    
    # Get number of components and columns from creating an instance of the model 
    smb_model = model(1,1,1)
    ncomp = smb_model.root.input.model.unit_000.ncomp
    # Find the number of columns 
    ncolumns = 0 
    i = 0
    while getattr(smb_model.root.input.model,f"unit_{i:03d}") != {}:
        unit = getattr(smb_model.root.input.model,f"unit_{i:03d}")
        if (unit.unit_type == 'LUMPED_RATE_MODEL_WITHOUT_PORES' or
            unit.unit_type == 'LUMPED_RATE_MODEL_WITH_PORES' or
            unit.unit_type == 'GENERAL_RATE_MODEL'):
            ncolumns += 1
    
            # Find the unit number which is a column
            if ncolumns == 1:
                unit_column_number = i
    
        # Iterate
        i += 1 
    
    transportModel = getattr(smb_model.root.input.model,f"unit_{unit_column_number:03d}").unit_type
    if smb_model.root.input.model.unit_004.discretization.SPATIAL_METHOD == 'DG':
        FV = False 
    else:
        FV = True

    if transportModel != "GENERAL_RATE_MODEL":
        iterPoly = polyDeg
    else:
        iterPoly = polyDegPore
    
    if FV:
        iterPoly = [0]
    # Run simulations
    for i in range(0, len(iterPoly)):
        for l in range(0, len(nCells)):
            print(f'Polynomial order {iterPoly[i]}')
            print(f'Column discretization {nCells[l]}')
            
            
            # Collocation, coupled
            print('Collocation, coupled')
            if transportModel == "GENERAL_RATE_MODEL":
                smb_model = model(nCells[l],polyDeg, polyDegPore[i],False)
            else:
                smb_model = model(nCells[l], polyDeg[i], 1, False)
            evals, l2norm, l2norm_analytical, maxE, rtime, rtime_tot, outletE, outletR = run_model_css(smb_model, model, max_iter, css_tol, c_analytical)
    
            # Storing data
            maxE_coupled_i.extend(maxE)
            L2_coupled_i.extend(l2norm)
            L2_coupled_analytical_i.extend(l2norm_analytical)
            iter_coupled_i.extend([evals] * len(maxE))
            runtime_iter_coupled_i.extend(rtime)
            runtime_iter_coupled_i_tot.extend(rtime_tot)



            # Collocation, OCA-FPI 
            # print('Collocation, OCA-FPI')
            # if transportModel == "GENERAL_RATE_MODEL":
            #     smb_model = model(nCells[l],polyDeg, polyDegPore[i],False)
            # else:
            #     smb_model = model(nCells[l], polyDeg[i], 1, False)
            # evals, l2norm, l2norm_analytical, maxE, rtime, outlet_E1, outlet_R1 = oca_fpi(smb_model, c_analytical, max_iter, css_tol)
            

            
            # Storing data
            # maxE_oca_fpi_i.extend(maxE)
            # L2_oca_fpi_i.extend(l2norm)
            # L2_oca_fpi_analytical_i.extend(l2norm_analytical)
            # iter_oca_fpi_i.extend([evals] * len(maxE))
            # runtime_iter_oca_fpi_i.extend(rtime)         
            
            
            nCellu.append(nCells[l])
            if FV:
                polyDegu.append(0)
            else:    
                polyDegu.append(smb_model.root.input.model.unit_004.discretization.polydeg)
                
            if transportModel == "GENERAL_RATE_MODEL":
                polyDegPoreu.append(iterPoly[i])
            else:
                polyDegPoreu.append(0)
            
            if FV:
                if transportModel == "LUMPED_RATE_MODEL_WITHOUT_PORES":
                    DOF.append( ncomp * nCells[l] * 2 * ncolumns)  # 2 phases
                elif transportModel == "LUMPED_RATE_MODEL_WITH_PORES":
                    DOF.append( ncomp * nCells[l] * 3 * ncolumns)   # 3 phases
                elif transportModel == "GENERAL_RATE_MODEL":
                    DOF.append( (ncomp * 2**nCells[l] + ncomp * 2 * 2**(nCells[l]-2) )* ncolumns)
            else:
                if transportModel == "LUMPED_RATE_MODEL_WITHOUT_PORES":                    
                    DOF.append(ncomp * nCells[l] * (polyDeg[i] + 1) * 2 * ncolumns)  # 2 phases
                elif transportModel == "LUMPED_RATE_MODEL_WITH_PORES":
                    DOF.append(ncomp * nCells[l] * (polyDeg[i] + 1) * 3 * ncolumns)  # 3 phases
                elif transportModel == "GENERAL_RATE_MODEL":
                    DOF.append(ncomp * nCells[l] * (polyDeg + 1) * ncolumns  + ncomp*2*nCells[l] * (polyDeg + 1)*(polyDegPore[i]+1) * nCellsPar * ncolumns)


            # Because the vectors have different lengths, they are adjusted with zeros to fill in to match the max length 
            max_length = max([len(maxE_coupled_i), len(runtime_iter_oca_fpi_i)]) + 1
            
            
            maxE_coupled_i.extend(np.zeros(abs(max_length - len(maxE_coupled_i))))
            L2_coupled_i.extend(np.zeros(abs(max_length - len(L2_coupled_i))))
            L2_coupled_analytical_i.extend(np.zeros(abs(max_length - len(L2_coupled_analytical_i))))
            runtime_iter_coupled_i.extend(np.zeros(abs(max_length - len(runtime_iter_coupled_i))))
            runtime_iter_coupled_i_tot.extend(np.zeros(abs(max_length - len(runtime_iter_coupled_i_tot))))
            iter_coupled_i.extend(np.zeros(abs(max_length - len(iter_coupled_i))))
            
            DOF.extend(np.ones(abs(max_length - len(DOF)))*DOF[-1])
            nCellu.extend(np.ones(abs(max_length - len(nCellu)))*nCellu[-1])
            polyDegu.extend(np.ones(abs(max_length - len(polyDegu)))*polyDegu[-1])
            polyDegPoreu.extend(np.ones(abs(max_length - len(polyDegPoreu)))*polyDegPoreu[-1])
            
            # maxE_oca_fpi_i.extend(np.zeros(abs(max_length - len(maxE_oca_fpi_i))))
            # L2_oca_fpi_i.extend(np.zeros(abs(max_length - len(L2_oca_fpi_i))))
            # L2_oca_fpi_analytical_i.extend(np.zeros(abs(max_length - len(L2_oca_fpi_analytical_i))))
            # runtime_iter_oca_fpi_i.extend(np.zeros(abs(max_length - len(runtime_iter_oca_fpi_i))))
            # iter_oca_fpi_i.extend(np.zeros(abs(max_length - len(iter_oca_fpi_i))))

            
            # Save results everytime a simulation as been carried out 
            convergenceDataDG = pd.DataFrame({  'DOF': DOF,
                                                'nCellu': nCellu,
                                                'polyDegu': polyDegu,
                                                'polyDegPoreu':polyDegPoreu,
                                                'maxE_coupled_i': maxE_coupled_i,
                                                'L2_coupled_i': L2_coupled_i,
                                                'L2_coupled_analytical_i': L2_coupled_analytical_i,
                                                'runtime_iter_coupled_i': runtime_iter_coupled_i,
                                                'runtime_iter_coupled_i_tot': runtime_iter_coupled_i_tot,
                                                'iter_coupled_i': iter_coupled_i,
                                                # 'maxE_oca_fpi_i': maxE_oca_fpi_i,
                                                # 'L2_oca_fpi_i': L2_oca_fpi_i,
                                                # 'L2_oca_fpi_analytical_i': L2_oca_fpi_analytical_i,
                                                # 'runtime_iter_oca_fpi_i': runtime_iter_oca_fpi_i,
                                                # 'iter_oca_fpi_i': iter_oca_fpi_i
                                                })
            #Save data in a CSV file
            # save results for in GSM results in case study folder
            if FV:
                convergenceDataDG.to_csv('CADETFVConvergenceCSS.csv')
            else:
                convergenceDataDG.to_csv('CADETDGConvergenceCSS.csv')
            
    
    
    

def oca_fpi(smb_model, c_analytical = [], max_iter = 1000, css_tol = 1e-4):
    # A function that evaluates CSS using the one-column analog with fixed point iteration
    # it takes the SMB model and solves the corresponding batch system
    
    # Number of components 
    ncomp = smb_model.root.input.model.unit_004.ncomp
    
    
    # Find the number of columns 
    ncolumns = 0 
    j = 0
    while getattr(smb_model.root.input.model,f"unit_{j:03d}") != {}:
        unit = getattr(smb_model.root.input.model,f"unit_{j:03d}")
        if (unit.unit_type == 'LUMPED_RATE_MODEL_WITHOUT_PORES' or
            unit.unit_type == 'LUMPED_RATE_MODEL_WITH_PORES' or
            unit.unit_type == 'GENERAL_RATE_MODEL'):
            ncolumns += 1
    
            # Find the unit number which is a column
            if ncolumns == 1:
                unit_column_number = j
    
        # Iterate
        j += 1 
    
    


    # Store the inlet flowrates and inlet concentrations
    # The inlets are from other columns and from feed/desorbent 
    # Cin are from feed/desorbent (constant)
    Cin = np.zeros([ncolumns,ncomp])
    ## Find ways to find these positions!!
    Cin[0,:] = smb_model.root.input.model.unit_000.sec_000.const_coeff                      # assuming you start with feed 
    Cin[int(ncolumns/2),:] = smb_model.root.input.model.unit_001.sec_000.const_coeff        # assuming desorbent is in the middle
    
    # Set up setups
    setups = [0]
    for j in range(1, ncolumns):
        setups.append(ncolumns - j)


    ## Flowrates into columns from inlets 
    Qin = np.zeros(ncolumns)
    
    connectionMatrix = np.array(smb_model.root.input.model.connections.switch_000.connections).reshape(-1, 5)
    feed_col = 4
    desorb_col = 4 + int(ncolumns/2)
    for row in connectionMatrix:
        if row[1] == feed_col and row[0] < feed_col:
            Qin[0] = row[4] # assuming you start with feed 
            break
    for row in connectionMatrix:
        if row[1] == desorb_col and row[0] < feed_col:
            Qin[int(ncolumns/2)] = row[4] # assuming desorbent is in the middle
            break
    
    # Flowrates into columns from other columns - assuming 4 column SMB written in specific format! 
    Qinrest = np.zeros(ncolumns)
    for j in range(ncolumns):
        connectionMatrixCol = np.array(getattr(smb_model.root.input.model.connections, f"switch_{setups[j]:03d}").connections).reshape(-1,5)
        for row in connectionMatrixCol:
            if row[1] == feed_col and row[0] > feed_col:
                Qinrest[j] = row[4]
    Qinrest = np.roll(Qinrest,1)
    outlet_R = np.empty(shape=(0,ncomp))
    outlet_E = np.empty(shape=(0,ncomp))
    # If we want to store the temporary concentrations 
    # For CSS, the number of cycles is unknown. 
  	# Hence, the length of the solution matrix is unknown and should be appended. 
  	# Therefore, the solution matrices are initially put with initial conditions. 
  	# Then, elements are appended using vcat matrix. 

      
    l2norm = []
    l2norm_analytical = []
    maxE = []
    rtime = []
    
    # Run initial simulation to get states
    j = 0
    times, sol_outlets = [0, smb_model.root.input.solver.sections.section_times[1]],np.zeros((2,ncomp))
    modelbatch = Flowrate_CADET_process(times, sol_outlets, ncomp, Cin[0,:])
    modelbatch = oca_fpi_batch_creator(modelbatch, smb_model, Qin[np.roll(setups,-j)[0]], Qinrest[np.roll(setups,-j)[0]], unit_column_number)       
    data = modelbatch.run()
    if data.returncode == 0:
        modelbatch.load()
    else:
        print(data)
        raise Exception("Simulation failed")
        
    # Update states
    state_y = modelbatch.root.output.last_state_y
    states_y = np.zeros([len(state_y), ncolumns])
    states_ydot = np.zeros([len(state_y), ncolumns])
    cycle_time = len(modelbatch.root.output.solution.solution_times)
    
    
    
    j = 0
    start = timeit.default_timer()
    while j < max_iter:
        # print(j)
        
        # Setup simulation
        # First get flowrates from piecewise polynomials
        modelbatch = Flowrate_CADET_process(times, sol_outlets, ncomp, Cin[np.roll(setups,-j)[0]])

        # setup simulation 
        modelbatch = oca_fpi_batch_creator(modelbatch, smb_model, Qin[np.roll(setups,-j)[0]], Qinrest[np.roll(setups,-j)[0]], unit_column_number, states_y[:,np.roll(setups,-j-1)[0]], states_ydot[:,np.roll(setups,-j-1)[0]])       


        # Run simulation
        data = modelbatch.run()
        if data.returncode == 0:
            modelbatch.load()
        else:
            print(data)
            raise Exception("Simulation failed")
            
        # Update states
        states_y[:,np.roll(setups,-j)[0]] = modelbatch.root.output.last_state_y
        states_ydot[:,np.roll(setups,-j)[0]] = modelbatch.root.output.last_state_ydot
        
        # section times and solution outlets for polynomial approximation
        times = modelbatch.root.output.solution.solution_times
        sol_outlets = modelbatch.root.output.solution.unit_002.solution_outlet
        
        

        # append outlets, if extract or raffinate, should be col_E and col_R
        if np.roll(setups,-j)[0] == 0:
            outlet_R = np.vstack((outlet_R, modelbatch.root.output.solution.unit_002.solution_outlet))
        elif np.roll(setups,-j)[0] == int(ncolumns/2):
            outlet_E = np.vstack((outlet_E, modelbatch.root.output.solution.unit_002.solution_outlet))
        
        # Evaluate CSS
        if j % ncolumns == 0 and j > (ncolumns + 1):
            
            l2norm_status, l2norm_analytical_status, maxE_status = eval_css(outlet_E, outlet_R, ncomp, cycle_time, c_analytical)
    
            # store data
            l2norm.append(l2norm_status)
            l2norm_analytical.append(l2norm_analytical_status)
            maxE.append(maxE_status)
            stop = timeit.default_timer()
            rtime.append(stop-start)
            # print(l2norm_status)
            if l2norm_status < css_tol:
                return int(j/ncolumns), l2norm, l2norm_analytical, maxE, rtime, outlet_E, outlet_R
        
        
        
        # increment 
        j += 1
    
    # If CSS has not been achieved by now, not enought iterations has been make, throw error 
    print(f"Max iterations of {max_iter} has been simulated without achieving CSS. Try increasing the max_iter") 
    return int(j/ncolumns), l2norm, l2norm_analytical, maxE, rtime, outlet_E, outlet_R
    
    


def Flowrate_CADET_process(times, sol_outlets, ncomp, Cin):
    
    component_system = ComponentSystem(ncomp)
    
    inlet_from_other_column = Inlet(component_system, 'inlet_from_other_column')
    
    inlet = Inlet(component_system, 'inlet')
    inlet.c = Cin
    
    # %
    
    flow_sheet = FlowSheet(component_system)
    
    flow_sheet.add_unit(inlet_from_other_column)
    flow_sheet.add_unit(inlet)
    
    # %
    process = Process(flow_sheet, 'inlet_profiles')
    process.cycle_time = times[-1]
    
    # %
    if len(times)>3:
        process.add_concentration_profile(inlet, times, sol_outlets, components=None,s=1e-6)
        process.add_concentration_profile(inlet_from_other_column, times, np.ones((len(times),ncomp))*Cin, components=None)


    # %
    simulator = CadetSimulator(install_path=os.path.dirname(os.path.dirname(Cadet.cadet_path)))
    cadet_process_model = simulator.get_process_config(process)
    
    return cadet_process_model
    

    



def oca_fpi_batch_creator(modelbatch, smb_model, Qin, Qinrest, unit_column_number, state_y = [0], state_ydot = [0]):
    # A function that takes the SMB setup and convertes it into a single column setup 
    # It gets the input from other columns as piecewise polynomials generated using CADET-Process 
    
    
    #Setting up the model
    model = Cadet()


    #Speciy number of unit operations: 2 input, column
    model.root.input.model.nunits = 3
    
    
    #First unit operation: inlet
    ## Inlet from Feed/desorbent
    model.root.input.model.unit_000 = modelbatch.input.model.unit_000
    
    
    ## Inlet from other columns - from modelBatch
    model.root.input.model.unit_001 = modelbatch.input.model.unit_001
 
    
    # Get column model 
    model.root.input.model.unit_002 = getattr(smb_model.root.input.model, f"unit_{unit_column_number:03d}")

    
    #Initial conditions
    if len(state_y) > 2:
        model.root.input.model.init_state_y = state_y
        model.root.input.model.init_state_ydot = state_ydot
        
        
    #Sections
    model.root.input.solver.sections.nsec = 1
    model.root.input.solver.sections.section_times = modelbatch.input.solver.sections.section_times   # s, section times
    model.root.input.solver.sections.section_continuity = [0]
     
    
    
    #Connections
    model.root.input.model.connections.nswitches = 1
    
    model.root.input.model.connections.switch_000.section = 0
    model.root.input.model.connections.switch_000.connections = [
        0, 2, -1, -1, Qin,                      # Feed/Desorbent, flowrates, Q, m3/s
        1, 2, -1, -1, Qinrest,                  # From other columns
    ]
    
    
    
    #solution times   
    filtered_times = [time for time in smb_model.root.input.solver.user_solution_times if time <= smb_model.root.input.solver.sections.section_times[1]]
    filtered_times[11:] = np.around(filtered_times[11:],6) # get rid of IDA 'too close to t0' error 
    model.root.input.solver.user_solution_times = filtered_times #np.unique(sorted(filtered_times + modelbatch.input.solver.sections.section_times))
    
    
    #To write out last output to check for steady state
    model.root.input['return'].WRITE_SOLUTION_LAST = True
    
    #Time 
    # Tolerances for the time integrator
    model.root.input.solver.time_integrator = smb_model.root.input.solver.time_integrator
    
    model.root.input.model.solver = smb_model.root.input.model.solver
    
    # Number of cores for parallel simulation
    model.root.input.solver.nthreads = 1
    
    #Specify which results we want to return
    # Return data
    model.root.input['return'].split_components_data = 0
    model.root.input['return'].split_ports_data = 0
    model.root.input['return'].unit_000.write_solution_bulk = 1
    model.root.input['return'].unit_000.write_solution_inlet = 1
    model.root.input['return'].unit_000.write_solution_outlet = 1
    
    
    # Copy settings to the other unit operations
    model.root.input['return'].unit_001 = model.root.input['return'].unit_000
    model.root.input['return'].unit_002 = model.root.input['return'].unit_000


        
    #Saving data - do not touch
    model.filename = 'Model.h5'
    model.save()
    
    return model


# def cube_interpolation(modelbatch, times, sol_outlets, Cin, splits = 20):
#     # Unit 000 corresponds to inlet from inlet
#     # Unit 001 corresponds to inlet from other columns 
    
#     # For the initial part
#     if len(times) < 2:
#         modelbatch.root.input.model.unit_000.sec_000.const_coeff = Cin
#         modelbatch.root.input.model.unit_001.sec_000.const_coeff = np.zeros(len(Cin))
#         return
    
#     ncomp = modelbatch.root.input.model.unit_000.ncomp
    
#     n_sections = len(times) //splits 
#     modelbatch.root.input.solver.sections.nsec = n_sections
    
#     section_times = [times[0]]
#     for i in range(n_sections):
#         poly_input = np.zeros([4, ncomp]) # From polynomial feed up to 3rd order
#         section_index = 'sec_{0:03d}'.format(i)
#         modelbatch.root.input.model.unit_000[section_index].const_coeff = Cin
        
#         for k in range(ncomp):
#             if i == range(n_sections)[-1]: # if at the end, take last split to interpolate
#                 remaining = len(times) % splits
#                 poly_input[:,k] = np.flip(np.polyfit(times[0 : splits + remaining], sol_outlets[-splits - remaining:, k], 3))
#             else:
#                 poly_input[:,k] = np.flip(np.polyfit(times[0 : splits], sol_outlets[i*splits : (i+1) * splits, k], 3))
                
#         modelbatch.root.input.model.unit_001[section_index].const_coeff = poly_input[0,:]    #Inlet flowrate concentration
#         modelbatch.root.input.model.unit_001[section_index].lin_coeff = poly_input[1,:]      #Inlet flowrate concentration
#         modelbatch.root.input.model.unit_001[section_index].quad_coeff = poly_input[2,:]     #Inlet flowrate concentration
#         modelbatch.root.input.model.unit_001[section_index].cube_coeff = poly_input[3,:]     #Inlet flowrate concentration
                
#         # Store time 
#         if i == range(n_sections)[-1]:
#             section_times.append(times[-1])
#         else:
#             section_times.append(times[(i+1) * splits - 1])
        
        
#     modelbatch.root.input.solver.sections.section_times = section_times

# def pp(poly_input,t):
#     return poly_input[0] + poly_input[1]*t + poly_input[2]*t**2 + poly_input[3]*t**3