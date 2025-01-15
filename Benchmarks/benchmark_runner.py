

import pandas as pd
	
	
def runCadetDG(transportModel, c_analytical, polyDeg, nCells, ncolumns, polyDegPore=0,nCellsPar=1):
    
    DOF = []
    nCellu = []
    polyDegu = []
    nCellsParu = []
    polyDegPoreu = []
    maxE_e = []
    maxE_i = []
    runtime_e = []
    runtime_e_tot = []
    runtime_i = []
    runtime_i_tot = []

    
    if transportModel != "GRM":
        iterPoly = polyDeg
    else:
        iterPoly = polyDegPore
    
    # Run simulations
    for i in range(0, len(iterPoly)):
        for l in range(0, len(nCells)):
            print(f'Polynomial order {iterPoly[i]}')
            print(f'Column discretization {nCells[l]}')
    
            if transportModel != "GRM":
                df, runtime, runtime_tot = run_simulation(transportModel, nCells[l], polyDeg[i], 1)
            else:
                df, runtime, runtime_tot = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i], 1, nCellsPar)
            
            if runtime != 666:
                runtime_e.append(runtime)
                runtime_e_tot.append(runtime_tot)
                err = 0
                ncomp = int((df.shape[1]-1)/2)
                for k in range(ncomp): #Number of components
                    idxxE = f'C{k}_E'
                    idxxR = f'C{k}_R'
                    err = max([err,abs(df[idxxE][:] - c_analytical[idxxE][:]).max()])
                    err = max([err,abs(df[idxxR][:] - c_analytical[idxxR][:]).max()])
                maxE_e.append(err)
            else:
                runtime_e.append(0)
                runtime_e_tot.append(0)
                maxE_e.append(0)
    
    
    
            if transportModel != "GRM":
                df, runtime, runtime_tot = run_simulation(transportModel, nCells[l], polyDeg[i], 0)
            else:
                df, runtime, runtime_tot = run_simulation_GRM(transportModel, nCells[l], polyDeg, polyDegPore[i], 0, nCellsPar)
    
            if runtime != 666:
                runtime_i.append(runtime)
                runtime_i_tot.append(runtime_tot)
                err = 0
                for k in range(ncomp):  # Number of components
                    idxxE = f'C{k}_E'
                    idxxR = f'C{k}_R'
                    err = max([err,abs(df[idxxE][:] - c_analytical[idxxE][:]).max()])
                    err = max([err,abs(df[idxxR][:] - c_analytical[idxxR][:]).max()])
                maxE_i.append(err)
            else:
                runtime_i.append(0)
                runtime_i_tot.append(0)
                maxE_i.append(0)
            
            
            nCellu.append(nCells[l])
            polyDegu.append(iterPoly[i])
            polyDegPoreu.append(iterPoly[i])
            nCellsParu.append(nCellsPar)
            
            if transportModel == "LRM":
                DOF.append(ncomp * nCells[l] * (polyDeg[i] + 1) * 2 * ncolumns)  # 2 phases
            elif transportModel == "LRMP":
                DOF.append(ncomp * nCells[l] * (polyDeg[i] + 1) * 3 * ncolumns)  # 3 phases
            elif transportModel == "GRM":
                DOF.append(ncomp * nCells[l] * (polyDeg + 1) * ncolumns  + ncomp*2*nCells[l] * (polyDeg + 1)*(polyDegPore[i]+1) * nCellsPar * ncolumns)

            
            # Save results everytime a simulation as been carried out 
            if transportModel != "GRM":
                convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegu': polyDegu,'runtime_e': runtime_e, 'runtime_e_tot': runtime_e_tot,'maxE_e': maxE_e,'runtime_i': runtime_i, 'runtime_i_tot': runtime_i_tot,'maxE_i': maxE_i,})
            elif transportModel == "GRM":
                convergenceDataDG = pd.DataFrame({'DOF': DOF, 'nCellu': nCellu,'polyDegPoreu': polyDegPoreu, 'polyDegPoreu' : polyDegPoreu, 'nCellsParu' : nCellsParu,'runtime_e': runtime_e, 'runtime_e_tot': runtime_e_tot,'maxE_e': maxE_e,'runtime_i': runtime_i, 'runtime_i_tot': runtime_i_tot,'maxE_i': maxE_i,})
            #Save data in a CSV file
            # save results for in GSM results in case study folder
            convergenceDataDG.to_csv('CADETDGConvergence.csv')
                   
        
        
        
# Run simulation for LRM and LRMP
def run_simulation(transportModel, ncol, polydeg, is_exact):
    rtimes = [0,0,0]
    for i in range(1): # run 3 simulations 
        df, rtime, rtime_tot = model(ncol, polydeg, is_exact)
        rtimes[i] = rtime
        time.sleep(3)
        if len(df)<5: # If simulation crashed, store the rtime as 666
            rtimes[i] = 666
        
    return df, rtimes[0], rtime_tot #min(rtimes)
    
    
# Run simulation GRM
def run_simulation_GRM(transportModel, ncol, polydeg,polyDegPore, is_exact, nCellsPar=1):
    rtimes = [0,0,0]
    for i in range(1): # run 3 simulations 
        df, rtime, rtime_tot = model(ncol, polydeg,polyDegPore, is_exact)
        rtimes[i] = rtime
        time.sleep(3)
        if len(df)<5: # If simulation crashed, store the rtime as 666
            rtimes[i] = 666
        
    return df, rtimes[0], rtime_tot #min(rtimes)
	
