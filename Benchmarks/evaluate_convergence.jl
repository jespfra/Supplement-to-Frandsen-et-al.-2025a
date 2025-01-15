
using Dates

function evaluate_convergence(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts, alg, c_analytical, nComp, nCell, polyDeg, polyDegPore::Union{Int64,Vector{Int64}}, transport_model, saveat)
	# A function to evaluate the convergence of the benchmarks. 

	# Check for correct typed transport model 
	println(transport_model)
	if transport_model != "LRM" && transport_model != "LRMP" && transport_model != "GRM"
		throw("Incorrect Transport model")
	end


	# Preload storage for data
	runtime_e = []
	maxE_e = []
	runtime_i = []
	maxE_i = []
	DOF = []

	nCellu = []
	polyDegu = []
	polyDegPoreu = []

		
	# Perform test runs and save data 
	if transport_model == "GRM"
		inlets1, outlets1, columns1, switches1, solverOptions1 = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], polyDegPore[1],1)
		solve_model(columns = columns1,switches = switches1,solverOptions = solverOptions1, outlets = outlets1, alg = alg)
		inlets2, outlets2, columns2, switches2, solverOptions2 = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[2],polyDeg[1], polyDegPore[1],0)
		solve_model(columns = columns2,switches = switches2,solverOptions = solverOptions2, outlets = outlets2, alg = alg)
	else 
		inlets1, outlets1, columns1, switches1, solverOptions1 = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1],1)
		solve_model(columns = columns1,switches = switches1,solverOptions = solverOptions1, outlets = outlets1, alg = alg) #alg = QNDF(autodiff=false)
		inlets2, outlets2, columns2, switches2, solverOptions2 = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[3],polyDeg[1],0)
		solve_model(columns = columns2,switches = switches2,solverOptions = solverOptions2, outlets = outlets2, alg = alg)
	end
	
	
	# Create a DataFrame with the details and output
	df1 = DataFrame(outlets1[1].solution_outlet,:auto)
	df1[!,"time"] = outlets1[1].solution_times
	df1[!, "nCell1"] .= nCell[1]
	df1[!, "polyDeg1"] .= polyDeg[1]
	df1[!, "polyDegPore1"] .= polyDegPore[1]
	for i=1:size(outlets2[1].solution_outlet)[2]
		label = "y$i"
		df1[!, label] = outlets2[1].solution_outlet[:,i]
		df1[!, "x$(i-1)_R"] = outlets1[2].solution_outlet[:,i] #raffinate
		df1[!, "y$(i-1)_R"] = outlets2[2].solution_outlet[:,i] #raffinate
	end
	df1[!, "nCell2"] .= nCell[2]
	df1[!, "polyDeg2"] .= polyDeg[1]
	df1[!, "polyDegPore2"] .= polyDegPore[1]

	# Write the DataFrame to a CSV file
	CSV.write((joinpath(saveat,"Profiles_data.csv")), df1)
	
	# Store three runtimes 
	rtime1 = zeros(Float64,2)

	for h = 1:length(polyDegPore)
		for i = 1:size(polyDeg)[1]
			for l=1:size(nCell)[1]
				println("polyDegPore = $(polyDegPore[h])")
				println("polyDeg = $(polyDeg[i])")
				println("nCell = $(nCell[l])")
			
				# Exact integration
				if transport_model == "GRM"
					inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 1) # must be predefined outside the loop to save modifications
					for m=1:2 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 1)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				else 
					inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 1) # must be predefined outside the loop to save modifications
					for m=1:2 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 1)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				end
				
				rtime = minimum(rtime1)
				err = 0.0
				for n = 0:size(outlets[1].solution_outlet)[2]-1
					err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$(n)_E"]))])
					err = maximum([err, maximum(abs.(outlets[2].solution_outlet[:,n+1]-c_analytical[:,"C$(n)_R"]))])
				end 
				
				#Storing data
				append!(runtime_e,rtime)
				append!(maxE_e,err[1])
				
				# Collocation 
				if transport_model == "GRM"
					for m=1:2 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 0)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				else 
					for m=1:2 # Run three times
						inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 0)
						rtime1[m] = @elapsed solve_model(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
					end
				end
				
				rtime = minimum(rtime1)
				err = 0.0
				for n = 0:size(outlets[1].solution_outlet)[2]-1
					err = maximum([err, maximum(abs.(outlets[1].solution_outlet[:,n+1]-c_analytical[:,"C$(n)_E"]))])
					err = maximum([err, maximum(abs.(outlets[2].solution_outlet[:,n+1]-c_analytical[:,"C$(n)_R"]))])
				end 
				
				# Storing data
				append!(runtime_i,rtime)
				append!(maxE_i,err[1])
				
				
				
				# Storing more data
				append!(nCellu,nCell[l])
				append!(polyDegu,polyDeg[i])
				append!(polyDegPoreu,polyDegPore[h])
				if transport_model == "LRM"
					append!(DOF,2*nCell[l]*nComp*(polyDeg[i]+1) * solverOptions.nColumns)
				elseif transport_model == "LRMP"
					append!(DOF,3*nCell[l]*nComp*(polyDeg[i]+1) * solverOptions.nColumns)
				elseif transport_model == "GRM"
					append!(DOF,nComp*(polyDeg[i]+1)*nCell[l] * solverOptions.nColumns + nComp*(polyDeg[i]+1)*nCell[l]*(polyDegPore[h]+1)*2 * solverOptions.nColumns) 
				end 
				
				
				df = DataFrame(runtime_e=runtime_e,maxE_e=maxE_e,runtime_i=runtime_i,maxE_i=maxE_i,DOF=DOF,nCellu=nCellu, polyDegu=polyDegu, polyDegPoreu=polyDegPoreu)
	
				# Write results to CSV
				CSV.write(joinpath(saveat,"CADETJuliaConvergence.csv"),df)
			
			end
		end 
	end


end

function evaluate_css(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts, alg, c_analytical, nComp, nCell, polyDeg, polyDegPore, transport_model, saveat, col_E = 0, col_R = 0)
	"""
	A function to evaluate convergence to CSS.
	It makes use of the solve_model_css which can solve the problem
	either as coupled or using the OCA-FPI approach. 
	Stores results as CSV file, located at the 'saveat' location 
	inputs
		- flowrates
		- ts 				- switching time
		- alg 				- algorithm to solve ODEs
		- c_analytical 		- the analytical solution if available. Is used to compute the MAE and l2 norm 
		- ncomp 			- number of components 
		- nCell 			- number of DG elements 
		- polydeg 			- polynomial degree on bulk phase 
		- polydegpore 		- polynomial degree in pore phase 
		- transport_model 	- transport model, either LRM, LRMP or GRM 
		- saveat 			- location to store results, typical @__DIR__ 
		- col_E 			- for OCA-FPI, the column number for which the extract can be found 
		- col_R 			- for OCA-FPI, the column number for which the raffinate can be found '

	outputs 
		- nothing

	"""
	
	

	# Check for correct typed transport model 
	println(transport_model)
	if transport_model != "LRM" && transport_model != "LRMP" && transport_model != "GRM"
		throw("Incorrect Transport model")
	end


	# Preload storage for data
	# storing maximum error, L2 norm with itself (stopping criterion), L2 with analytical solution, 
	# iterations, runtime total, runtime per iteration (average)
	# _e = exact DGSEM (dont run for now)
	# _i = collocation (inexact) DGSEM

	# maxE_coupled_e = []
	# L2_coupled_e = []
	# L2_coupled_analytical_e = []
	# runtime_tot_coupled_e = []
	# runtime_iter_coupled_e = []
	# iter_coupled_e = []
	
	maxE_coupled_i = []
	L2_coupled_i = []
	L2_coupled_analytical_i = []
	runtime_iter_coupled_i = []
	iter_coupled_i = []

	# maxE_oca_fpi_e = []
	# L2_oca_fpi_e = []
	# L2_oca_fpi_analytical_e = []
	# runtime_tot_oca_fpi_e = []
	# runtime_iter_oca_fpi_e = []
	# iter_oca_fpi_e = []

	maxE_oca_fpi_i = []
	L2_oca_fpi_i = []
	L2_oca_fpi_analytical_i = []
	runtime_iter_oca_fpi_i = []
	iter_oca_fpi_i = []
	

	DOF = Float64[]
	nCellu = Float64[]
	polyDegu = Float64[]
	polyDegPoreu = Float64[]

	
	
	# Perform test runs and save data 
	println("Running initial tests")
	if transport_model == "GRM"
		# exact DGSEM, coupled 
		# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], polyDegPore[1], 1, false, true) 
		# solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
		
		# inexact DGSEM, coupled
		inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], polyDegPore[1], 0, false, true)  
		solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
		
		# exact DGSEM, OCA FPI
		# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], polyDegPore[1], 1, true, true) 	 
		# solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical, col_E = col_E, col_R = col_R) 
		
		# inexact DGSEM, OCA FPI
		inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], polyDegPore[1], 0, true, true)	 
		solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical, col_E = col_E, col_R = col_R)
	else 
		# exact DGSEM, coupled 
		# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], 1, false, true) 
		# solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg)
		
		# inexact DGSEM, coupled
		inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], 0, false, true)  
		solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
		
		# exact DGSEM, OCA FPI
		# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], 1, true, true) 	 
		# solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical) 
		
		# inexact DGSEM, OCA FPI
		inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[1],polyDeg[1], 0, true, true)	 
		solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical, col_E = col_E, col_R = col_R)
	end
	
	
	
	for h = 1:length(polyDegPore)
		for i = 1:size(polyDeg)[1]
			for l=1:size(nCell)[1]
				println("polyDegPore = $(polyDegPore[h])")
				println("polyDeg = $(polyDeg[i])")
				println("nCell = $(nCell[l])")
			
				# # Exact integration, coupled 
				# if transport_model == "GRM"
					# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 1, false, true) 
					# rtime1 = @elapsed evals, l2norm, l2norm_anzalytical, maxE = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
				# else 
					# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 1, false, true) 
					# rtime1 = @elapsed evals, l2norm, l2norm_anzalytical, maxE = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
				# end
				
				
				# #Storing data
				# append!(maxE_coupled_e, maxE)
				# append!(L2_coupled_e, l2norm)								
				# append!(L2_coupled_analytical_e, l2norm_anzalytical)
				# append!(iter_coupled_e, ones(length(maxE))*evals)
				# append!(runtime_tot_coupled_e, ones(evals)*rtime1)
				# append!(runtime_iter_coupled_e, rtime1/evals * collect(1:evals))

				
				
				# Collocation, coupled
				if transport_model == "GRM"
					inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 0, false, true) 
					evals, l2norm, l2norm_anzalytical, maxE, rtime = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
				else 
					inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 0, false, true) 
					evals, l2norm, l2norm_anzalytical, maxE, rtime = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
				end
				
				
				#Storing data
				append!(maxE_coupled_i, maxE)
				append!(L2_coupled_i, l2norm)								
				append!(L2_coupled_analytical_i, l2norm_anzalytical)
				append!(iter_coupled_i, ones(length(maxE))*evals)
				append!(runtime_iter_coupled_i, rtime)
				
				
				
				
				# # Exact integration, OCA-FPI 
				# if transport_model == "GRM"
					# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 1, true, true) 
					# rtime1 = @elapsed evals, l2norm, l2norm_anzalytical, maxE = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
				# else 
					# inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 1, true, true) 
					# rtime1 = @elapsed evals, l2norm, l2norm_anzalytical, maxE = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical)
				# end
				
				
				# #Storing data
				# append!(maxE_oca_fpi_e, maxE)
				# append!(L2_oca_fpi_e, l2norm)								
				# append!(L2_oca_fpi_analytical_e, l2norm_anzalytical)
				# append!(iter_oca_fpi_e, ones(length(maxE))*evals)
				# append!(runtime_tot_oca_fpi_e, ones(evals)*rtime1)
				# append!(runtime_iter_oca_fpi_e, rtime1/evals * collect(1:evals))

				
				
				# Collocation, OCA-FPI 
				if transport_model == "GRM"
					inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], polyDegPore[h], 0, true, true) 
					evals, l2norm, l2norm_anzalytical, maxE, rtime = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical, col_E = col_E, col_R = col_R)
				else 
					inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCell[l],polyDeg[i], 0, true, true) 
					evals, l2norm, l2norm_anzalytical, maxE, rtime = solve_model_css(columns = columns,switches = switches,solverOptions = solverOptions, outlets = outlets, alg = alg, c_analytical = c_analytical, col_E = col_E, col_R = col_R)
				end
				
				
				#Storing data
				append!(maxE_oca_fpi_i, maxE)
				append!(L2_oca_fpi_i, l2norm)								
				append!(L2_oca_fpi_analytical_i, l2norm_anzalytical)
				append!(iter_oca_fpi_i, ones(length(maxE))*evals)
				append!(runtime_iter_oca_fpi_i, rtime)


				# Because the vectors have different lengths, they are adjusted with zeros to fill in to match the max length 
				max_length = maximum([length(runtime_iter_coupled_i), length(runtime_iter_oca_fpi_i)])
				
				
				# Storing more data
				append!(nCellu,nCell[l] * ones(max_length-length(nCellu)))
				append!(polyDegu,polyDeg[i] * ones(max_length-length(polyDegu)))
				append!(polyDegPoreu,polyDegPore[h] * ones(max_length-length(polyDegPoreu)))
				if transport_model == "LRM"
					append!(DOF,2*nCell[l]*nComp*(polyDeg[i]+1) * solverOptions.nColumns * ones(max_length-length(DOF)))
				elseif transport_model == "LRMP"
					append!(DOF,3*nCell[l]*nComp*(polyDeg[i]+1) * solverOptions.nColumns * ones(max_length-length(DOF)))
				elseif transport_model == "GRM"
					append!(DOF,(nComp*(polyDeg[i]+1)*nCell[l] * solverOptions.nColumns + nComp*(polyDeg[i]+1)*nCell[l]*(polyDegPore[h]+1)*2 * solverOptions.nColumns) * ones(max_length-length(DOF))) 
				end 
				
				
				# constructing the dict to fill in the zeros	
				# maxE_coupled_e = vcat(maxE_coupled_e, zeros(Int, abs(max_length - length(maxE_coupled_e))))
				# L2_coupled_e = vcat(L2_coupled_e, zeros(Int, abs(max_length - length(L2_coupled_e))))
				# L2_coupled_analytical_e = vcat(L2_coupled_analytical_e, zeros(Int, abs(max_length - length(L2_coupled_analytical_e))))
				# runtime_tot_coupled_e = vcat(runtime_tot_coupled_e, zeros(Int, abs(max_length - length(runtime_tot_coupled_e))))
				# runtime_iter_coupled_e = vcat(runtime_iter_coupled_e, zeros(Int, abs(max_length - length(runtime_iter_coupled_e))))
				# iter_coupled_e = vcat(iter_coupled_e, zeros(Int, abs(max_length - length(iter_coupled_e))))
				
				maxE_coupled_i = vcat(maxE_coupled_i, zeros(Int, abs(max_length - length(maxE_coupled_i))))
				L2_coupled_i = vcat(L2_coupled_i, zeros(Int, abs(max_length - length(L2_coupled_i))))
				L2_coupled_analytical_i = vcat(L2_coupled_analytical_i, zeros(Int, abs(max_length - length(L2_coupled_analytical_i))))
				runtime_iter_coupled_i = vcat(runtime_iter_coupled_i, zeros(Int, abs(max_length - length(runtime_iter_coupled_i))))
				iter_coupled_i = vcat(iter_coupled_i, zeros(Int, abs(max_length - length(iter_coupled_i))))

				# maxE_oca_fpi_e = vcat(maxE_oca_fpi_e, zeros(Int, abs(max_length - length(maxE_oca_fpi_e))))
				# L2_oca_fpi_e = vcat(L2_oca_fpi_e, zeros(Int, abs(max_length - length(L2_oca_fpi_e))))
				# L2_oca_fpi_analytical_e = vcat(L2_oca_fpi_analytical_e, zeros(Int, abs(max_length - length(L2_oca_fpi_analytical_e))))
				# runtime_tot_oca_fpi_e = vcat(runtime_tot_oca_fpi_e, zeros(Int, abs(max_length - length(runtime_tot_oca_fpi_e))))
				# runtime_iter_oca_fpi_e = vcat(runtime_iter_oca_fpi_e, zeros(Int, abs(max_length - length(runtime_iter_oca_fpi_e))))
				# iter_oca_fpi_e = vcat(iter_oca_fpi_e, zeros(Int, abs(max_length - length(iter_oca_fpi_e))))


				maxE_oca_fpi_i = vcat(maxE_oca_fpi_i, zeros(Int, abs(max_length - length(maxE_oca_fpi_i))))
				L2_oca_fpi_i = vcat(L2_oca_fpi_i, zeros(Int, abs(max_length - length(L2_oca_fpi_i))))
				L2_oca_fpi_analytical_i = vcat(L2_oca_fpi_analytical_i, zeros(Int, abs(max_length - length(L2_oca_fpi_analytical_i))))
				runtime_iter_oca_fpi_i = vcat(runtime_iter_oca_fpi_i, zeros(Int, abs(max_length - length(runtime_iter_oca_fpi_i))))
				iter_oca_fpi_i = vcat(iter_oca_fpi_i, zeros(Int, abs(max_length - length(iter_oca_fpi_i))))





				# write into dataframe
				df = DataFrame(DOF=DOF, nCellu=nCellu, polyDegu=polyDegu, polyDegPoreu=polyDegPoreu, 
								# maxE_coupled_e=maxE_coupled_e, L2_coupled_e=L2_coupled_e, L2_coupled_analytical_e=L2_coupled_analytical_e, runtime_tot_coupled_e=runtime_tot_coupled_e, runtime_iter_coupled_e=runtime_iter_coupled_e, iter_coupled_e=iter_coupled_e, 
								maxE_coupled_i=maxE_coupled_i, L2_coupled_i=L2_coupled_i, L2_coupled_analytical_i=L2_coupled_analytical_i, runtime_iter_coupled_i=runtime_iter_coupled_i, iter_coupled_i=iter_coupled_i,
								# maxE_oca_fpi_e=maxE_oca_fpi_e, L2_oca_fpi_e=L2_oca_fpi_e, L2_oca_fpi_analytical_e=L2_oca_fpi_analytical_e, runtime_tot_oca_fpi_e=runtime_tot_oca_fpi_e, runtime_iter_oca_fpi_e=runtime_iter_oca_fpi_e, iter_oca_fpi_e=iter_oca_fpi_e,
								maxE_oca_fpi_i=maxE_oca_fpi_i, L2_oca_fpi_i=L2_oca_fpi_i, L2_oca_fpi_analytical_i=L2_oca_fpi_analytical_i, runtime_iter_oca_fpi_i=runtime_iter_oca_fpi_i, iter_oca_fpi_i=iter_oca_fpi_i							
								)
				
				# Write results to CSV
				CSV.write(joinpath(saveat,"CADETJuliaConvergenceCSS.csv"),df)
				
			
			end
		end 
	end
	
end
	
	
	
	
	
