




#  
function solve_model_css(; columns, switches::Switches, solverOptions, outlets=(0,), alg=QNDF(autodiff=false), c_analytical = nothing, col_E = 0, col_R = 0)
	"""
	Solve the differential equations using the ODE solver
	For this case, only the CSS is of interest of cyclic setups
	inputs
		- columns
		- switches
		- solverOptions
		- outlets
		- alg 
		- c_analytical 		- the analytical solution if available. Is used to compute the MAE and l2 norm 
		- col_E 			- for OCA-FPI, the column number for which the extract can be found 
		- col_R 			- for OCA-FPI, the column number for which the raffinate can be found '

	outputs 
		- iterations
		- l2norm 
		- l2norm_analytical - zeros if c_analytical is not provided 
		- maxE				- zeros if c_analytical is not provided 

	"""

	# preallocate return variables 
	l2norm = Float64[]
	l2norm_analytical = Float64[]
	maxE = Float64[]
	rtime = Float64[]

	# Align columns in the order of ncolumns
	# For a 4C SMB, setups = [1,4,3,2] i.e., counter SMB flow direction 
	# config = [1,2,..nColumns]
	setups = [1]
	config = [1]
	for i = 2:solverOptions.nColumns
		append!(setups, (solverOptions.nColumns - i + 2) )
		append!(config, i)
	end

	# If the outlet columns are specified for the SMB unit, overwrite the existing outlets. 
	# This can be necessary for the OCA-FPI for columns > 4. 
	if col_E != 0 && col_R != 0
		outlets[1].idx_unit[1] = col_E
		outlets[2].idx_unit[1] = col_R
	end
		
	# Define initial conditions 
	# if oca_fpi is activated, the initial conditions for all columns are split in columns
	if solverOptions.oca_fpi == true 
		x0 = zeros(columns[1].adsStride + columns[1].bindStride * columns[1].nComp * 2, solverOptions.nColumns)
		for i = 1:solverOptions.nColumns
			x0[:,i] = solverOptions.x0[1 + solverOptions.idx_units[i] : columns[i].adsStride + columns[i].bindStride * columns[i].nComp * 2 + solverOptions.idx_units[i]] 
		end

		# Define inlet as continuous interpolating polynomials
		# This is required to approximate the continuous input from one column to the other 
		# Should resemble the resulting chromatogram 
		# The initial one is trivial 
		poly_input = Array{Any}(undef, columns[1].nComp)
		for i = 1:columns[1].nComp
			poly_input[i] = LinearInterpolation(ones(length(solverOptions.solution_times)) .* x0[1], solverOptions.solution_times)
		end

		# For CSS, the number of cycles is unknown. 
		# Hence, the length of the solution matrix is unknown and should be appended. 
		# Therefore, the solution matrices are initially put with initial conditions. 
		# Then, elements are appended using vcat matrix. 
		for i in eachindex(columns)
			columns[i].solution_outlet = zeros(1, columns[1].nComp)
			for j = 1:columns[1].nComp
				columns[i].solution_outlet[1,j] = x0[j*columns[1].ConvDispOpInstance.nPoints, i]
			end
		end 
		for i in eachindex(outlets)
			outlets[i].solution_outlet = zeros(1, columns[1].nComp)
			for j = 1:columns[1].nComp
				outlets[i].solution_outlet[1,j] = x0[j*columns[1].ConvDispOpInstance.nPoints, i]
			end
		end 

		# The problem to be solved, either oca_fpi or normal 
		ode_prob = problem_oca_fpi!
	else 
		# For simulation the full SMB system, no changes should be made 
		x0 = solverOptions.x0
		

		# For CSS, the number of cycles is unknown. 
		# Hence, the length of the solution matrix is unknown and should be appended. 
		# Therefore, the solution matrices are initially put with initial conditions. 
		# Then, elements are appended using vcat matrix. 
		for i in eachindex(columns)
			columns[i].solution_outlet = hcat(columns[i].solution_outlet[1,:]')
		end 

		for i in eachindex(outlets)
			outlets[i].solution_outlet = hcat(outlets[i].solution_outlet[1,:]')
		end

		# The problem to be solved, either oca_fpi or normal 
		ode_prob = problem!
	end
	

	
	
	
	#running simulations
	i = 1
	start = now()
	while i < solverOptions.css_max_iter
		# Here different sections refer to different positions of the single column in the cyclic setup

		# Set up parameter vector and empty elements 
		jacProto = nothing
		p_jac = nothing
		analytical_jac = nothing
		
		if solverOptions.oca_fpi == true 
			p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, circshift(setups, 1-i)[1], 1, solverOptions.idx_units, switches, poly_input, p_jac) 
		else 
			p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, circshift(config, 1-i)[1], solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
		end


		# If Analytical Jacobian == yes, set analytical Jacobian
		# Is only supported for batch operation! 
		if solverOptions.analyticalJacobian == true
			if solverOptions.nColumns > 1 
				throw("Analytical Jacobian only supported for single column batch setups")
			end
			# determine static jacobian and allocation matrices that are stored in p_jac
			p_jac = jac_static(columns[1], switches.ConnectionInstance.u_tot[switches.switchSetup[i], 1], p) 
			p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
			analytical_jac = analytical_jac! #refers to the function
		end

		# If jacobian prototype, compute at every section time as switches might change Jacobian 
		if solverOptions.prototypeJacobian == true

			# determine jacobian prototype using finite differences - saves compilation time but perhaps change
			if solverOptions.oca_fpi == true
				jacProto = sparse(jac_finite_diff(ode_prob,p, x0[:,circshift(setups, 1-i)[1]] .+ 1e-6, 1e-8))
			else 
				jacProto = sparse(jac_finite_diff(ode_prob,p, x0 .+ 1e-6, 1e-8))
			end			

			# set dq0dq to zero when computing SMA w. formulation 1 as dq/dt is not needed to be computed
			# This makes it slightly faster
			# if typeof(bind)==SMA 
			# 	@. @views jacProto[1 +columns[h].adsStride +columns[h].nComp*units[h].bindStride :columns[h].bindStride +columns[h].adsStride +columns[h].nComp*units[h].bindStride] = 0
			# end
		end

		# update the tspan and the inlets through i to the system
		tspan = (switches.section_times[1], switches.section_times[2]) # assuming constant ts! 
		fun = ODEFunction(ode_prob; jac_prototype = jacProto, jac = analytical_jac)
		if solverOptions.oca_fpi == true
			prob = ODEProblem(fun, x0[:,circshift(setups, -i)[1]], (0, tspan[2]-tspan[1]), p)
		else 
			prob = ODEProblem(fun, x0, (0, tspan[2]-tspan[1]), p)
		end
		sol = solve(prob, alg, saveat=solverOptions.solution_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) 
		
		# For storing data
		vcatmatrix = ones(Float64, length(sol.t[2:end]), columns[1].nComp)
		
		if solverOptions.oca_fpi == true 

			# New initial conditions
			x0[:, circshift(setups, 1-i)[1]] = sol.u[end]

			# Solution should be stored depending on the column that is being simulated 
			# Hence the solution_outlet is written for i that corresponds to the column 
			# Extract solution in solution_outlet in each unit 

			for k = 1:columns[1].nComp 
				vcatmatrix[:,k] = sol(sol.t[2:end], idxs=k*columns[1].ConvDispOpInstance.nPoints + solverOptions.idx_units[1]).u
			end 
			columns[circshift(config, 1-i)[1]].solution_outlet = vcat(columns[circshift(setups, 1-i)[1]].solution_outlet, vcatmatrix)
			append!(columns[circshift(config, 1-i)[1]].solution_times, sol.t[2:end] .+ tspan[1])

			# Write outlets - if specified 	
			if outlets != (0,) 
				for j in eachindex(outlets)
					if outlets[j].idx_outlet != [-1] && outlets[j].idx_unit[1] == circshift(setups, 1-i)[1]
						for k = 1:columns[1].nComp
							
							vcatmatrix[:,k] = sol(sol.t[2:end], idxs=k*columns[1].ConvDispOpInstance.nPoints + solverOptions.idx_units[1]).u
						end 
						outlets[j].solution_outlet = vcat(outlets[j].solution_outlet, vcatmatrix)
						append!(outlets[j].solution_times,sol.t[2:end] .+ tspan[1])
					end
				end
			end

			# Write new polynomial functions #k=1
			for k = 1:columns[1].nComp
				poly_input[k] = CubicSpline(
				columns[circshift(config, 1-i)[1]].solution_outlet[end - length(sol.t) + 2: end, k],
				sol.t[2:end], 
				extrapolate = true
				)
			end

			# Check for CSS at outlet every ncolumns iterations after the first full cycle
			if i % solverOptions.nColumns != 0
				i += 1 
				continue
			end

		else 
			# New initial conditions
			x0 = sol.u[end]			

			# Store solution 
			#Extract solution in solution_outlet in each unit 
			for j = 1: solverOptions.nColumns
				for k = 1:columns[j].nComp 
					vcatmatrix[:,k] = sol(sol.t[2:end], idxs=k*columns[j].ConvDispOpInstance.nPoints + solverOptions.idx_units[j]).u
				end 
				columns[j].solution_outlet = vcat(columns[j].solution_outlet, vcatmatrix)
				append!(columns[j].solution_times, sol.t[2:end] .+ tspan[1])
			end

			# Write outlets - if specified 
			if outlets != (0,) 
				for j in eachindex(outlets)
					if outlets[j].idx_outlet != [-1]
						for k = 1:columns[1].nComp
							vcatmatrix[:,k] = sol(sol.t[2:end], idxs=k*columns[outlets[j].idx_unit[switches.switchSetup[circshift(config, 1-i)[1]]]].ConvDispOpInstance.nPoints + outlets[j].idx_outlet[switches.switchSetup[circshift(config, 1-i)[1]]]).u
						end 
						outlets[j].solution_outlet = vcat(outlets[j].solution_outlet, vcatmatrix)
						append!(outlets[j].solution_times,sol.t[2:end] .+ tspan[1])
					end
				end
			end
		end
		
		

		# Write to HDF5 using a function if relevant 
		

		# Check for CSS at outlets 
		if i > solverOptions.nColumns + 1
			css_status = 0.0
			l2norm_analytical_status = 0.0
			maxE_status = 0.0
			
			for j in eachindex(outlets) # j = 1
				
				for k = 1:columns[1].nComp # k = 1
					css_status += norm(abs.(outlets[j].solution_outlet[2 + end - length(sol.t) : end, k] - outlets[j].solution_outlet[3 + end - 2 * length(sol.t) : 1 + end - length(sol.t), k]), 2)

					# evaluate against analytical solution is provided 
					if c_analytical !== nothing
						if j == 1
							idxx = "C$(k-1)_E" # Extract 
						elseif j == 2
							idxx = "C$(k-1)_R" # Raffinate
						end
						l2norm_analytical_status += norm(abs.(outlets[j].solution_outlet[2 + end - length(sol.t) : end, k] - c_analytical[2 + end - length(sol.t) : end, idxx]), 2)
						maxE_status = maximum([maxE_status, maximum(abs.(outlets[j].solution_outlet[2 + end - length(sol.t) : end, k] - c_analytical[2 + end - length(sol.t) : end, idxx]))])
					end
				end
			end
			append!(l2norm_analytical, l2norm_analytical_status)
			append!(l2norm,css_status)
			append!(maxE,maxE_status)
			stop = now()
			append!(rtime,Dates.value(stop-start)/1000)
			# println(css_status)
			if css_status < solverOptions.css_tol 
				if solverOptions.oca_fpi == true
					return Int(i/solverOptions.nColumns), l2norm, l2norm_analytical, maxE, rtime
				else 
					return Int(i), l2norm, l2norm_analytical, maxE, rtime
				end
			end
		end

		# Increment
		i += 1
	end

	# If CSS has not been achieved by now, not enought iterations has been make, throw error 
	println("Max iterations of $(solverOptions.css_max_iter) has been simulated without achieving CSS. Try increasing the max_iter")
	return Int(i), l2norm, l2norm_analytical, maxE, rtime
end


# Define the function representing the differential equations for transport and binding
function problem_oca_fpi!(RHS, x, p, t)
    columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches, poly_input = p
	# i corresponds to section 
	
	@inbounds for h = 1:nColumns 
		# Compute binding term. 
		# The cpp, qq and rhs_q are set up to ease code reading
		cpp = @view x[1 + columns[h].adsStride + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h]]
		qq = @view x[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		RHS_q = @view RHS[1 + columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		compute_binding!(RHS_q, cpp, qq, columns[h].bind, columns[h].nComp, columns[h].bindStride, t)

		# Compute transport term
		compute_transport_oca_fpi!(RHS, RHS_q, cpp, x, columns[h], t, i, h, switches, idx_units, poly_input)

	end
	nothing
end