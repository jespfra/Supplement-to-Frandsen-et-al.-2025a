#Run a convergence test where the results are compared to the semi-analytical solution.


# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 4)..., "include.jl"))
include(joinpath(@__DIR__, fill("..", 3)..., "evaluate_convergence.jl"))

# Specify number of cells, polynomial degree and number of components 
nCell =  [2,4,8,16]
polyDeg = [4,5,6]
nComp = 2

# Load semi analytical solution
using CSV, DataFrames
c_analytical = CSV.read((joinpath(@__DIR__,"Semi-analytical_LRM_Linear.csv")),DataFrame)

function model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts, nCells, polyDeg, exactInt, oca_fpi=false, css = false, analJac=false)

	nComp = 2
	model = OrderedDict(
		"root" => OrderedDict(
			"input" => OrderedDict(
				"model" => OrderedDict()
			)
		)
	)


	# Set elements sequentially for unit_000
    model["root"]["input"]["model"]["unit_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["unit_type"] = "INLET"
    model["root"]["input"]["model"]["unit_000"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_000"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

    model["root"]["input"]["model"]["unit_000"]["sec_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [2.78,2.78]


    model["root"]["input"]["model"]["unit_001"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["unit_type"] = "INLET"
    model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_001"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

    model["root"]["input"]["model"]["unit_001"]["sec_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["sec_000"]["const_coeff"] = [0,0]
	
	# Set elements for outlet 
    model["root"]["input"]["model"]["unit_002"] = OrderedDict()
    model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp

    # Set elements for outlet
    model["root"]["input"]["model"]["unit_003"] = OrderedDict()
    model["root"]["input"]["model"]["unit_003"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_003"]["ncomp"] = nComp


	# Set elements sequentially for unit_001
	model["root"]["input"]["model"]["unit_004"] = OrderedDict()
	model["root"]["input"]["model"]["unit_004"]["unit_type"] = "LUMPED_RATE_MODEL_WITHOUT_PORES"
	model["root"]["input"]["model"]["unit_004"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_004"]["col_porosity"] = 0.38
	model["root"]["input"]["model"]["unit_004"]["col_dispersion"] = 3.81e-6
	model["root"]["input"]["model"]["unit_004"]["col_length"] = 5.36e-1
	model["root"]["input"]["model"]["unit_004"]["cross_section_area"] = 5.31e-4
	model["root"]["input"]["model"]["unit_004"]["adsorption_model"] = "LINEAR"

	model["root"]["input"]["model"]["unit_004"]["adsorption"] = OrderedDict()
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["is_kinetic"] = false
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["LIN_KA"] = [0.54, 0.28] 
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["LIN_KD"] = [1.0, 1.0]

	model["root"]["input"]["model"]["unit_004"]["init_c"] = [0, 0]
	model["root"]["input"]["model"]["unit_004"]["init_q"] = [0, 0]

	model["root"]["input"]["model"]["unit_004"]["discretization"] = OrderedDict()
	model["root"]["input"]["model"]["unit_004"]["discretization"]["polyDeg"] = polyDeg
	model["root"]["input"]["model"]["unit_004"]["discretization"]["ncol"] = nCells
	model["root"]["input"]["model"]["unit_004"]["discretization"]["exact_integration"] = exactInt
	model["root"]["input"]["model"]["unit_004"]["discretization"]["nbound"] = ones(Bool, nComp)
	model["root"]["input"]["model"]["unit_004"]["discretization"]["use_analytic_jacobian"] = analJac
	
	# Copy to remaining units
    model["root"]["input"]["model"]["unit_005"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
    model["root"]["input"]["model"]["unit_006"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
    model["root"]["input"]["model"]["unit_007"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
	model["root"]["input"]["model"]["unit_008"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
	model["root"]["input"]["model"]["unit_009"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
	model["root"]["input"]["model"]["unit_010"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
	model["root"]["input"]["model"]["unit_011"] = deepcopy(model["root"]["input"]["model"]["unit_004"])



	# Set elements for solver
    n_cycles = 10
    switch_time = ts #s
    model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
    model["root"]["input"]["solver"]["sections"]["nsec"] = 8*n_cycles
    # Define switch times where cycle change
    section_times =  Float64[]
    push!(section_times, 0)
    for i in 0:n_cycles-1
        push!(section_times, Int64((8*i+1)*ts))
        push!(section_times, Int64((8*i+2)*ts))
        push!(section_times, Int64((8*i+3)*ts))
        push!(section_times, Int64((8*i+4)*ts))
		push!(section_times, Int64((8*i+5)*ts))
		push!(section_times, Int64((8*i+6)*ts))
		push!(section_times, Int64((8*i+7)*ts))
		push!(section_times, Int64((8*i+8)*ts))
    end

    model["root"]["input"]["solver"]["sections"]["section_times"] = section_times
    model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


	# Set elements for connections
    model["root"]["input"]["model"]["connections"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["nswitches"] = 8
	
    model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
    model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] =[
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

    model["root"]["input"]["model"]["connections"]["switch_001"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_001"]["section"] = 1
    model["root"]["input"]["model"]["connections"]["switch_001"]["connections"] =[
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

    model["root"]["input"]["model"]["connections"]["switch_002"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_002"]["section"] = 2
    model["root"]["input"]["model"]["connections"]["switch_002"]["connections"] =[
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

    model["root"]["input"]["model"]["connections"]["switch_003"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_003"]["section"] = 3
    model["root"]["input"]["model"]["connections"]["switch_003"]["connections"] =[
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
	
	model["root"]["input"]["model"]["connections"]["switch_004"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_004"]["section"] = 4
    model["root"]["input"]["model"]["connections"]["switch_004"]["connections"] =[
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
	
	model["root"]["input"]["model"]["connections"]["switch_005"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_005"]["section"] = 5
    model["root"]["input"]["model"]["connections"]["switch_005"]["connections"] =[
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
	
	model["root"]["input"]["model"]["connections"]["switch_006"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_006"]["section"] = 6
    model["root"]["input"]["model"]["connections"]["switch_006"]["connections"] =[
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
	
	model["root"]["input"]["model"]["connections"]["switch_007"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_007"]["section"] = 7
    model["root"]["input"]["model"]["connections"]["switch_007"]["connections"] =[
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
	
	


	# Set elements for user_solution_times
    model["root"]["input"]["solver"]["user_solution_times"] = collect(0: 1.0: n_cycles*8*ts) 

    # Set elements for time_integrator
    model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
    model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
    model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
    model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10

    # Enable One-column-analog using fixed-point-iteration (OCA-FPI) for SMB setups
    model["root"]["input"]["solver"]["css"] = OrderedDict() 
    model["root"]["input"]["solver"]["css"]["activate"] = css
    model["root"]["input"]["solver"]["css"]["oca_fpi"] = oca_fpi        # Activate OCA-FPI - should run tests for cyclic SMB setups 
    model["root"]["input"]["solver"]["css"]["css_max_iter"] = 1000      # Set max number of iterations 
    model["root"]["input"]["solver"]["css"]["css_tol"] = 1e-6           # Set cyclic steady state tolerance between two cycles  



	inlets, outlets, columns, switches, solverOptions = create_units(model)
	

	return inlets, outlets, columns, switches, solverOptions

	
end


# Determining flows 
ts = 1552
QD = 4.14e-8
QE = 3.48e-8
QF = 2.00e-8
QR = 2.66e-8
Q2 = 1.05e-7
Q3 = Q2 + QF
Q4 = Q3 - QR
Q1 = Q4 + QD

inlets, outlets, columns, switches, solverOptions = model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts, 4, 4, 1)

# Solve model 
solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			alg = QNDF(autodiff=false), # Defaults to alg = QNDF(autodiff=false)
			)

# Evaluate the convergence using the evaluate_convergence function 
evaluate_convergence(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,QNDF(autodiff=false), c_analytical, nComp, nCell, polyDeg, 1, "LRM", @__DIR__)



