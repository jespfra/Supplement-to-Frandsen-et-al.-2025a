#Run a convergence test where the results are compared to the semi-analytical solution.


# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 5)..., "include.jl"))
include(joinpath(@__DIR__, fill("..", 4)..., "evaluate_convergence.jl"))
include(joinpath(@__DIR__, fill("..", 5)..., "src\\utils","solver_css.jl"))
using DataInterpolations

# Specify number of cells, polynomial degree and number of components 
nCell = [2,4]
polyDeg = [4]
polyDegPore = [4,6] #[4,5,6,8]
nComp = 4


# Load semi analytical solution
using CSV, DataFrames



function model_setup(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts,nCells, polyDeg, polyDegPore, exactInt, oca_fpi=false, css = false, analJac=false)


	nComp = 4
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
    model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [282,1,1,1]


    model["root"]["input"]["model"]["unit_001"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["unit_type"] = "INLET"
    model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_001"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

    model["root"]["input"]["model"]["unit_001"]["sec_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["sec_000"]["const_coeff"] = [282,0,0,0]
	
	# Set elements for outlet 
    model["root"]["input"]["model"]["unit_002"] = OrderedDict()
    model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp

    # Set elements for outlet
    model["root"]["input"]["model"]["unit_003"] = OrderedDict()
    model["root"]["input"]["model"]["unit_003"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_003"]["ncomp"] = nComp


	# Set elements sequentially for unit_004
	model["root"]["input"]["model"]["unit_004"] = OrderedDict()
	model["root"]["input"]["model"]["unit_004"]["unit_type"] = "GENERAL_RATE_MODEL"
	model["root"]["input"]["model"]["unit_004"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_004"]["col_porosity"] = 0.37
	model["root"]["input"]["model"]["unit_004"]["col_dispersion"] = 5.75e-8
	model["root"]["input"]["model"]["unit_004"]["col_length"] = 0.014
	model["root"]["input"]["model"]["unit_004"]["cross_section_area"] = 1
	model["root"]["input"]["model"]["unit_004"]["par_porosity"] = 0.75
	model["root"]["input"]["model"]["unit_004"]["par_radius"] = 4.5e-5
	model["root"]["input"]["model"]["unit_004"]["par_coreradius"] = 0
	model["root"]["input"]["model"]["unit_004"]["par_diffusion"] = [70e-11, 6.07e-11, 6.07e-11, 6.07e-11]
	model["root"]["input"]["model"]["unit_004"]["film_diffusion"] = [6.9e-6, 6.9e-6, 6.9e-6, 6.9e-6]
	model["root"]["input"]["model"]["unit_004"]["adsorption_model"] = "STERIC_MASS_ACTION"

	model["root"]["input"]["model"]["unit_004"]["adsorption"] = OrderedDict()
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["is_kinetic"] = true
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["SMA_KA"] = [0.0, 35.5e-3, 1.59e-3, 7.70e-3]
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["SMA_KD"] = [1.0, 1, 1, 1]
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["SMA_LAMBDA"] = 1200.0
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["SMA_NU"] = [0.0, 4.7, 5.29, 3.7 ]
	model["root"]["input"]["model"]["unit_004"]["adsorption"]["SMA_SIGMA"] = [0.0, 11.83, 10.6, 10.0]

	model["root"]["input"]["model"]["unit_004"]["init_c"] = [282, 0, 0, 0]
	model["root"]["input"]["model"]["unit_004"]["init_q"] = [1200.0, 0, 0, 0]

	model["root"]["input"]["model"]["unit_004"]["discretization"] = OrderedDict()
	model["root"]["input"]["model"]["unit_004"]["discretization"]["polyDeg"] = polyDeg
	model["root"]["input"]["model"]["unit_004"]["discretization"]["polyDegPore"] = polyDegPore
	model["root"]["input"]["model"]["unit_004"]["discretization"]["ncol"] = nCells
	model["root"]["input"]["model"]["unit_004"]["discretization"]["exact_integration"] = exactInt
	model["root"]["input"]["model"]["unit_004"]["discretization"]["nbound"] = ones(Bool, nComp)
	model["root"]["input"]["model"]["unit_004"]["discretization"]["use_analytic_jacobian"] = analJac

	# Copy to remaining units
    model["root"]["input"]["model"]["unit_005"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
    model["root"]["input"]["model"]["unit_006"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
    model["root"]["input"]["model"]["unit_007"] = deepcopy(model["root"]["input"]["model"]["unit_004"])


	# Set elements for solver
    n_cycles = 2
    switch_time = ts #s
    model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
    model["root"]["input"]["solver"]["sections"]["nsec"] = 4*n_cycles
    # Define switch times where cycle change
    section_times =  Float64[]
    push!(section_times, 0)
    for i in 0:n_cycles-1
        push!(section_times, Int64((4*i+1)*ts))
        push!(section_times, Int64((4*i+2)*ts))
        push!(section_times, Int64((4*i+3)*ts))
        push!(section_times, Int64((4*i+4)*ts))
    end

    model["root"]["input"]["solver"]["sections"]["section_times"] = section_times
    model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


	# Set elements for connections
    model["root"]["input"]["model"]["connections"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["nswitches"] = 4
    model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
    model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] =[
        4, 5, -1, -1, Q4,#flowrates, Q, m3/s
        5, 6, -1, -1, Q4,
        6, 7, -1, -1, Q2,
        7, 4, -1, -1, Q2,
        0, 4, -1, -1, QF,
        1, 6, -1, -1, QD,
        4, 3, -1, -1, QR,
        6, 2, -1, -1, QE
    ]

    model["root"]["input"]["model"]["connections"]["switch_001"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_001"]["section"] = 1
    model["root"]["input"]["model"]["connections"]["switch_001"]["connections"] =[
        4,	5,	-1,	-1,	Q2,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q4,
        6,	7,	-1,	-1,	Q4,
        7,	4,	-1,	-1,	Q2,
        0,	5,	-1,	-1,	QF,
        1,	7,	-1,	-1,	QD,
        5,	3,	-1,	-1,	QR,
        7,	2,	-1,	-1,	QE
    ]

    model["root"]["input"]["model"]["connections"]["switch_002"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_002"]["section"] = 2
    model["root"]["input"]["model"]["connections"]["switch_002"]["connections"] =[
        4,	5,	-1,	-1,	Q2,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q2,
        6,	7,	-1,	-1,	Q4,
        7,	4,	-1,	-1,	Q4,
        0,	6,	-1,	-1,	QF,
        1,	4,	-1,	-1,	QD,
        6,	3,	-1,	-1,	QR,
        4,	2,	-1,	-1,	QE
    ]

    model["root"]["input"]["model"]["connections"]["switch_003"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_003"]["section"] = 3
    model["root"]["input"]["model"]["connections"]["switch_003"]["connections"] =[
        4,	5,	-1,	-1,	Q4,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q2,
        6,	7,	-1,	-1,	Q2,
        7,	4,	-1,	-1,	Q4,
        0,	7,	-1,	-1,	QF,
        1,	5,	-1,	-1,	QD,
        7,	3,	-1,	-1,	QR,
        5,	2,	-1,	-1,	QE
    ]


	# Set elements for user_solution_times
    tt = [0.0]
    for i=1:4
        ttt = tt[end] .+ [1e-10, 5e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3] 
        append!(tt,vcat(ttt, collect(0.1 + (i-1)*ts: 0.1: i*ts)))
    end
    model["root"]["input"]["solver"]["user_solution_times"] = tt #vcat([0, 1e-10, 5e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3], LinRange(1e-1, 1*4*ts, 1*4*ts*10))


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

# Isotherm parameters
ka = [0.0, 35.5e-3, 1.59e-3, 7.70e-3]; # [-]
ionicCapacity = 1200.0; # [mol / m^3 PV]
v = [0.0, 4.7, 5.29, 3.7 ]; # [-], charge
Cs = 282

# Specified flowrates - based on linear approximation
m1 = ka[3]*(ionicCapacity/Cs)^v[3]*1.5 
m2 = 5.0 #
m3 = 10.0 # 
m4 = 0.7256 #
ColLength = 0.014
ts = 180 #s

eps_c = 0.37
eps_p = 0.75
eps_t = eps_c + (1-eps_c) * eps_p

Q1 = -(m1*eps_t - m1 - eps_t)*ColLength/ts
Q2 = -(m2*eps_t - m2 - eps_t)*ColLength/ts
Q3 = -(m3*eps_t - m3 - eps_t)*ColLength/ts
Q4 = -(m4*eps_t - m4 - eps_t)*ColLength/ts  
QD = -ColLength*(m1*eps_t - m4*eps_t - m1 + m4)/ts
QE = -ColLength*(m1*eps_t - m2*eps_t - m1 + m2)/ts
QF = ColLength*(m2*eps_t - m3*eps_t - m2 + m3)/ts
QR = -ColLength*(m3*eps_t - m4*eps_t - m3 + m4)/ts


# Evaluate CSS convergence 
c_analytical = CSV.read(joinpath(@__DIR__,"Semi-analytical_GRM_SMA_CSS.csv"), DataFrame)
evaluate_css(Q1,Q2,Q3,Q4,QF,QD,QE,QR,ts, QNDF(autodiff=false), c_analytical, nComp, nCell, polyDeg, polyDegPore, "GRM", @__DIR__)
