#Run a convergence test where the results are compared to the semi-analytical solution.


# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 4)..., "include.jl"))

# Specify number of cells, polynomial degree and number of components 
nCell =  [1,2,4,8,16]
polyDeg = [2]
nComp = 1

# Load semi analytical solution
using CSV, DataFrames
c_analytical = CSV.read((joinpath(@__DIR__,"Semi-analytical_cyclic.csv")),DataFrame)

function model_setup(nCells, polyDeg, exactInt)

	nComp = 1
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
    model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [1]
    model["root"]["input"]["model"]["unit_000"]["sec_001"] = OrderedDict()
    model["root"]["input"]["model"]["unit_000"]["sec_001"]["const_coeff"] = [0] 
 


	
	# Set elements for outlet 
    model["root"]["input"]["model"]["unit_003"] = OrderedDict()
    model["root"]["input"]["model"]["unit_003"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_003"]["ncomp"] = nComp



	# Set elements sequentially for unit_001
	model["root"]["input"]["model"]["unit_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["unit_type"] = "LUMPED_RATE_MODEL_WITH_PORES"
	model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_001"]["col_porosity"] = 0.37
	model["root"]["input"]["model"]["unit_001"]["col_dispersion"] = 2e-7
	model["root"]["input"]["model"]["unit_001"]["col_length"] = 1.4e-2
	model["root"]["input"]["model"]["unit_001"]["cross_section_area"] = 1
	model["root"]["input"]["model"]["unit_001"]["par_porosity"] = 0.75
	model["root"]["input"]["model"]["unit_001"]["film_diffusion"] = [6.9e-6]
	model["root"]["input"]["model"]["unit_001"]["par_radius"] = 45e-6
    LRMP_Q3 = 3.45*1e-2 / 60 * 0.37
	model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "LINEAR"

	model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = true
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["LIN_KA"] = [3.55] 
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["LIN_KD"] = [0.1]


	model["root"]["input"]["model"]["unit_001"]["init_c"] = [0]
	model["root"]["input"]["model"]["unit_001"]["init_q"] = [0]

	model["root"]["input"]["model"]["unit_001"]["discretization"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDeg"] = polyDeg
	model["root"]["input"]["model"]["unit_001"]["discretization"]["ncol"] = nCells
	model["root"]["input"]["model"]["unit_001"]["discretization"]["exact_integration"] = exactInt
	model["root"]["input"]["model"]["unit_001"]["discretization"]["nbound"] = ones(Bool, nComp)

    
	# Copy to remaining units
    model["root"]["input"]["model"]["unit_002"] = deepcopy(model["root"]["input"]["model"]["unit_001"])


    # Unit LRMP2 
    model["root"]["input"]["model"]["unit_002"]["adsorption"]["is_kinetic"] = false
	model["root"]["input"]["model"]["unit_002"]["adsorption"]["LIN_KA"] = [35.5] 
	model["root"]["input"]["model"]["unit_002"]["adsorption"]["LIN_KD"] = [1]
    
    
 
    
    #solution times
    model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
    model["root"]["input"]["solver"]["sections"]["nsec"] = 2
    model["root"]["input"]["solver"]["sections"]["section_times"] = [0., 100, 6000]  


	# Set elements for connections
    model["root"]["input"]["model"]["connections"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["nswitches"] = 1

    model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
    model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] =[
        0, 1, -1, -1, LRMP_Q3/2,#flowrates, Q, m3/s
        1, 2, -1, -1, LRMP_Q3,
        2, 1, -1, -1, LRMP_Q3/2,
        2, 3, -1, -1, LRMP_Q3/2,
    ]

	# Set elements for user_solution_times
    model["root"]["input"]["solver"]["user_solution_times"] = collect(0: 1.0: 6000) 

    # Set elements for time_integrator
    model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
    model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
    model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
    model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10


	inlets, outlets, columns, switches, solverOptions = create_units(model)
	

	return inlets, outlets, columns, switches, solverOptions

	
end


MAE = []
DOF = []
rtime = []
polyDegu = []


for NElem in nCell #NElem = 4
    # set up the model 
    inlets, outlets, columns, switches, solverOptions = model_setup(NElem, polyDeg[1], true)

    # Solve model 
    rtime1 = @elapsed solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			alg = QNDF(autodiff=false), # Defaults to alg = QNDF(autodiff=false)
			)
    
    # Determine error 
    err = maximum(abs.(outlets[1].solution_outlet[2:end] .- c_analytical[:,1]))


    append!(MAE, err)
    append!(rtime, rtime1)
    append!(DOF,3*NElem*nComp*(polyDeg[1]+1) * solverOptions.nColumns)
    append!(polyDegu, polyDeg[1])

end

df = DataFrame(rtime = rtime,MAE = MAE, nCellu=nCell, DOF = DOF, polyDegu=polyDegu)
	
# Write results to CSV
CSV.write(joinpath(@__DIR__,"CADETJuliaConvergence.csv"),df)


