using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CSV
using DataFrames
using Test
using ProgressMeter
using LinearAlgebra
using Statistics
using JLD2
using PythonPlot
using Turing

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

#Load Data
cell = 1
cycle = 2
vahnum = lpad(cell, 2, "0")
df = CSV.read("/Users/abills/Datasets/cycle_individual_data/VAH$(vahnum)_$cycle.csv", DataFrame)
df.times = df.times .- df.times[1]
filter!(row->row.Ns>=4,df)

#set up cycle arrays
cycle_array_vector = CellFitElectrolyte.load_airbus_cyclearrays()
cycle_array_vec = cycle_array_vector["cycle_array_vector"][cell]
cycle_array = cycle_array_vec[cycle]
cycle_array = CellFitElectrolyte.truncate_cycle_array(5, cycle_array)
num_steps = Int(cycle_array[1])
times = cycle_array[2:num_steps+1]
types = cycle_array[num_steps+2:num_steps+1+num_steps]
values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

function cycle_array_evaluator(p::ComponentVector{T}) where {T}
    #Prep Cycle Array
    num_steps = Int(cycle_array[1])
    times = cycle_array[2:num_steps+1]
    types = cycle_array[num_steps+2:num_steps+1+num_steps]
    values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

    #Initialize data
    tstops = df.times
    endV = Array{T, 1}(undef, length(df.times))

    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,8)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p = input_type,input_value
    Temp = df.TemperatureC[1] .+ 273.15
    @pack! p = Temp
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    du = similar(u)

    #Mass Matrix
    vars = ones(8)
    vars[8] = 0
    mm = diagm(vars)

    #Create Function and Initialize Integrator
    func = ODEFunction(( du, u, p, t)->CellFitElectrolyte.equations_electrolyte_life_allocating(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv), mass_matrix=mm)
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, tstops=tstops)
    default_dtmax = integrator.opts.dtmax
    integrator.opts.dtmax = 5.0
    integrator.opts.maxiters = 1e7

    #Iterate through each data point
    num_steps = Int(length(tstops))

    #Turn CCCV off (outside the loop)
    integrator.p.cccv_switch_2 = false
    integrator.p.cccv_switch = false


    if num_steps == -1.0
        error("lifetime is not supported by this simulator")
    else
        for step::Int in 1:num_steps-1
            #Calculate Voltage and Save
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            endV[step] = Voltage

            #Set Temperature
            Temp = df.TemperatureC[step+1].+273.15
            @pack! p = Temp

            #Get Input type
            end_time = tstops[step+1]
            start_time = tstops[step]
            cycle_array_idx = searchsortedlast(times, start_time)
            input_type = types[cycle_array_idx]
            input_value = values[cycle_array_idx]

            #Ensure integrator is going to find a happy place
            @pack! p = input_type,input_value
            if (input_type == 0.0) | (input_type == 5.0)
                u = deepcopy(integrator.u)
                u[8] = input_value
                OrdinaryDiffEq.set_u!(integrator, u)
                #integrator.opts.dtmax = 1.0
            elseif input_type == 1.0
                u = deepcopy(integrator.u)
                u[8] = input_value/Voltage
                OrdinaryDiffEq.set_u!(integrator, u)
            end

            #Simulate to the next timestop
            while integrator.t < end_time
                step!(integrator)
                savevalues!(integrator, true)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv, integrator.u[8])
                if ((Voltage >= p.vfull) & ((integrator.p.cccv_switch != true) & (input_type ==5)))
                    integrator.p.cccv_switch = true
                end
                if ((((integrator.u[8] >= integrator.p.ifull) & ((integrator.p.cccv_switch == true)) & (input_type ==5))) & (integrator.p.cccv_switch_2 != true))
                    u = deepcopy(integrator.u)
                    u[8] = 0
                    OrdinaryDiffEq.set_u!(integrator, u)
                    integrator.p.cccv_switch_2 = true
                end
            end
        end

        #Save final voltage
        Voltage = CellFitElectrolyte.calc_voltage(integrator.u, integrator.p, integrator.t, cache, cellgeometry, cathodeocv, anodeocv, integrator.u[8])
        endV[end] = Voltage
    end
    return endV
end

#Default Parameter Value
@model function fit_cfe(voltage)
    # Prior distributions.
    ω ~ truncated(Normal(0.02, 0.001), 0.00, 0.04)

    p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.75, εₛ⁺ = 0.6, εᵧ⁺ = 0, εᵧ⁻ = 0, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch_2 = false, cccv_switch = false, vfull = 4.2, ifull = -0.05)
    
    
    predicted = cycle_array_evaluator(p)

    # Observations.
    voltage ~ MvNormal(predicted, 0.1)

    return nothing
end

model = fit_cfe(df.EcellV)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, NUTS(0.65), 1000)






