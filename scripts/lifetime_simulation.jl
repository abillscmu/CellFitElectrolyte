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

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
cellgeometry = CellFitElectrolyte.cell_geometry()

cell = 2

#set up cycle arrays
println("$(pwd())")
@load "src/cycle_array_vector_new.jld2"
cycle_array_vec = cycle_array_vector[cell][2:2]


function lifetime_evaluator(p::ComponentVector{T}) where {T}
    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,8)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    u[8] = 1.0
    Temp = 320.0
    @pack! p = Temp
    du = similar(u)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p = input_type,input_value

    #Create Function and Initialize Integrator
    func = ODEFunction((du::Array{T,1},u::Array{T,1},p::ComponentVector{T},t::T)->CellFitElectrolyte.equations_electrolyte_life(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=true)

    #we're really only interested in temperature and voltage, so we'll just save those
    endV::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)
    endT::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)


    for (i, cycle_array) in enumerate(cycle_array_vec)
        println(i)
        num_steps = Int(cycle_array[1])
        if num_steps == -1.0
            error("not yet")
        else
            times = cycle_array[2:num_steps+1]
            types = cycle_array[num_steps+2:num_steps+1+num_steps]
            values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
        end
        for step::Int in 1:num_steps-1
            println(step)
            @pack! p = Temp
            input_type = types[step]
            input_value = values[step]
            end_time::T = times[step+1]
            @pack! p = input_type,input_value
            try
                step!(integrator,end_time-integrator.t,true)
            catch
                return integrator
            end
            endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
            endT[step] = Temp
        end
    end

    V_rmse::T = sqrt.(mean((endV.-interpolated_voltage[1:end-1]).^2))
    T_rmse::T = sqrt.(mean((endT.-interpolated_temperature[1:end-1]).^2))
    return endV,endT,integrator.sol.t
end

p = ComponentVector{Float64}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.7500861154034334, εₛ⁺ = 0.45039623350809316, δ⁻ = 3.815600768773315e-8, δ⁺ = 4.170570135429523e-8, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, E = 5000.0, ω = 0.01, n_li = 0.21)
params = CSV.read("results/outputs1211_vah02_2/outputs1211_vah02/VAH02_10_PARAM.csv",DataFrame)
param_sym = Symbol.(names(params))

for param in param_sym[1:end-3]
    if param in keys(p)
        p[param] = params[!,param][1]
    end
end
p.ω = 0

integrator = lifetime_evaluator(p)