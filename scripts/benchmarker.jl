using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CSV
using DataFrames
using Test
using ProgressMeter
using LinearAlgebra
using Statistics
using JLD2
using BenchmarkTools



VAH = "VAH05_1501"


function f(VAH)
T = Float64
cache = CellFitElectrolyte.initialize_cache(Float64)
mass_mat = Matrix{Float64}(I,14,14)
mass_mat[14,14] = 0.0

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p::ComponentVector{T} = CellFitElectrolyte.p_transport()

split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("data/cycle_individual_data/$(VAH).csv",DataFrame)
df.times::Array{T,1} = df.times.-df.times[1]
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
current::Array{T,1} = -df.ImA./1000
p.Tamb = df.TemperatureC[1].+273.15
TData::Array{T,1} = df.TemperatureC[1:end-1]::Array{T,1}
VData::Array{T,1} = df.EcellV[1:end-1]::Array{T,1}


cycle_array::Array{T,1} = CellFitElectrolyte.current_profile(current,df.times)
num_steps::Int = Int(cycle_array[1])
times::Array{T,1} = cycle_array[2:num_steps+1]
types::Array{T,1} = cycle_array[num_steps+2:num_steps+1+num_steps]
values::Array{T,1} = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

cellgeometry::ComponentVector{T} = CellFitElectrolyte.cell_geometry()





    # Handle Initial Conditions
    u::Array{T,1}= Array{Float64,1}(undef,14)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    u[13] = p.Tamb
    du = similar(u)
    p.input_type = types[1]
    p.input_value = values[1]
    #Create Function and Initialize Integrator
    func = ODEFunction((du::Array{Float64,1},u::Array{Float64,1},p::ComponentVector{Float64},t::Float64)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv),mass_matrix=mass_mat)
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    #we're really only interested in temperature and voltage, so we'll just save those
    endV = Array{T,1}(undef,length(VData))
    endT = Array{T,1}(undef,length(VData))

function evaluator(p::ComponentVector{T},prob) where {T}
    input_type::T = types[1]
    input_value::T = values[1]
    prob = remake(prob,p=p)
    integrator = init(prob,Rosenbrock23(autodiff=false),save_everystep=false,tstops=df.times)



    for step::Int in 1:num_steps::Int-1
        input_type = types[step]::T
        input_value = values[step]::T
        end_time::T = times[step+1]::T
        dt::T = end_time-integrator.t
        p.input_value::T = input_value
        step!(integrator)
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
        endT[step] = integrator.u[13]  
    end

    V_rmse::T = @fastmath sqrt.(mean((endV.-VData).^2))
    T_rmse::T = @fastmath sqrt.(mean((endT.-TData.-273.15).^2))
    return V_rmse,T_rmse
end

@time evaluator(p,prob)
@time evaluator(p,prob)
@time evaluator(p,prob)
@time evaluator(p,prob)
V_RMSE,T_RMSE = evaluator(p,prob)
println("Voltage RMSE: $V_RMSE")
println("Temperature RMSE: $T_RMSE")

return 4

end



