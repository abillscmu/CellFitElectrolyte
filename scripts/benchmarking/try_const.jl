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

const cache = CellFitElectrolyte.initialize_cache(Float64)
myvec = ones(14)
myvec[end] = 0.0
const mass_mat = diagm(myvec)


const cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

VAH = "VAH05_1501"
const split1 = split(VAH,['H','_'])
const cell = parse(Int,split1[2])
const cycle = parse(Int,split1[3])


const df = CSV.read("data/cycle_individual_data/$(VAH).csv",DataFrame)
const df.times = df.times.-df.times[1]
const initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
const current = -df.ImA./1000
p.Tamb = df.TemperatureC[1].+273.15
const TData = df.TemperatureC[1:end-1]
const VData = df.EcellV[1:end-1]


const cycle_array = CellFitElectrolyte.current_profile(current,df.times)
const num_steps = Int(cycle_array[1])
const times = cycle_array[2:num_steps+1]
const types = cycle_array[num_steps+2:num_steps+1+num_steps]
const values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

const cellgeometry = CellFitElectrolyte.cell_geometry()








function evaluator(p::ComponentVector{T}) where {T}

    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,14)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    u[13] = p.Tamb
    du = similar(u)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p = input_type,input_value

    #Create Function and Initialize Integrator
    func = ODEFunction((du::Array{T,1},u::Array{T,1},p::ComponentVector{T},t::T)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv),mass_matrix=mass_mat)
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,Rodas4(autodiff=false),save_everystep=false)

    #we're really only interested in temperature and voltage, so we'll just save those
    endV::Array{T,1} = Array{T,1}(undef,length(df.EcellV)-1)
    endT::Array{T,1} = Array{T,1}(undef,length(df.EcellV)-1)

    for step::Int in 1:num_steps-1
        input_type = types[step]
        input_value = values[step]
        if input_type==5
            input_type = 3
            input_value = values[step]
            @pack! p = input_type,input_value
            V::T = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
            while V<=initialcond["Starting Voltage[V]"]
                step!(integrator)
                V = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
            end
            input_type = 2
            input_value = initialcond["Starting Voltage[V]"]
            @pack! p = input_type,input_value
            while integrator.u[14]<=-0.1
                step!(integrator)
            end
            input_type = 0
            input_value = 0
            @pack! p = input_type,input_value
            while integrator.t<=times[step+1]
                step!(integrator)
            end
        else
            end_time::T = times[step+1]
            @pack! p = input_type,input_value
            step!(integrator,end_time-integrator.t,true)
        end
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
        endT[step] = integrator.u[13]  
    end

    V_rmse::T = sqrt.(mean((endV.-VData).^2))
    T_rmse::T = sqrt.(mean((endT.-TData.-273.15).^2))
    return endV,endT
end



params = CSV.read("../CellFitElectrolyteData/PARAM/$(VAH)_PARAM.csv",DataFrame)
param_sym = Symbol.(names(params))
for param in param_sym[1:end-3]
    p[param] = params[!,param][1]
end
V,t = evaluator(p)


bench = @benchmark evaluator($p)

display(bench)
