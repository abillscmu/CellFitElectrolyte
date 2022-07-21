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

cache = CellFitElectrolyte.initialize_cache(Float64)
mass_mat = Matrix{Float64}(I,14,14)
mass_mat[14,14] = 0.0

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

VAH = "VAH05_1501"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("data/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
current = -df.ImA./1000
p.Tamb = df.TemperatureC[1].+273.15
TData = df.TemperatureC[1:end-1]
VData = df.EcellV[1:end-1]


cycle_array = CellFitElectrolyte.current_profile(current,df.times)
num_steps = Int(cycle_array[1])
times = cycle_array[2:num_steps+1]
types = cycle_array[num_steps+2:num_steps+1+num_steps]
values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

cellgeometry = CellFitElectrolyte.cell_geometry()








function evaluator(p::ComponentVector{T}) where {T}

    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,14)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    u[13] = p.Tamb
    du = similar(u)

    #Create Function and Initialize Integrator
    func = ODEFunction((du::Array{T,1},u::Array{T,1},p::ComponentVector{T},t::T)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv),mass_matrix=mass_mat)
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,Rodas4(autodiff=false),save_everystep=false)

    #we're really only interested in temperature and voltage, so we'll just save those
    endV::Array{T,1} = Array{T,1}(undef,length(df.EcellV)-1)
    endT::Array{T,1} = Array{T,1}(undef,length(df.EcellV)-1)

    for step::Int in 1:num_steps-1
        input_type::T = types[step]
        input_value::T = values[step]
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


using MATLABPlots
figure(1)
clf()
plot(df.times,df.EcellV)
hold_on()
plot(df.times[1:end-1],V)
xlabel("Time[s]")
ylabel("Voltage[V]")
legend(["Data","Model"])
title("VAH $cell cycle $cycle. Error: $(params.fval) mV")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))




