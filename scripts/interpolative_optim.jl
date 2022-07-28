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

cache = CellFitElectrolyte.initialize_cache(Float64)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

VAH = ARGS[1]
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
current = -df.ImA./1000


#Get and set up interpolants
current_interpolant = LinearInterpolation(current,df.times)
voltage_interpolant = LinearInterpolation(df.EcellV,df.times)
temperature_interpolant = LinearInterpolation(df.TemperatureC.+273.15,df.times)

interpolated_time = collect(range(df.times[1],stop=df.times[end],step=1.0))


interpolated_current = current_interpolant.(interpolated_time)
interpolated_voltage = voltage_interpolant.(interpolated_time)
interpolated_temperature = temperature_interpolant.(interpolated_time)


#set up cycle arrays
cycle_array = CellFitElectrolyte.current_profile(interpolated_current,interpolated_time)
num_steps = Int(cycle_array[1])
times = cycle_array[2:num_steps+1]
types = cycle_array[num_steps+2:num_steps+1+num_steps]
values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]



cellgeometry = CellFitElectrolyte.cell_geometry()








function evaluator(p::ComponentVector{T}) where {T}

    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,7)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    Temp = df.TemperatureC[1].+273.15
    @pack! p = Temp
    du = similar(u)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p = input_type,input_value

    #Create Function and Initialize Integrator
    func = ODEFunction((du::Array{T,1},u::Array{T,1},p::ComponentVector{T},t::T)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,Rodas4(autodiff=false),save_everystep=false)

    #we're really only interested in temperature and voltage, so we'll just save those
    endV::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)
    endT::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)

    for step::Int in 1:num_steps-1
        Temp = interpolated_temperature[step]
        @pack! p = Temp
        input_type = types[step]
        input_value = values[step]
        end_time::T = times[step+1]
        @pack! p = input_type,input_value
        step!(integrator,end_time-integrator.t,true)
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        endT[step] = Temp
    end

    V_rmse::T = sqrt.(mean((endV.-interpolated_voltage[1:end-1]).^2))
    T_rmse::T = sqrt.(mean((endT.-interpolated_temperature[1:end-1]).^2))
    return endV
end


params = keys(p)
params = filter(e->e∉[:Tamb,:c,:h,:x⁻₀,:β⁻,:β⁺,:βˢ,:Temp,:input_value,:input_type],params)
parent = zeros(length(params))
ub = zeros(length(params))
lb = zeros(length(params))




function loss(vec)
    for (i,param) in enumerate(params)
        p[param] = vec[i]
    end
    p.input_type = 3
    p.input_value = 4.2
    l=10e6
    try
     l = evaluator(p)
    catch
        l=10e6
    end

    return l
end

options = CellFitElectrolyte.annealOptions( 1., 1000, 2000, 3000, 1e-4, -Inf, 2)

p.δ⁺ = 2.2e-8
p.δ⁻ = 2.2e-8
p.E = 5000

for (i,param) in enumerate(params)
    parent[i] = p[param]
    if param == :t⁺
        ub[i] = 0.75
        lb[i] = 0.25
    elseif param in [:εₛ⁻,:εₑ⁻,:εₑˢ,:εₑ⁺,:εₛ⁺]
        ub[i] = 1.0
        lb[i] = 0.0
    elseif param in [:R⁺,:θₛ⁺]
        ub[i] = p[param]*100.0
        lb[i] = 0.1*p[param]
    elseif param==:θₑ
        ub[i] = p[param]*10.0
        lb[i] = p[param]*0.001
    else
        ub[i] = p[param]*10.0
        lb[i] = p[param]*0.1
    end
end



loss(parent)

param = CellFitElectrolyte.anneal(loss,parent,ub,lb,options=options)

minimum = param[1]
fval = param[2]


minimum=vcat(minimum,fval)
minimum=vcat(minimum,cell)
minimum=vcat(minimum,cycle)
newminimum = [[x] for x in minimum]


newparams = [p for p in params]
newparams = vcat(newparams,:fval)
newparams = vcat(newparams,:cell)
newparams = vcat(newparams,:cycle)

df = DataFrame(newminimum,newparams)
CSV.write("$(VAH)_PARAM.csv",df)
sleep(5)


