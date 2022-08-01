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

VAH = "VAH02_500"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("data/cycle_individual_data/$(VAH).csv",DataFrame)
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
    return endV,endT
end


params = CSV.read("data/outputs/outputs0728_interpolative/$(VAH)_PARAM.csv",DataFrame)
param_sym = Symbol.(names(params))
for param in param_sym[1:end-3]
    if param in keys(p)
        p[param] = params[!,param][1]
    end
end

p.θₛ⁻ = 2e-8
p.ω⁻ = 0.0001
p.δ⁻ = 3.14e-9
p.δ⁺ = 0
p.εₛ⁻ = 0.6


V,t = evaluator(p)

V_rmse = sqrt.(mean((V.-interpolated_voltage[1:end-1]).^2))

println("V_rmse:",V_rmse)


using MATLABPlots
figure(1)
clf()
plot(df.times,df.EcellV)
hold_on()
plot(interpolated_time[1:end-1],V)
xlabel("Time[s]")
ylabel("Voltage[V]")
legend(["Data","Model"])
title("VAH $cell cycle $cycle. Error: $(params.fval) mV")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))



