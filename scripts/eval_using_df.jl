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

VAH = "VAH02_200"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("data/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
idx = findfirst(isequal(0),df.Ns)
df = df[1:idx,:]


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




p = p = ComponentVector{Float64}(θₛ⁻ = 3.238105814128935e-8, θₑ = 1.6464068552786306, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.7500861154034334, εₛ⁺ = 0.45039623350809316, δ⁻ = 3.815600768773315e-8, δ⁺ = 4.170570135429523e-8, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.6, input_type = 3.0, input_value = 4.2, E = 5000.0, ω = 0.01)


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
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false)

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
    return endV,endT,integrator.sol.t
end

p = ComponentVector{Float64}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.7500861154034334, εₛ⁺ = 0.45039623350809316, δ⁻ = 3.815600768773315e-8, δ⁺ = 4.170570135429523e-8, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, E = 5000.0, ω = 0.0001, n_li = 0.16)
params = CSV.read("results/outputs1210_full/$(VAH)_PARAM.csv",DataFrame)
param_sym = Symbol.(names(params))
for param in param_sym[1:end-3]
    if param in keys(p)
        p[param] = params[!,param][1]
    end
end

#=
p.θₛ⁻ = 2e-8
p.ω⁻ = 0.0001
p.δ⁻ = 3.14e-9
p.δ⁺ = 0
p.εₛ⁻ = 0.6
=#

V,t = evaluator(p)

V_rmse = sqrt.(mean((V.-interpolated_voltage[1:end-1]).^2))

println("V_rmse:",V_rmse)


using PythonPlot
figure(1)
clf()
plot(interpolated_time[1:end-1],V,"--")
plot(df.times,df.EcellV)
PythonPlot.yticks(2.4:0.4:4.2)
legend(["Model", "Experiment"])
PythonPlot.grid()