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
p = ComponentVector{Float64}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.7500861154034334, εₛ⁺ = 0.45039623350809316, δ⁻ = 3.815600768773315e-8, δ⁺ = 4.170570135429523e-8, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, E = 5000.0, ω = 0.01, n_li = 0.21)
VAH = ARGS[1]
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
filter!(row->row.Ns>=4,df)

vfull = 4.2
if cell==7
	vfull=4.0
elseif cell ==23
	vfull=4.1
end

initialcond = Dict("Starting Voltage[V]"=>vfull,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
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
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, verbose=false)

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

    V_rmse::T = mean(abs.((endV.-interpolated_voltage[1:end-1])))
    min_V_error = abs(minimum(endV) .- minimum(interpolated_voltage))
    V_T = (min_V_error + V_rmse)/2
    T_rmse::T = sqrt.(mean((endT.-interpolated_temperature[1:end-1]).^2))
    return V_T
end


params = [:εₛ⁻,:εₛ⁺,:δ⁻,:δ⁺, :n_li, :ω]
#params = filter(e->e∉[:Tamb,:c,:h,:x⁻₀,:β⁻,:β⁺,:βˢ,:Temp,:input_value,:input_type,:t⁺,:ω⁺,:E,:εₑˢ],params)
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



for (i,param) in enumerate(params)
    parent[i] = p[param]
    if param in [:εₛ⁻,:εₑ⁻,:εₑˢ,:εₑ⁺,:εₛ⁺]
        ub[i] = 1.0
        lb[i] = 0.0
    elseif param in [:R⁺,:θₛ⁺]
        ub[i] = p[param]*100.0
        lb[i] = 0.1*p[param]
    elseif param==:θₑ
        ub[i] = p[param]*10.0
        lb[i] = p[param]*0.001
    elseif param in [:δ⁺,:δ⁻]
        ub[i] = 1e-6
        lb[i] = 0.0
    elseif param in [:n_li]
	    ub[i] = 0.23
	    lb[i] = 0.12
    elseif param in [:ω]
        ub[i] = 0.05
        lb[i] = 0.0
    else
        ub[i] = parent[i].*10.0
        lb[i] = parent[i].*0.1
    end
end



loss(parent)
evaluator(p)

param = CellFitElectrolyte.anneal(loss,parent,ub,lb,options=options)

minimum_ps = param[1]
fval = param[2]


minimum_ps=vcat(minimum_ps,fval)
minimum_ps=vcat(minimum_ps,cell)
minimum_ps=vcat(minimum_ps,cycle)
newminimum_ps = [[x] for x in minimum_ps]


newparams = [p for p in params]
newparams = vcat(newparams,:fval)
newparams = vcat(newparams,:cell)
newparams = vcat(newparams,:cycle)

df = DataFrame(newminimum_ps,newparams)
CSV.write("$(VAH)_PARAM.csv",df)
sleep(5)


