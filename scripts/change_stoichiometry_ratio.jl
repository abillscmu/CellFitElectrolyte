using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using Turing
using DynamicHMC
using CSV
using DataFrames
using Test
using ProgressMeter
using LinearAlgebra
using Statistics
using JLD2
using PythonPlot
pygui(true)
using KernelDensity

cache = CellFitElectrolyte.initialize_cache(Float64)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

cyc = ENV["SLURM_ARRAY_TASK_ID"]

VAH = "VAH01_$(cyc)"
#split1 = split(VAH,['H','_'])
#cell = parse(Int,split1[2])
#cycle = parse(Int,split1[3])


df = CSV.read("~/cycle_individual_data/$(VAH).csv",DataFrame)
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
    func = ODEFunction((du, u, p, t)->CellFitElectrolyte.equations_electrolyte_allocating(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, tstops = times)

    #we're really only interested in temperature and voltage, so we'll just save those
    endV::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)
    endt::Array{T,1} = Array{T,1}(undef,length(interpolated_voltage)-1)

    for step::Int in 1:num_steps-1
        Temp = interpolated_temperature[step]
        @pack! p = Temp
        input_type = types[step]
        input_value = values[step]
        end_time::T = times[step+1]
        @pack! p = input_type,input_value
        while integrator.t < times[step+1]
            step!(integrator)
        end
        if any(integrator.u .< 0)
            endV[step] = integrator.u[findfirst(integrator.u .< 0)]
            continue
        end
        cₛˢ⁺ = u[7]
        cₛˢ⁻ = u[1]
        x⁺ = (cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min)
        x⁻ = calcocv(anodeocv,(cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),Temp)
        if (x⁺ >= 1)
            endV[step] = 50*x⁺
            continue
        elseif (x⁻ >= 1)
            endV[step] = 50*x⁻
            continue
        end
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        endt[step] = integrator.t
    end
    return endV
end




@model function fit_cfe(interpolated_voltage)
    # Prior distributions.
    #n_li ~ truncated(Normal(0.2, 0.01), 0.16, 0.22)
    ω ~ truncated(Normal(0.02, 0.001), 0.01, 0.05)
    x⁻₀ ~ Uniform(0.5, 0.6)
    
    #coordinate transforms to stay in a nice area
    εₑ⁺ ~ Uniform(0.05, 0.5)
    εₑ⁻ ~ Uniform(0.05, 0.5)

    frac_sol_am_pos ~ Uniform(0.5, 1.0)
    frac_sol_am_neg ~ Uniform(0.5, 1.0)

    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)
    
    
    predicted = evaluator(p)

    # Observations.
    interpolated_voltage[1:end-1] ~ MvNormal(predicted, 0.1)

    return nothing
end


model = fit_cfe(interpolated_voltage)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, NUTS(0.65), MCMCSerial(), 1000, 3; progress=true)
d = Dict("chain" => chain)


save("results_1215_stoichs/$VAH.jld2", d)
sleep(5)

#=
εₛ⁻ = kde(chain[:εₛ⁻].data[:,1])
εₛ⁺ = kde(chain[:εₛ⁺].data[:,1])
εᵧ⁺ = kde(chain[:εᵧ⁺].data[:,1])
εᵧ⁻ = kde(chain[:εᵧ⁻].data[:,1])
ω = kde(chain[:ω].data[:,1])
n_li = kde(chain[:n_li].data[:,1])



figure(1);
clf()
subplot(321)
plot(εₛ⁻.x, εₛ⁻.density)
xlabel("εₛ⁻")
ylabel("Density")
grid()
subplot(322)
plot(ω.x, ω.density)
ylabel("Density")
xlabel("SEI Resistance")
grid()
subplot(323)
plot(n_li.x,n_li.density)
ylabel("Density")
xlabel("Moles Li")
grid()
subplot(324)
plot(εₛ⁺.x,εₛ⁺.density)
ylabel("Density")
xlabel("εₛ⁺")
grid()
subplot(325)
plot(εᵧ⁺.x,εᵧ⁺.density)
ylabel("Density")
xlabel("εᵧ⁺")
grid()
subplot(326)
plot(εᵧ⁻.x,εᵧ⁻.density)
ylabel("Density")
xlabel("εᵧ⁻")
grid()
=#
