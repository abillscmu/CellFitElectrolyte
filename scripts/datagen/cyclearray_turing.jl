using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.Turing
using CellFitElectrolyte.DynamicHMC
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
cellgeometry = CellFitElectrolyte.cell_geometry()

#Parse Arguments
VAH = "VAH01_200"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])

#Load Data
df = CSV.read("/Users/abills/Datasets/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
filter!(row->row.Ns>=4,df)

#Set Full Voltage (100%SOC)
vfull = 4.2
if cell==7
	vfull=4.0
elseif cell ==23
	vfull=4.1
end
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
vfull = initialcond["Starting Voltage[V]"]

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
    u[8] = 10.0

    #Mass Matrix
    vars = ones(8)
    vars[8] = 0
    mm = diagm(vars)

    #Create Function and Initialize Integrator
    func = ODEFunction((du, u, p, t)->CellFitElectrolyte.equations_electrolyte_life_allocating(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv), mass_matrix=mm)
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = try
        integrator = init(prob,QNDF(autodiff=false),save_everystep=false, tstops=tstops, verbose=false)
    catch
        endV = 100 .*ones(length(endV))
        return endV
    end
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
            println(step)
            #Calculate Voltage and Save
            cₛˢ⁺ = integrator.u[7]
            cₛˢ⁻ = integrator.u[1]
            x⁺ = (cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min)
            x⁻ = (cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min)
            if (x⁺ >= 1)
                endV[step] = 50*x⁺
                continue
            elseif (x⁻ >= 1)
                endV[step] = 50*x⁻
                continue
            end
            Voltage = try
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                endV[step] = Voltage
            catch
                endV[step] = -50
                continue
            end

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
            #try
                while integrator.t < end_time
                    step!(integrator)
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
            #catch
            #        continue
            #end
        end

        #Save final voltage
        try
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u, integrator.p, integrator.t, cache, cellgeometry, cathodeocv, anodeocv, integrator.u[8])
            endV[end] = Voltage
        catch
            endV[end] = -50
        end
    end
    return endV
end

#Default Parameter Value
@model function fit_cfe(voltage)
    # Prior distributions.
    ω ~ truncated(Normal(0.02, 0.001), 0.0, 0.035)
    x⁻₀ =  0.6
     
     #coordinate transforms to stay in a nice area
     εₛ⁺ ~ Uniform(0.4, 1.0)
     εₛ⁻ ~ Uniform(0.4, 1.0)
 
     frac_rest_elec_pos ~ Uniform(0.0, 1.0)
     frac_rest_elec_neg ~ Uniform(0.0, 1.0)
 
     εₑ⁺ = (1 - εₛ⁺)*frac_rest_elec_pos
     εₑ⁻ = (1 - εₛ⁻)*frac_rest_elec_neg
 
     εᵧ⁺ = 1 - εₑ⁺ - εₛ⁺
     εᵧ⁻ = 1 - εₑ⁻ - εₛ⁻

    p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch_2 = false, cccv_switch = false, vfull = 4.2, ifull = -0.05)
    
    
    predicted = cycle_array_evaluator(p)

    # Observations.
    voltage ~ MvNormal(predicted, 0.1)

    return nothing
end

model = fit_cfe(df.EcellV)

p = ComponentVector{Float64}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.7500861154034334, εₛ⁺ = 0.45039623350809316, δ⁻ = 3.815600768773315e-8, δ⁺ = 4.170570135429523e-8, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, n_li = 0.21, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, εᵧ⁺=0, εᵧ⁻=0, E_SEI = 1000, cccv_switch=false, cccv_switch_2 = false, vfull = 4.2, ifull = -0.01)
params = CSV.read("results/outputs1211_full/$(VAH)_PARAM.csv",DataFrame)
param_sym = Symbol.(names(params))

for param in param_sym[1:end-3]
    if param in keys(p)
        p[param] = params[!,param][1]
    end
end



p[:εᵧ⁺] = 0.1
p[:εᵧ⁻] = 0.1
p.ω = 0.026

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
#chain = sample(model, MH(), 1000)


#d = Dict("chain" => chain)


#save("$(VAH)_HMC.jld2", d)
#sleep(5)
V = cycle_array_evaluator(p)

figure(1)
pygui(true)
clf()
plot(df.times, V)
plot(df.times, df.EcellV)



