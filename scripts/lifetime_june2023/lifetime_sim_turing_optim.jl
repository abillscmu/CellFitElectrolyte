using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.DiffEqFlux
using CellFitElectrolyte.Turing
using BenchmarkTools
using StaticArrays
using CSV
using DataFrames
using Test
using ProgressMeter
using LinearAlgebra
using Statistics
using JLD2
using KernelDensity
using KDEDistributions
using Turing
using PythonPlot
using Optim

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

#Problem is now built. Next step is to load data. This section roughly corresponds to batch settings.
cell = 1
first_cycle = 2
last_cycle = 100


#Load data
FOLDERNAME = "results/0302/"
vah = lpad(cell, 2, "0")
first_chain = load(FOLDERNAME*"VAH$(vah)_$(first_cycle)_HMC.jld2")["chain"]
last_chain = load(FOLDERNAME*"VAH$(vah)_$(last_cycle)_HMC.jld2")["chain"]

#Build Distributions
distribution_dict = Dict()

predicted_states = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]

for sym in predicted_states
    first_distribution_data = first_chain[sym].data[:,1]
    last_distribution_data = last_chain[sym].data[:,1]
    first_distribution_kde = kde(first_distribution_data)
    last_distribution_kde = kde(last_distribution_data)
    first_distribution_bw = KernelDensity.default_bandwidth(first_distribution_data)
    last_distribution_bw = KernelDensity.default_bandwidth(last_distribution_data)
    first_distribution = KDEDistribution(first_distribution_data, first_distribution_kde, first_distribution_bw)
    last_distribution = KDEDistribution(last_distribution_data, last_distribution_kde, last_distribution_bw)
    distribution_dict[sym] = Dict()
    distribution_dict[sym]["Initial"] = first_distribution
    distribution_dict[sym]["Final"] = last_distribution
end

ω = mean(distribution_dict[:ω]["Initial"].data)
εₑ⁻ = mean(distribution_dict[:εₑ⁻]["Initial"].data)
εₑ⁺ = mean(distribution_dict[:εₑ⁺]["Initial"].data)
frac_sol_am_neg = mean(distribution_dict[:frac_sol_am_neg]["Initial"].data)
frac_sol_am_pos = mean(distribution_dict[:frac_sol_am_pos]["Initial"].data)

εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺



p_phys = ComponentVector{Float64}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0,εₑ⁻ = 0, εₑ⁺ = 0, frac_sol_am_pos=0, frac_sol_am_neg=0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)

p_phys = ComponentVector{Float64}(θₛ⁻ = 6.130391775012598e-10, θₑ = 3e-6, θₛ⁺ = 3.728671559985511e-8, R⁺ = 4.2902932816468984e-6, R⁻ = 6.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 1e-1, k₀⁻ = 1e-1, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 2000.0, κ = 0.25, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)

#Load CycleArrays
cycle_array_vec = CellFitElectrolyte.load_airbus_cyclearrays()["cycle_array_vector"][cell][first_cycle:last_cycle]

function lifetime_evaluator(p, cycle_array_vec, u, k_resistance)
    #Create Function (with neural network)
    println(k_resistance)
    function f(du, u, p, t)
        @unpack p_nn, p_phys = p
        CellFitElectrolyte.equations_electrolyte_allocating_new_withvoltage(du, u, p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
        du[9:end] .= 0
        du[9] = k_resistance
        return nothing
    end

    #Construct Problem
    integrator = CellFitElectrolyte.construct_odeproblem(f, u, p)

    #simulate
    for (i, cycle_array) in enumerate(cycle_array_vec)
        num_steps = Int(cycle_array[1])
        #Reset cccv switch
        integrator.p.p_phys.cccv_switch_2 = false
        integrator.p.p_phys.cccv_switch = false
        if num_steps == -1.0
            success = CellFitElectrolyte.simulate_rpt!(integrator, cycle_array, cache, cellgeometry, cathodeocv, anodeocv)
            if !(success)
                return zeros(size(integrator.u))
            end
        else
            success = CellFitElectrolyte.simulate_normal_cycle!(integrator, cycle_array, cache, cellgeometry, cathodeocv, anodeocv, num_steps)
            if !(success)
                return zeros(size(integrator.u))
            end
        end
    end
    return integrator.u
end

@model function turing_life_fit(distribution_dict, cycle_array_vec, p, lifetime_evaluator, ω_data)
    ω = mean(distribution_dict[:ω]["Initial"].data)
    k_resistance ~ Normal(0, 1e-9)
    u_top = vcat(ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos)
    u_bot = zeros(typeof(k_resistance), 8)
    p_phys = p.p_phys
    CellFitElectrolyte.initial_conditions!(u_bot,p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
    u = vcat(u_bot, u_top)
    sol = try 
        sol = lifetime_evaluator(p, cycle_array_vec, u, k_resistance)
    catch
        sol = zeros(length(u))
    end
    ω_end = sol[9]
    Turing.@addlogprob! loglikelihood(distribution_dict[:ω]["Initial"], ω_end)
    return nothing
end
    



# Handle Initial Conditions
p_nn = zeros(10)
T = eltype(p_nn)
p = ComponentArray{Float64}(p_phys = p_phys, p_nn = p_nn)

u_top = vcat(ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos)
u_bot = zeros(eltype(u_top), 8)
CellFitElectrolyte.initial_conditions!(u_bot,p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)

    
u = vcat(u_bot,u_top)
        
Temp = 320.0
@pack! p.p_phys = Temp
input_type = 0
input_value = 0
@pack! p.p_phys = input_type,input_value

k_resistance = 1e-10


sol = lifetime_evaluator(p, cycle_array_vec, u, k_resistance)

ω_data = distribution_dict[:ω]["Final"]


model = turing_life_fit(distribution_dict, cycle_array_vec, p, lifetime_evaluator, ω_data)



opt = optimize(model, MLE(), SimulatedAnnealing())


