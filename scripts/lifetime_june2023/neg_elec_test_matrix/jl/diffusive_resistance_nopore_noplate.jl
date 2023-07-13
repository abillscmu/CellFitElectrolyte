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
using MultiKDE

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

@load "files_to_use.jld2"
dataset = ARGS[1]
cells_to_use = files_to_use[dataset]
#cells_to_use = ["VAH01",]
fitting_cycles = Dict(k => files_to_use[k] for k in cells_to_use)


predicted_states = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]
param_to_idx = Dict(
    :ω => 9,
    :εₑ⁻=>10,
    :εₑ⁺=>11,
    :frac_sol_am_neg=>12,
    :frac_sol_am_pos=>13
)

params = (:ω, :εₑ⁻, :frac_sol_am_neg)

function chainoperator(chain, param, operator)
    return operator(chain[param].data[:,1])
end

function scottoperator(data)
    σ = std(data)
    n = length(data)
    return 1.06*σ*(n^(-1/5))
end


# Load data and build kernel density estimates
FOLDERNAME = "results/0302/"
distribution_dict = Dict()
dims = [ContinuousDim() for param in params]
for vah in keys(fitting_cycles)
    distribution_dict[vah] = Dict()
    for cycle in fitting_cycles[vah]
        #Load All Raw Data
        chain = load(FOLDERNAME*"$(vah)_$(cycle)_HMC.jld2")["chain"]
        N = length(chain)
        distribution_dict[vah][cycle] = Dict()
        for sym in predicted_states
            distribution_dict[vah][cycle][sym] = chain[sym].data[:,1]
        end
        #Get observations for the multivariate KDE
        observations = Vector{Vector{Float64}}(undef, N)
        for i in 1:N
            observation_vec = zeros(length(params))
            for (j,param) in enumerate(params)
                observation_vec[j] = chain[param].data[i,1]
            end
            observations[i] = observation_vec
        end
        distribution_dict[vah][cycle]["observations"] = observations
        # Now calculate the bandwidths for the multivariate KDE
        distribution_dict[vah][cycle]["bandwidths"] = Dict()
        for (j, param) in enumerate(params)
            distribution_dict[vah][cycle]["bandwidths"][param] = chainoperator(chain, param, scottoperator)
        end
        # Now do the kde estimation
        bws = [distribution_dict[vah][cycle]["bandwidths"][param] for param in params]
        kde = KDEMulti(dims, bws, observations)
        distribution_dict[vah][cycle]["kde"] = kde
    end
end

#Load CycleArrays
cycle_array_vec = CellFitElectrolyte.load_airbus_cyclearrays()["cycle_array_vector"]



function lifetime_evaluator(p, cycle_array_vec, u, saves, k_sei, E_sei, k_resistance)
    #Create Function (with neural network)
    println((k_sei, E_sei, k_resistance))
    function f(du, u, p, t)
        @unpack p_nn, p_phys, p_thermal = p
        Voltage, U⁺, U⁻, η⁺, η⁻ = CellFitElectrolyte.equations_electrolyte_allocating_new_withvoltage!(du, u, p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
        du[9:end] .= 0
        j_sei = CellFitElectrolyte.sei_diffusive(k_sei, u[9], E_sei, u[14])
        du[9] = k_resistance * j_sei
        du[12] = - j_sei
        # Thermal model
        OCV = U⁺ - U⁻
        CellFitElectrolyte.thermal_model(du, u, p_thermal, Voltage, OCV)
        return nothing
    end

    #Construct Problem
    integrator = CellFitElectrolyte.construct_odeproblem(f, u, p)

    #simulate
    for (i, cycle_array) in enumerate(cycle_array_vec)
        p.p_thermal.T_amb = cycle_array[end]
        num_steps = Int(cycle_array[1])
        #Reset cccv switch
        integrator.p.p_phys.cccv_switch_2 = false
        integrator.p.p_phys.cccv_switch = false
        if num_steps == -1.0
            success = CellFitElectrolyte.simulate_rpt!(integrator, cycle_array, cache, cellgeometry, cathodeocv, anodeocv)
            if !(success)
                println("Failed at step $i due to $(integrator.sol.retcode)")
                println(integrator.u)
                return 0
            end
        else
            success = CellFitElectrolyte.simulate_normal_cycle!(integrator, cycle_array, cache, cellgeometry, cathodeocv, anodeocv, num_steps)
            if !(success)
                println("Failed at step $i due to $(integrator.sol.retcode)")
                println(integrator.u)
                return 0
            end
        end
        if i in saves
            savevalues!(integrator, true)
        end
    end
    return Array(integrator.sol)
    #return integrator.sol
end

@model function turing_life_fit(distribution_dict, cycle_array_vec, lifetime_evaluator, fitting_cycles)
    k_sei_exponent ~ Uniform(-16, -14)
    k_sei = 1*10. ^k_sei_exponent
    E_sei ~ Uniform(0, 10000)
    k_resistance_exponent ~ Uniform(-2, 1)
    k_resistance = 1*10. ^k_resistance_exponent


    for vah in keys(fitting_cycles)
        cycles = fitting_cycles[vah]
        first_cycle = cycles[1]
        last_cycle = cycles[end]
        cellnum = parse(Int, vah[4:end])

        ω = mean(distribution_dict[vah][first_cycle][:ω])
        εₑ⁻ = mean(distribution_dict[vah][first_cycle][:εₑ⁻])
        εₑ⁺ = mean(distribution_dict[vah][first_cycle][:εₑ⁺])
        frac_sol_am_neg = mean(distribution_dict[vah][first_cycle][:frac_sol_am_neg])
        frac_sol_am_pos = mean(distribution_dict[vah][first_cycle][:frac_sol_am_pos])

        εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
        εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
        εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
        εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

        p_phys = ComponentVector{Float64}(θₛ⁻ = 6.130391775012598e-10, θₑ = 3e-6, θₛ⁺ = 3.728671559985511e-8, R⁺ = 4.2902932816468984e-6, R⁻ = 6.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 1e-1, k₀⁻ = 1e-1, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 2000.0, κ = 0.25, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)
        p_thermal = ComponentVector(T_amb = cycle_array_vec[cellnum][first_cycle][end], C = 48.52559886527583, ha=0.005659969188357159, hb = 0.011298387604847037)


        # Handle Initial Conditions
        p_nn = zeros(10)
        T = eltype(p_nn)
        p = ComponentArray{Float64}(p_phys = p_phys, p_nn = p_nn, p_thermal=p_thermal)
        
        Temp = 320.0
        @pack! p.p_phys = Temp
        input_type = 0
        input_value = 0
        @pack! p.p_phys = input_type,input_value
        saves = fitting_cycles[vah] .- first_cycle .+ 1

        u_top = vcat(ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos)
        u_bot = zeros(typeof(k_resistance), 8)
        p_phys = p.p_phys
        CellFitElectrolyte.initial_conditions!(u_bot,p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
        u = vcat(u_bot, u_top, cycle_array_vec[cellnum][first_cycle][end])
        sol = try 
            sol = lifetime_evaluator(p, cycle_array_vec[cellnum][first_cycle:last_cycle], u, saves, k_sei, E_sei, k_resistance)
        catch
            @info "Solver errored"
            sol = 0
        end
        if sol == 0
            Turing.@addlogprob! -Inf
            return nothing
        end
        ω_end = sol[12,end]
        println("ending frac: $ω_end")
        for (i,cycle) in enumerate(cycles)
            ending_vals = [sol[param_to_idx[p],i+1] for (j,p) in enumerate(params)]
            logprob = log(MultiKDE.pdf(distribution_dict[vah][cycle]["kde"], ending_vals))
            if (logprob == -Inf) & (i > 1)
                println("Became unfeasible at step $i")
                println("Ending vals: $(ending_vals)")
                Turing.@addlogprob! -1000
            else
                println("Log probability of step $i $logprob")
                Turing.@addlogprob! logprob
            end
        end
    end
    return nothing
end


model = turing_life_fit(distribution_dict, cycle_array_vec, lifetime_evaluator, fitting_cycles)

opt = optimize(model, MAP(),[-15, 250.0, -0.75], NelderMead())
@save "lifetime_results/diffusive_resistance_nopore_noplate_checkpoint_$(dataset).jld2" opt

function run_thru(distribution_dict, cycle_array_vec, lifetime_evaluator, fitting_cycles, k_sei_exponent, E_sei, k_resistance_exponent)

    k_resistance = 1*10. ^ k_resistance_exponent
    k_sei = 1*10. ^ k_sei_exponent
    sol_dict = Dict()


    for vah in keys(fitting_cycles)
        cycles = fitting_cycles[vah]
        first_cycle = cycles[1]
        last_cycle = cycles[end]
        cellnum = parse(Int, vah[4:end])

        ω = mean(distribution_dict[vah][first_cycle][:ω])
        εₑ⁻ = mean(distribution_dict[vah][first_cycle][:εₑ⁻])
        εₑ⁺ = mean(distribution_dict[vah][first_cycle][:εₑ⁺])
        frac_sol_am_neg = mean(distribution_dict[vah][first_cycle][:frac_sol_am_neg])
        frac_sol_am_pos = mean(distribution_dict[vah][first_cycle][:frac_sol_am_pos])

        εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
        εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
        εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
        εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

        p_phys = ComponentVector{Float64}(θₛ⁻ = 6.130391775012598e-10, θₑ = 3e-6, θₛ⁺ = 3.728671559985511e-8, R⁺ = 4.2902932816468984e-6, R⁻ = 6.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 1e-1, k₀⁻ = 1e-1, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 2000.0, κ = 0.25, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)
        p_thermal = ComponentVector(T_amb = cycle_array_vec[cellnum][first_cycle][end], C = 48.52559886527583, ha=0.005659969188357159, hb = 0.011298387604847037)


        # Handle Initial Conditions
        p_nn = zeros(10)
        T = eltype(p_nn)
        p = ComponentArray{Float64}(p_phys = p_phys, p_nn = p_nn, p_thermal=p_thermal)
        
        Temp = 320.0
        @pack! p.p_phys = Temp
        input_type = 0
        input_value = 0
        @pack! p.p_phys = input_type,input_value
        saves = fitting_cycles[vah] .- first_cycle .+ 1

        u_top = vcat(ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos)
        u_bot = zeros(typeof(k_resistance), 8)
        p_phys = p.p_phys
        CellFitElectrolyte.initial_conditions!(u_bot,p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
        u = vcat(u_bot, u_top, cycle_array_vec[cellnum][first_cycle][end])
        sol = try 
            sol = lifetime_evaluator(p, cycle_array_vec[cellnum][first_cycle:last_cycle], u, saves, k_sei, E_sei, k_resistance)
        catch
            sol = zeros(length(u), length(saves))
        end
        ω_end = sol[12,end]
        println("ending frac: $ω_end")
        sol_dict[vah] = Dict()
        for p in params
            sol_dict[vah][p] = [sol[param_to_idx[p],:]]
        end
    end
    return sol_dict
end

k_sei_exponent = opt.values[1]
E_sei = opt.values[2]
k_resistance_exponent = opt.values[3]

my_sol = run_thru(distribution_dict, cycle_array_vec, lifetime_evaluator, fitting_cycles, k_sei_exponent, E_sei, k_resistance_exponent)

cell_to_plot = "VAH01"

@save "lifetime_results/diffusive_resistance_nopore_noplate_$(dataset).jld2" my_sol distribution_dict fitting_cycles opt
