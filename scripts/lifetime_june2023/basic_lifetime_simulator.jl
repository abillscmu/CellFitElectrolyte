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
using PythonPlot

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
last_cycle = 845


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

function lifetime_evaluator(p::ComponentVector{T}, cycle_array_vec, u) where {T}
    #Create Function (with neural network)
    Temp = 320
    c_s_p_max = cathodeocv.c_s_max
    c_s_n_max = anodeocv.c_s_max
    
    function f(du, u, p, t)
        @unpack p_nn, p_phys = p
        @unpack ω = p_phys
        ω = u[9]
        εₑ⁻ = u[10]
        εₑ⁺ = u[11]
        frac_sol_am_neg = u[12]
        frac_sol_am_pos = u[13]
        εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
        εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
        εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
        εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺
        @pack! p_phys = ω, εₛ⁺,εₛ⁻, εᵧ⁻, εᵧ⁺ 
        CellFitElectrolyte.equations_electrolyte_allocating_new_withvoltage((@view du[1:8]), (@view u[1:8]), p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
        du[9:end] .= 0
        return nothing
    end

    #Construct Problem
    vars = ones(length(u))
    vars[8] = 0
    mm = diagm(vars)
    func = ODEFunction(f, mass_matrix=mm)
    prob = ODEProblem(func,u,(0.0,Inf),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, verbose=false)
    integrator.opts.maxiters = 1e7
    #integrator.opts.dtmax = 1.0

    #simulate
    for (i, cycle_array) in enumerate(cycle_array_vec)
        println(i)
        num_steps = Int(cycle_array[1])
        #Reset cccv switch (note that p will not be correct with NN)
        integrator.p.p_phys.cccv_switch_2 = false
        integrator.p.p_phys.cccv_switch = false
        if num_steps == -1.0
            integrator.opts.dtmax = 10.0
            #RPT
            input_type = 3.0
            input_value = cycle_array[2]
            @pack! integrator.p.p_phys = input_type, input_value
            Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            u = copy(integrator.u)
            u[8] = input_value
            OrdinaryDiffEq.set_u!(integrator, u)
            while Voltage >= 2.5
                step!(integrator)
                if integrator.sol.retcode != :Default
                    return zeros(length(integrator.u))
                end
                Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            end
            input_type = 0.0
            input_value = 0.0
            @pack! integrator.p.p_phys = input_type, input_value
            end_time = integrator.t + cycle_array[4]
            add_tstop!(integrator, end_time)
            while integrator.t < end_time
                step!(integrator)
                if integrator.sol.retcode != :Default
                    return zeros(length(integrator.u))
                end
                Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            end
            input_type = 5.0
            input_value = cycle_array[3]
            @pack! integrator.p.p_phys = input_type, input_value
            u = copy(integrator.u)
            u[8] = input_value
            Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            OrdinaryDiffEq.set_u!(integrator, u)
            while integrator.u[8] < -0.01
                step!(integrator)
                Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                if integrator.sol.retcode != :Default
                    return zeros(length(integrator.u))
                end
                if ((Voltage >= integrator.p.p_phys.vfull) & ((integrator.p.p_phys.cccv_switch != true)))
                    integrator.p.p_phys.cccv_switch = true
                end
            end
            end_time = integrator.t + cycle_array[5]
            add_tstop!(integrator, end_time)
            input_type = 0.0
            input_value = 0.0
            u = deepcopy(integrator.u)
            u[8] = input_value
            OrdinaryDiffEq.set_u!(integrator, u)
            @pack! integrator.p.p_phys = input_type, input_value
            while integrator.t < end_time
                step!(integrator)
                if integrator.sol.retcode != :Default
                    return zeros(length(integrator.u))
                end
            end
        else
            #Normal array
            times = cycle_array[2:num_steps+1]
            times .+= integrator.t
            types = cycle_array[num_steps+2:num_steps+1+num_steps]
            values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
            #by definition, times[1] == integrator.t
            for t in (@view times[2:end])
                add_tstop!(integrator,t)
            end
            #Step through the cycle
            for step::Int in 1:num_steps-1
                @pack! p.p_phys = Temp
                input_type = types[step]
                input_value = values[step]
                end_time::T = times[step+1]
                @pack! p.p_phys = input_type,input_value
                #make sure we can deal with the initial condition
                if (input_type == 0.0) | (input_type == 5.0)
                    u = copy(integrator.u)
                    u[8] = input_value
                    OrdinaryDiffEq.set_u!(integrator, u)
                elseif input_type == 1.0
                    u = copy(integrator.u)
                    Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                    u[8] = input_value/Voltage
                    OrdinaryDiffEq.set_u!(integrator, u)
                end
                while integrator.t < end_time
                    #if CCCV, have to do some weird control stuff.
                    if input_type == 5
                        Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                        if (Voltage < p.p_phys.vfull * 0.8)
                            integrator.opts.dtmax = 10.0
                        else
                            integrator.opts.dtmax = 5.0
                        end
                        step!(integrator)
                        if ((Voltage >= p.p_phys.vfull) & ((integrator.p.p_phys.cccv_switch != true)))
                            integrator.p.p_phys.cccv_switch = true
                        end
                        if ((((integrator.u[8] >= integrator.p.p_phys.ifull) & ((integrator.p.p_phys.cccv_switch == true)) & (input_type ==5))) & (integrator.p.p_phys.cccv_switch_2 != true))
                            u = copy(integrator.u)
                            #going to rest, have to reset integrator val :(
                            u[8] = 0
                            OrdinaryDiffEq.set_u!(integrator, u)
                            integrator.p.p_phys.cccv_switch_2 = true
                        end
                    else
                        integrator.opts.dtmax = end_time - integrator.t
                        step!(integrator)
                    end
                    if integrator.sol.retcode != :Default
                        return zeros(length(integrator.u))
                    end
                end
            end
        end
    end
    return integrator.u
end



# Handle Initial Conditions
p_nn = zeros(10)
T = eltype(p_nn)
p = ComponentArray{Float64}(p_phys = p_phys, p_nn = p_nn)
u = zeros(8)

CellFitElectrolyte.initial_conditions!(u,p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)

    
u = vcat(u, ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos)
        
Temp = 320.0
@pack! p.p_phys = Temp
input_type = 0
input_value = 0
@pack! p.p_phys = input_type,input_value

sol = lifetime_evaluator(p, cycle_array_vec, u)

