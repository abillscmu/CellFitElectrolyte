using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.DiffEqFlux
using CellFitElectrolyte.Turing
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

#set up neural network (inputs)
num_layers = 1
activation = tanh
width = 2
num_augmented_states = 0
predicted_states = [:ω]

#set up neural network (derived)
num_predicted_states = length(predicted_states)
input_dim = 8 + num_predicted_states + num_augmented_states
output_dim = num_predicted_states + num_augmented_states
nn = CellFitElectrolyte.build_nn(num_layers, activation, width, input_dim, output_dim)

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

p = ComponentVector{Any}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)
p_nn_init = DiffEqFlux.initial_params(nn)

#Problem is now built. Next step is to load data. This section roughly corresponds to batch settings.
cell = 1
first_cycle = 2
last_cycle = 4


#Load data
FOLDERNAME = "results/outputs1214_nuts_vah01/"
vah = lpad(cell, 2, "0")
first_chain = load(FOLDERNAME*"VAH$(vah)_$(first_cycle)_HMC.jld2")["chain"]
last_chain = load(FOLDERNAME*"VAH$(vah)_$(last_cycle)_HMC.jld2")["chain"]

#Build Distributions
first_distribution_data = first_chain[:ω].data[:,1]
last_distribution_data = last_chain[:ω].data[:,1]
first_distribution_kde = kde(first_distribution_data)
last_distribution_kde = kde(last_distribution_data)
first_distribution_bw = KernelDensity.default_bandwidth(first_distribution_data)
last_distribution_bw = KernelDensity.default_bandwidth(last_distribution_data)
first_distribution = KDEDistribution(first_distribution_data, first_distribution_kde, first_distribution_bw)
last_distribution = KDEDistribution(last_distribution_data, last_distribution_kde, last_distribution_bw)

#Load CycleArrays
cycle_array_vec = CellFitElectrolyte.load_airbus_cyclearrays()["cycle_array_vector"][cell][first_cycle:last_cycle]

function lifetime_evaluator(p::ComponentVector{T}, cycle_array_vec, u) where {T}
    #Create Function (with neural network)
    Temp = 320
    function f(du, u, p, t)
        @unpack p_nn, p_phys = p
        for (i, param) in enumerate(predicted_states)
            p_phys[param] = u[8+i]
        end
        CellFitElectrolyte.equations_electrolyte_life_allocating((@view du[1:8]), (@view u[1:8]), p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
        new_u = nn(u, p_nn)
        du[9:end] .= new_u
        return nothing
    end

    #Construct Problem
    vars = ones(length(u))
    vars[8] = 0
    mm = diagm(vars)
    func = ODEFunction(f, mass_matrix=mm)
    prob = ODEProblem(func,u,(0.0,Inf),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false)
    integrator.opts.maxiters = 1e7
    integrator.opts.dtmax = 1.0

    #simulate
    for (i, cycle_array) in enumerate(cycle_array_vec)
        num_steps = Int(cycle_array[1])
        #Reset cccv switch (note that p will not be correct with NN)
        integrator.p.p_phys.cccv_switch_2 = false
        integrator.p.p_phys.cccv_switch = false
        if num_steps == -1.0
            #RPT
            input_type = 3.0
            input_value = cycle_array[2]
            @pack! integrator.p.p_phys = input_type, input_value
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phs,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            u = deepcopy(integrator.u)
            u[8] = input_value
            OrdinaryDiffEq.set_u!(integrator, u)
            while Voltage >= 2.5
                step!(integrator)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            end
            input_type = 0.0
            input_value = 0.0
            @pack! integrator.p.p_phys = input_type, input_value
            end_time = integrator.t + cycle_array[4]
            add_tstop!(integrator, end_time)
            while integrator.t < end_time
                step!(integrator)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            end
            input_type = 5.0
            input_value = cycle_array[3]
            @pack! integrator.p.p_phys = input_type, input_value
            u = deepcopy(integrator.u)
            u[8] = input_value
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            OrdinaryDiffEq.set_u!(integrator, u)
            while integrator.u[8] < -0.01
                step!(integrator)
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
            end
        else
            #Normal array
            times = cycle_array[2:num_steps+1]
            times .+= integrator.t
            types = cycle_array[num_steps+2:num_steps+1+num_steps]
            values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
            #by definition, times[1] == integrator.t
            for t in times[2:end]
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
                    u = deepcopy(integrator.u)
                    u[8] = input_value
                    OrdinaryDiffEq.set_u!(integrator, u)
                elseif input_type == 1.0
                    u = deepcopy(integrator.u)
                    Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                    u[8] = input_value/Voltage
                    OrdinaryDiffEq.set_u!(integrator, u)
                end
                while integrator.t < end_time
                    step!(integrator)
                    #if CCCV, have to do some weird control stuff.
                    if input_type == 5
                        Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                        if ((Voltage >= p.p_phys.vfull) & ((integrator.p.p_phys.cccv_switch != true)))
                            integrator.p.p_phys.cccv_switch = true
                        end
                        if ((((integrator.u[8] >= integrator.p.p_phys.ifull) & ((integrator.p.p_phys.cccv_switch == true)) & (input_type ==5))) & (integrator.p.p_phys.cccv_switch_2 != true))
                            u = deepcopy(integrator.u)
                            #going to rest, have to reset integrator val :(
                            u[8] = 0
                            OrdinaryDiffEq.set_u!(integrator, u)
                            integrator.p.p_phys.cccv_switch_2 = true
                        end
                    end
                end
            end
        end
    end
    return integrator.u
end




@model function fit_bw(first_distribution, last_distribution, cycle_array_vec, p_phys, nn)
        # Handle Initial Conditions
        u = zeros(8)
        CellFitElectrolyte.initial_conditions!(u,p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
        augmented_states = zeros(num_augmented_states)
        ω ~ first_distribution
        p_nn ~ MvNormal(zeros(length(initial_params(nn))), 1e-10)
        p = ComponentArray(p_phys = p_phys, p_nn = p_nn)
        u = vcat(u, ω, augmented_states)
        
        Temp = 320.0
        @pack! p.p_phys = Temp
        input_type = 0
        input_value = 0
        @pack! p.p_phys = input_type,input_value
        
        #Simulate
        new_u = lifetime_evaluator(p, cycle_array_vec, u)

        #Evaluate
        ω_out = new_u[9]
        probability = Distributions.logpdf(last_distribution, ω_out)
        Turing.Turing.@addlogprob! probability
        
        return nothing
end

model = fit_bw(first_distribution, last_distribution, cycle_array_vec, p, nn)

chain = sample(model, PG(10), 1000)


