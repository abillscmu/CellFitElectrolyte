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

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

function ood(u, p, t)
    return any(u .< 0)
end

cell = 1

#set up cycle arrays
println("$(pwd())")
@load "src/cycle_array_vector_new.jld2"
cycle_array_vec = cycle_array_vector[cell]
cycle_array = cycle_array_vec[1]
num_steps = Int(cycle_array[1])
times = cycle_array[2:num_steps+1]
types = cycle_array[num_steps+2:num_steps+1+num_steps]
g_values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]




function lifetime_evaluator(p::ComponentVector{T}) where {T}
    Varr = Float64[]
    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,8)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    u[8] = 1.0
    Temp = 320.0
    @pack! p = Temp
    du = similar(u)
    input_type::T = 0
    input_value::T = 0
    @pack! p = input_type,input_value

    vars = ones(8)
    vars[8] = 0
    mm = diagm(vars)

    #Create Function and Initialize Integrator
    func = ODEFunction((du::Array{T,1},u::Array{T,1},p::ComponentVector{T},t::T)->CellFitElectrolyte.equations_electrolyte_life(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv), mass_matrix=mm)
    prob = ODEProblem(func,u,(0.0,Inf),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false)
    default_dtmax = integrator.opts.dtmax
    integrator.opts.dtmax = 1.0
    integrator.opts.maxiters = 1e7



    @showprogress "simulating..." for (i, cycle_array) in enumerate(cycle_array_vec)
        num_steps = Int(cycle_array[1])
        if num_steps == -1.0
            integrator.p.cccv_switch = false
            input_type = 3.0
            input_value = cycle_array[2]
            @pack! integrator.p = input_type, input_value
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            u = deepcopy(integrator.u)
            u[8] = input_value
            OrdinaryDiffEq.set_u!(integrator, u)
            while Voltage >= 2.5
                step!(integrator)
                savevalues!(integrator, true)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                push!(Varr, Voltage)
            end
            input_type = 0.0
            input_value = 0.0
            @pack! integrator.p = input_type, input_value
            end_time = integrator.t + cycle_array[4]
            add_tstop!(integrator, end_time)
            while integrator.t < end_time
                step!(integrator)
                savevalues!(integrator, true)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                push!(Varr, Voltage)
            end
            input_type = 5.0
            input_value = cycle_array[3]
            @pack! integrator.p = input_type, input_value
            u = deepcopy(integrator.u)
            u[8] = input_value
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            OrdinaryDiffEq.set_u!(integrator, u)
            while integrator.u[8] < -0.01
                step!(integrator)
                savevalues!(integrator, true)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                if ((Voltage >= integrator.p.vfull) & ((integrator.p.cccv_switch != true) & (input_type ==5)))
                    integrator.p.cccv_switch = true
                end
                push!(Varr, Voltage)
            end
            end_time = integrator.t + cycle_array[5]
            add_tstop!(integrator, end_time)
            input_type = 0.0
            input_value = 0.0
            u = deepcopy(integrator.u)
            u[8] = input_value
            OrdinaryDiffEq.set_u!(integrator, u)
            @pack! integrator.p = input_type, input_value
            while integrator.t < end_time
                step!(integrator)
                savevalues!(integrator, true)
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                push!(Varr, Voltage)
            end
        else
            times = cycle_array[2:num_steps+1]
            times .+= integrator.t
            types = cycle_array[num_steps+2:num_steps+1+num_steps]
            values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
            for t in times[2:end]
                add_tstop!(integrator,t)
            end
            for step::Int in 1:num_steps-1
                integrator.p.cccv_switch_2 = false
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                @pack! p = Temp
                input_type = types[step]
                input_value = values[step]
                end_time::T = times[step+1]
                @pack! p = input_type,input_value
                if (input_type == 0.0) | (input_type == 5.0)
                    u = deepcopy(integrator.u)
                    u[8] = input_value
                    OrdinaryDiffEq.set_u!(integrator, u)
                    #integrator.opts.dtmax = 1.0
                    integrator.p.cccv_switch = false
                elseif input_type == 1.0
                    u = deepcopy(integrator.u)
                    u[8] = input_value/Voltage
                    OrdinaryDiffEq.set_u!(integrator, u)
                    integrator.p.cccv_switch = false
                end
    
                #try
                    while integrator.t < end_time
                        step!(integrator)
                        savevalues!(integrator, true)
                        Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                        if ((Voltage >= p.vfull) & ((integrator.p.cccv_switch != true) & (input_type ==5)))
                            integrator.p.cccv_switch = true
                        end
                        if ((((integrator.u[8] >= integrator.p.ifull) & ((integrator.p.cccv_switch == true)) & (input_type ==5))) & (integrator.p.cccv_switch_2 != true))
                            u = deepcopy(integrator.u)
                            u[8] = 0
                            OrdinaryDiffEq.set_u!(integrator, u)
                            integrator.p.cccv_switch_2 = true
                        end
                        push!(Varr, Voltage)
                    end
                #catch
                #    @warn "integrator error"
                #    return integrator
                #end
            end
        end
    end
    return integrator.sol, Varr
end

p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)
params = CSV.read("results/outputs1211_vah02_2/outputs1211_vah02/VAH02_10_PARAM.csv",DataFrame)
param_sym = Symbol.(names(params))

for param in param_sym[1:end-3]
    if param in keys(p)
        p[param] = params[!,param][1]
    end
end
p.ω = 0
asdf = lifetime_evaluator(p);
figure(1)
clf()
#plot(asdf[1].t, asdf[1][8,:])
#grid()
#twinx()
plot(asdf[1].t[2:end], asdf[2],"g")




