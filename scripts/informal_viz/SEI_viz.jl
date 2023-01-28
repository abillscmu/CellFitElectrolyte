using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.DiffEqFlux
using CellFitElectrolyte.Turing
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

p_phys = ComponentVector{Any}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0,εₑ⁻ = 0, εₑ⁺ = 0, frac_sol_am_pos=0, frac_sol_am_neg=0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)

cells = ["VAH01"]

#Load data
FOLDERNAME = "results/outputs0117_elec/"

#Build Distributions
data_dict = Dict()
predicted_states = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]
for file in readdir(FOLDERNAME)
        vah = split(file, "_")[1]
        if !(vah in cells)
            continue
        end 
        filename = FOLDERNAME*file
        chain = try
            d = load(filename)
            chain = d["chain"]
        catch
            @warn "problem with file $filename"
            continue
        end
        cycle_num = parse(Int, split(file, "_")[2])
        if vah in keys(data_dict)
            append!(data_dict[vah]["cycles"], cycle_num)
        else
            data_dict[vah] = Dict(
                "cycles" => [cycle_num,],
                "distributions" => Dict(),
                "mean" => Dict()
            )
        end
        frac_data = chain[:frac_sol_am_neg].data[:,1]
        eps_data = chain[:εₑ⁻].data[:,1]
        eps_g = (1 .- frac_data).*(1 .- eps_data)
        eps_g_kde = kde(eps_g)
        eps_g_distribution_bw = KernelDensity.default_bandwidth(eps_g)
        eps_g_distribution = KDEDistribution(eps_g, eps_g_kde, eps_g_distribution_bw) 


        data_dict[vah]["distributions"][cycle_num] = Dict(
            :εᵧ⁻ => eps_g_distribution
        )
        data_dict[vah]["mean"][cycle_num] = Dict(
            :εᵧ⁻ => mean(eps_g)
        )
        for sym in predicted_states
            distribution_data = chain[sym].data[:,1]
            distribution_kde = kde(distribution_data)
            distribution_bw = KernelDensity.default_bandwidth(distribution_data)
            distribution = KDEDistribution(distribution_data, distribution_kde, distribution_bw)
            data_dict[vah]["distributions"][cycle_num][sym] = distribution
            data_dict[vah]["mean"][cycle_num][sym] = mean(distribution_data)
        end
end

#Load CycleArrays
cycle_array_vec = CellFitElectrolyte.load_airbus_cyclearrays()["cycle_array_vector"]

#u[13] => δ⁻
#u[14] => δ⁺

function lifetime_evaluator(p::ComponentVector{T}, cycle_array_vec, cycles_to_save, eps_0, f_0) where {T}
    # Handle Initial Conditions
    u::Array{T,1}  = zeros(T, 10)
    CellFitElectrolyte.initial_conditions!(u,p.p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
    Temp = 300
    @unpack R⁻ = p.p_phys
    d_0 = R⁻*(f_0 + 0.03)^(-1/3) - R⁻

    u[9] = d_0
    u[10] = eps_0

    cycle_array = cycle_array_vec[1]
    num_steps = Int(cycle_array[1])
    times = cycle_array[2:num_steps+1]
    types = cycle_array[num_steps+2:num_steps+1+num_steps]
    values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

    @pack! p.p_phys = Temp
    du = similar(u)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p.p_phys = input_type,input_value
#Create Function (with neural network)
Temp = 320
c_s_p_max = cathodeocv.c_s_max
c_s_n_max = anodeocv.c_s_max

function f(du, u, p, t)
    #p_phys = @view p[:p_phys]
    #p_deg = @view p[:p_deg]
    @unpack p_phys, p_deg = p
    #COUPLING
    #for (i, param) in enumerate(predicted_states)
    #    p_phys[param] = u[8+i]
    #end
    CellFitElectrolyte.equations_electrolyte_life_allocating((@view du[1:8]), (@view u[1:8]), p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
    @unpack εₛ⁻,εₛ⁺,εₑˢ,Temp = p_phys
    @unpack R⁻ = p_phys
    @unpack εₑ⁻ = p_phys
    #get stuff for SEI model out of state vector
    cₑ⁻ = u[3]
    cₛˢ⁻ = u[1]
    δ = u[9]
    Iapp=u[8]
    xₛˢ⁻ = (cₛˢ⁻)/(c_s_n_max)
    #OCV
    U⁻ = calcocv(anodeocv, xₛˢ⁻, Temp)
    #Geometry 
    a⁻::T = 3 *εₛ⁻/R⁻
    A⁻::T = 2 *cellgeometry.Vₛ⁻*a⁻
    J⁻ = Iapp/A⁻
    J₀⁻ = CellFitElectrolyte.exchange_current_density(cₛˢ⁻,cₑ⁻,c_s_n_max,p_phys.k₀⁻)
    η⁻ = CellFitElectrolyte.butler_volmer(J₀⁻,J⁻,Temp)
    surface_potential = U⁻ + η⁻
    @unpack εᵧ⁻ = p_phys
    a_sei = 3*(εᵧ⁻ + εₛ⁻)/(R⁻ + δ)
    #unpack 
    @unpack c_s_bulk_neg, k_sei_neg, α_sei_neg, D_sei, U_sei, Vₘ = p_deg
    #SEI OP
    η_sei = surface_potential - U_sei
    
    #Calculate SEI model(J is molar current)!!
    J_sei = CellFitElectrolyte.sei_growth(c_s_bulk_neg, k_sei_neg, α_sei_neg, Temp, η_sei, D_sei, δ)
    δ̇ = CellFitElectrolyte.sei_delta_growth(J_sei, δ, R⁻, εₑ⁻, Vₘ)
    #grow δ
    du[9] = δ̇
    du[10] = a_sei*δ̇
    #remove Li
    du[1] = du[1] - J_sei/(2*CellFitElectrolyte.F*cellgeometry.Vₛ⁻*εₛ⁻)
    return nothing
end

#f(du, u, p, 0.0)
#jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> f(du, u, p, 0.0), du, u)

#Construct Problem
vars = ones(length(u))
vars[8] = 0
mm = diagm(vars)
func = ODEFunction(f, mass_matrix=mm)
prob = ODEProblem(func,u,(0.0,Inf),p)

integrator = init(prob,QNDF(),save_everystep=false)
integrator.opts.maxiters = 1e7
#integrator.opts.dtmax = 1.0

#simulate
for (i, cycle_array) in enumerate(cycle_array_vec)
    num_steps = Int(cycle_array[1])
    #Reset cccv switch (note that p will not be correct with NN)
    integrator.p.p_phys.cccv_switch_2 = false
    integrator.p.p_phys.cccv_switch = false
    if num_steps == -1.0
        integrator.opts.dtmax = 15.0
        #RPT
        input_type = 3.0
        input_value = cycle_array[2]
        @pack! integrator.p.p_phys = input_type, input_value
        Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
        u = deepcopy(integrator.u)
        u[8] = input_value
        OrdinaryDiffEq.set_u!(integrator, u)
        while Voltage >= 2.5
            step!(integrator)
            #if integrator.sol.retcode != :Default
            #    return zeros(length(integrator.u))
            #end
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
        end
        input_type = 0.0
        input_value = 0.0
        @pack! integrator.p.p_phys = input_type, input_value
        end_time = integrator.t + cycle_array[4]
        add_tstop!(integrator, end_time)
        while integrator.t < end_time
            step!(integrator)
            #if integrator.sol.retcode != :Default
            #    return zeros(length(integrator.u))
            #end
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
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
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
            #if integrator.sol.retcode != :Default
            #    return zeros(length(integrator.u))
            #end
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
        for t in times[2:end]
            add_tstop!(integrator,t)
        end
        #Step through the cycle
        for step::Int in 1:num_steps-1
            @pack! (@view p[:p_phys]) = Temp
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
                Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                u[8] = input_value/Voltage
                OrdinaryDiffEq.set_u!(integrator, u)
            end
            
                #if CCCV, have to do some weird control stuff.
                if input_type == 5
                    Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                    while integrator.t < end_time
                        #println(Voltage)
                        if Voltage < 3.9
                            integrator.opts.dtmax = 10.0
                        else
                            integrator.opts.dtmax = 10
                        end
                            step!(integrator)
                            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
                        if ((Voltage >= p.p_phys.vfull) & ((integrator.p.p_phys.cccv_switch != true)))
                            integrator.p.p_phys.cccv_switch = true
                        end
                        if ((((integrator.u[8] <= integrator.p.p_phys.ifull) & ((integrator.p.p_phys.cccv_switch == true)) & (input_type ==5))) & (integrator.p.p_phys.cccv_switch_2 != true))
                            u = deepcopy(integrator.u)
                            #going to rest, have to reset integrator val :(
                            u[8] = 0
                            OrdinaryDiffEq.set_u!(integrator, u)
                            integrator.p.p_phys.cccv_switch_2 = true
                        end
                    end
                else
                    integrator.opts.dtmax = end_time - integrator.t
                    while integrator.t < end_time
                        step!(integrator)
                    end
                end
                #if integrator.sol.retcode != :Default
                #    return zeros(length(integrator.u), length(cycles_to_save))
                #end
            
        end
    end
    if i in cycles_to_save
        savevalues!(integrator, true)
    end
end
return integrator.sol
end

function fit_cfe_degradation(cycle_array_vec, data_dict, cells, p_deg)
    #(this would be rand for SG)
    err = 0.0
    num = 0
    for cell in cells
        cell_integer = parse(Int, split(cell, "VAH")[2])
        cycle_array = cycle_array_vec[cell_integer][2:end]

        ω = mean(data_dict[cell]["distributions"][2][:ω].data)
        εₑ⁻ = mean(data_dict[cell]["distributions"][2][:εₑ⁻].data)
        εₑ⁺ = mean(data_dict[cell]["distributions"][2][:εₑ⁺].data)
        frac_sol_am_neg = mean(data_dict[cell]["distributions"][2][:frac_sol_am_neg].data)
        frac_sol_am_pos = mean(data_dict[cell]["distributions"][2][:frac_sol_am_pos].data)
    
        εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
        εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
        εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
        εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

        eps_0 = εᵧ⁻
        f_0 = frac_sol_am_neg

        p_phys = ComponentVector{Any}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0,εₑ⁻ = 0, εₑ⁺ = 0, frac_sol_am_pos=0, frac_sol_am_neg=0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)

        ics = ComponentVector(ω=ω, εₑ⁻=εₑ⁻, εₑ⁺=εₑ⁺, frac_sol_am_neg=frac_sol_am_neg, frac_sol_am_pos=frac_sol_am_pos)
        p = ComponentVector(p_deg=p_deg, p_phys = p_phys, ics=ics)

        u_out = lifetime_evaluator(p, cycle_array, data_dict[cell]["cycles"], eps_0, f_0)
        sorted_cycles = sort(data_dict[cell]["cycles"])
        return u_out
        for (i, cycle) in enumerate(sorted_cycles)
            #fit based on volume fraction
            err += (data_dict[cell]["mean"][cycle][:εᵧ⁻] - u_out[10, i]).^2
            num += 1
        end
    end
    return sqrt(err/num)
end

cycles_to_save = 

model(p_deg) = fit_cfe_degradation(cycle_array_vec, data_dict, cells, p_deg)


k_sei_neg = 1e-15
D_sei = 1e-22
c_s_bulk_neg = 1000
Vₘ = 1.3e-5
U_sei = 3.5

p_deg = ComponentArray(
    c_s_bulk_neg = c_s_bulk_neg,
    k_sei_neg = k_sei_neg,
    α_sei_neg = 0.5,
    D_sei = D_sei,
    U_sei = U_sei,
    Vₘ = Vₘ
)

sol = model(p_deg)
