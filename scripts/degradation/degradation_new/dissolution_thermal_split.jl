using Pkg
Pkg.activate(".")
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
const cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
const cellgeometry = CellFitElectrolyte.cell_geometry()

p_phys = ComponentVector{Any}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0,εₑ⁻ = 0, εₑ⁺ = 0, frac_sol_am_pos=0, frac_sol_am_neg=0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)

cell = "VAH"*lpad(parse(Int, ENV["SLURM_ARRAY_TASK_ID"]),2,"0")
cells = [cell,]


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

function lifetime_evaluator(p::ComponentVector{T}, cycle_array_vec, cycles_to_save, eps_0) where {T}
    # Handle Initial Conditions
    u::Array{T,1}  = zeros(T, 10)
    CellFitElectrolyte.initial_conditions!(u,p.p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
    Temp = p.p_thermal.T_amb
    @unpack R⁻ = p.p_phys
    @unpack frac_sol_am_pos = p.p_phys

    u[9] = Temp
    u[10] = eps_0
    εₑ⁺ = eps_0
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁺ = 1 - εₑ⁺ - εₛ⁺
    

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
    c_s_p_max = cathodeocv.c_s_max
    c_s_n_max = anodeocv.c_s_max

    function f(du, u, p, t)
        #p_phys = @view p[:p_phys]
        #p_deg = @view p[:p_deg]
        @unpack p_phys, p_deg, p_thermal = p
        Temp = u[9]
        Iapp = u[8]
        @pack! p_phys = Temp
        u[10] = εₑ⁺
        εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
        εᵧ⁺ = 1 - εₑ⁺ - εₛ⁺
        p.p_phys.εᵧ⁺ = εᵧ⁺
        p.p_phys.εₛ⁺ = εₛ⁺
        p.p_phys.εₑ⁺ = εₑ⁺

        #COUPLING
        #for (i, param) in enumerate(predicted_states)
        #    p_phys[param] = u[9+i]
        #end
        CellFitElectrolyte.equations_electrolyte_life_allocating((@view du[1:8]), (@view u[1:8]), p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
        @unpack εₛ⁻,εₛ⁺,εₑˢ = p_phys
        @unpack R⁻,R⁺, t⁺, ω, β⁻, β⁺, κ = p_phys
        @unpack εₛ⁻,εₛ⁺,εₑˢ = p.p_phys
        @unpack εᵧ⁻ = p.p_phys
        
        εₑ⁻ = 1-εₛ⁻-εᵧ⁻
        εₑ⁺ = 1-εₛ⁺-εᵧ⁺
        #get stuff for SEI model out of state vector
        cₑ⁻ = u[3]
        cₛˢ⁻ = u[1]
        cₑ⁺ = u[5]
        cₛˢ⁺ = u[7]
        cₑˢ = u[4]
    
    
        #Thermals
        a⁻ = 3 *εₛ⁻/R⁻
        A⁻ = 2 *cellgeometry.Vₛ⁻*a⁻
        a⁺ = 3 *εₛ⁺/R⁺
        A⁺ = 2 *cellgeometry.Vₛ⁺*a⁺
        J⁻ = Iapp/A⁻
        J⁺ = Iapp/A⁺
        @unpack C, ha, hb, T_amb = p_thermal
    
        J₀⁻ = CellFitElectrolyte.exchange_current_density(cₛˢ⁻,cₑ⁻,c_s_n_max,p_phys.k₀⁻)
        J₀⁺ = CellFitElectrolyte.exchange_current_density(cₛˢ⁺,cₑ⁺,c_s_p_max,p_phys.k₀⁺)
    
        η₊ = CellFitElectrolyte.butler_volmer(J₀⁺,J⁺,Temp)
        η₋ = CellFitElectrolyte.butler_volmer(J₀⁻,J⁻,Temp)
    
        ηc₋ = CellFitElectrolyte.concentration_overpotential(cₑ⁻,cₑˢ,t⁺,Temp,cellgeometry.T⁻)
        ηc₊ = CellFitElectrolyte.concentration_overpotential(cₑˢ,cₑ⁺,t⁺,Temp,cellgeometry.T⁺)
    
        ηₒ₋ = CellFitElectrolyte.electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,cellgeometry.T⁻)
        ηₒ₊ = CellFitElectrolyte.electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,cellgeometry.T⁺)
    
        ηₛ = CellFitElectrolyte.sei_ohmic(ω, Iapp)
    
        η = η₊+η₋+ηc₋+ηc₊+ηₒ₋+ηₒ₊+ηₛ
    
        #println(η)
    
        Q_gen = Iapp*η
        h = ha*(Temp - T_amb) + hb
        Q_cool = h*(Temp - T_amb)
        du[9] = (Q_gen - Q_cool)/C
        
        #DISSOLUTION MODEL
        @unpack i0_diss, U, α_diss  = p_deg
        xₛˢ⁺ = (cₛˢ⁺)/(c_s_p_max)
        q = c_s_p_max*cellgeometry.T⁺
        U⁺ = calcocv(cathodeocv, xₛˢ⁺, Temp)
        surface_potential = U⁺ - η₊ - U

        i_diss = CellFitElectrolyte.active_material_dissolution(i0_diss, a⁺, α_diss, Temp, surface_potential)
        du[10] = CellFitElectrolyte.active_material_dissolution_volfrac(i_diss, q, 1)

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

integrator = init(prob,QNDF(autodiff=false),save_everystep=false)
integrator.opts.maxiters = 1e7
#integrator.opts.dtmax = 1.0

#simulate
for (i, cycle_array) in enumerate(cycle_array_vec)
    T_amb = cycle_array[end]
    @pack! p.p_thermal = T_amb
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
            if integrator.sol.retcode != :Default
                return zeros(length(integrator.u))
            end
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,@view(integrator.p[:p_phys]),integrator.t,cache,cellgeometry,cathodeocv,anodeocv,integrator.u[8])
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
                            if integrator.sol.retcode !=:Default
                                return zeros(length(integrator.u))
                            end
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
                        if integrator.sol.retcode != :Default
                            return zeros(length(integrator.u), length(cycles_to_save))
                        end
                    end
                end
            
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
    err = degradation_simulator(cell, p_deg)
    num = length(data_dict[cell]["cycles"])
    return sqrt(err/num)
end

function degradation_simulator(cell, p_deg)
    err = 0
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

    if cell == "VAH07"
        vfull = 4.0
    elseif cell == "VAH23"
        vfull = 4.1
    else
        vfull = 4.2
    end

    eps_0 = εₑ⁺

    p_phys = ComponentVector{Any}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻,εₑ⁻ = εₑ⁻, εₑ⁺ = εₑ⁺, frac_sol_am_pos=frac_sol_am_pos, frac_sol_am_neg=frac_sol_am_neg, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)
    p_thermal = ComponentVector(T_amb = cycle_array[1][end], C = 48.52559886527583, ha=0.005659969188357159, hb = 0.011298387604847037)

    ics = ComponentVector(ω=ω, εₑ⁻=εₑ⁻, εₑ⁺=εₑ⁺, frac_sol_am_neg=frac_sol_am_neg, frac_sol_am_pos=frac_sol_am_pos)
    p = ComponentVector(p_deg=p_deg, p_phys = p_phys, ics=ics, p_thermal=p_thermal)

    u_out = lifetime_evaluator(p, cycle_array, data_dict[cell]["cycles"], eps_0)
    sorted_cycles = sort(data_dict[cell]["cycles"])
    for (i, cycle) in enumerate(sorted_cycles)
        #fit based on volume fraction
        err += (data_dict[cell]["mean"][cycle][:εₑ⁺] - u_out[10, i]).^2
    end
    return err
end


function loss(p_deg::T) where {T <: ComponentArray}
    try
        sol = fit_cfe_degradation(cycle_array_vec, data_dict,cells, p_deg)
        return sol
    catch
        return 1e12
    end
end

function loss(p_deg)
    i0_diss, U, α_diss  = p_deg
    p_deg = ComponentArray(
        i0_diss = i0_diss,
        U = U,
        α_diss = α_diss

    )
    return loss(p_deg)
end

i0_diss = 1e-18
U = 4.6
α_diss = 0.5

p_deg = ComponentArray(
    i0_diss = i0_diss,
    U = U,
    α_diss = α_diss,
)

ub = ComponentArray(
    i0_diss = 1e-17,
    U = 5.0,
    α_diss = 1.0
)


lb = ComponentArray(
    i0_diss = 0,
    U = 3.0,
    α_diss = 0.0
)

sol = CellFitElectrolyte.anneal(loss, p_deg, ub, lb)

#@save "sei_$(cell)_fullycoupled.jld2" sol