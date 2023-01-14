using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.Turing
using CellFitElectrolyte.DynamicHMC
using KernelDensity
using KDEDistributions
using DiffEqFlux
using StaticArrays
using CSV
using DataFrames
using Test
using ProgressMeter
using LinearAlgebra
using Statistics
using JLD2

#set up neural network (inputs)
num_layers = 2
activation = tanh
width = 10
num_augmented_states = 3
#legacy
predicted_states = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]


#set up neural network (derived)
num_predicted_states = length(predicted_states)
const input_dim = 7 + num_predicted_states + num_augmented_states
output_dim = num_predicted_states + num_augmented_states
nn = CellFitElectrolyte.build_nn(num_layers, activation, width, input_dim, output_dim)

cache = CellFitElectrolyte.initialize_cache(Float64)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

VAH = "VAH01_2"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cell_padded = lpad(cell, 2, "0")
cycle = parse(Int,split1[3])
first_cycle = cycle
num_dfs = 500
last_cycle = cycle + num_dfs

#idf : initial data frame

idf = CSV.read("/home/abills/Data/cycle_individual_data/VAH$(cell_padded)_$cycle.csv", DataFrame)
df_vec = [idf,]
if num_dfs >=2
for i in 2:num_dfs
    local_cycle = cycle + i - 1
    ndf = CSV.read("/home/abills/Data/cycle_individual_data/VAH$(cell_padded)_$local_cycle.csv", DataFrame)
    ndf.times = ndf.times .- idf.times[1]
    push!(df_vec, ndf)
end
end
idf.times = idf.times.-idf.times[1]
#filter!(row->row.Ns>=4,df)

vfull = 4.2
if cell==7
	vfull=4.0
elseif cell ==23
	vfull=4.1
end


initialcond = Dict("Starting Voltage[V]"=>vfull,"Ambient Temperature[K]" => idf.TemperatureC[1].+273.15)
current = -idf.ImA./1000


#Get and set up interpolants
cycle_array_vec = []
interpolated_vectors_vec = []
for (i, df) in enumerate(df_vec)
    current_interpolant = LinearInterpolation(-df.ImA./1000,df.times)
    voltage_interpolant = LinearInterpolation(df.EcellV,df.times)
    temperature_interpolant = LinearInterpolation(df.TemperatureC.+273.15,df.times)

    interpolated_time = collect(range(df.times[1],stop=df.times[end],step=1.0))


    interpolated_current = current_interpolant.(interpolated_time)
    interpolated_voltage = voltage_interpolant.(interpolated_time)
    interpolated_temperature = temperature_interpolant.(interpolated_time)

    interpolated_vectors = Dict(
        "Current" => interpolated_current,
        "Temperature" => interpolated_temperature,
        "Voltage" => interpolated_voltage,
        "Time" => interpolated_time,

    )
    push!(interpolated_vectors_vec, interpolated_vectors)


    #set up cycle arrays
    cycle_array = CellFitElectrolyte.current_profile(interpolated_current,interpolated_time)
    num_steps = Int(cycle_array[1])
    times = cycle_array[2:num_steps+1]
    types = cycle_array[num_steps+2:num_steps+1+num_steps]
    values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
    push!(cycle_array_vec, cycle_array) 
end



cellgeometry = CellFitElectrolyte.cell_geometry()

#ML Data
#Load data
FOLDERNAME = "results/outputs0106_fullcyc/"
vah = lpad(cell, 2, "0")
first_chain = load(FOLDERNAME*"VAH$(vah)_$(first_cycle)_HMC.jld2")["chain"]
last_chain = load(FOLDERNAME*"VAH$(vah)_$(last_cycle)_HMC.jld2")["chain"]

#Build Distributions
distribution_dict = Dict()

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

function evaluator(p::ComponentVector{T}, cycle_array_vec, interpolated_vectors_vec) where {T}
    # Handle Initial Conditions
    u::Array{T,1}  = Array{T, 1}(undef, input_dim)
    CellFitElectrolyte.initial_conditions!(u,p.p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
    Temp = idf.TemperatureC[1].+273.15

    u[8] = p.ics.ω
    u[9] = p.ics.εₑ⁻
    u[10] = p.ics.εₑ⁺
    u[11] = p.ics.frac_sol_am_neg
    u[12] = p.ics.frac_sol_am_pos
    
    if length(u) >= 13
        u[13:end] .= 0
    end

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

    #Create Function and Initialize Integrator
    function func(du, u, p, t)
        @unpack p_nn, p_phys = p
        #for (i, param) in enumerate(predicted_states)
        #    p_phys[param] = u[8+i]
        #end
        CellFitElectrolyte.equations_electrolyte_allocating((@view du[1:8]), (@view u[1:8]), p_phys, t, cache, cellgeometry, cathodeocv, anodeocv)
        
        xₛˢ⁻ = (u[1].-anodeocv.c_s_min)./(anodeocv.c_s_max-anodeocv.c_s_min)
        xₛᵇ⁻ = (u[2].-anodeocv.c_s_min)./(anodeocv.c_s_max-anodeocv.c_s_min)
        xₑ⁻ = u[3] ./ p_phys.cₑ₀
        xₑˢ = u[4] ./ p_phys.cₑ₀
        xₑ⁺ = u[5] ./ p_phys.cₑ₀
        xₛᵇ⁺ = (u[6].-anodeocv.c_s_min)./(anodeocv.c_s_max-anodeocv.c_s_min)
        xₛˢ⁺ = (u[7].-anodeocv.c_s_min)./(anodeocv.c_s_max-anodeocv.c_s_min)

        ω = (u[8] .- 0.01)./0.04
        εₑ⁻ = (u[9] .- 0.05)./0.45
        εₑ⁺ = (u[10] .- 0.05)./0.45
        frac_sol_am_neg = (u[11] .- 0.5)./0.5
        frac_sol_am_pos = (u[12] .- 0.5)./0.5
        if input_dim < 13
            u_in = SVector(xₛˢ⁻, xₛᵇ⁻, xₑ⁻, xₑˢ, xₑ⁺, xₛᵇ⁺, xₛˢ⁺, ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos)
        else
            u_in = SVector(xₛˢ⁻, xₛᵇ⁻, xₑ⁻, xₑˢ, xₑ⁺, xₛᵇ⁺, xₛˢ⁺, ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos, u[13:end]...)
        end
        new_u = nn(u_in, p_nn)
        du[8:end] .= new_u
        return nothing
    end
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,Rodas5(autodiff=false),save_everystep=false)
    integrator.opts.maxiters=1e7

    #we're really only interested in temperature and voltage, so we'll just save those
    total_length = sum([length(interpolated_vector["Voltage"]) for interpolated_vector in interpolated_vectors_vec])
    length_vec = vcat(0, cumsum([length(interpolated_vector["Voltage"]) for interpolated_vector in interpolated_vectors_vec]))

    endV::Array{T,1} = Array{T,1}(undef,total_length-1)
    endt::Array{T,1} = Array{T,1}(undef,total_length-1)

    for v in 1:length(cycle_array_vec)
        new_u = deepcopy(integrator.u)
        CellFitElectrolyte.initial_conditions!(new_u,p.p_phys,cellgeometry,initialcond,cathodeocv,anodeocv)
        set_u!(integrator, new_u)

        cycle_array = cycle_array_vec[v]
        interpolated_vectors = interpolated_vectors_vec[v]

        num_steps = Int(cycle_array[1])
        times = cycle_array[2:num_steps+1]
        for t in times[2:end]
            add_tstop!(integrator, t)
        end
        types = cycle_array[num_steps+2:num_steps+1+num_steps]
        values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]
        
        interpolated_temperature = interpolated_vectors["Temperature"]
        set_t!(integrator, times[1])

        for step::Int in 1:num_steps-1
            Temp = interpolated_temperature[step]
            @pack! p.p_phys = Temp
            input_type = types[step]
            input_value = values[step]
            end_time::T = times[step+1]
            @pack! p.p_phys = input_type,input_value
            while integrator.t < times[step+1]
                step!(integrator)
                if integrator.sol.retcode != :Default
                    return zeros(length(integrator.u))
                end
            end
            if any(integrator.u[1:7] .< 0)
                endV[step + length_vec[v]] = integrator.u[findfirst(integrator.u .< 0)]
                continue
            end
            cₛˢ⁺ = integrator.u[7]
            cₛˢ⁻ = integrator.u[1]
            x⁺ = (cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min)
            x⁻ = (cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min)
            if (x⁺ >= 1)
                endV[step + length_vec[v]] = 50*x⁺
                continue
            elseif (x⁻ >= 1)
                endV[step + length_vec[v]] = 50*x⁻
                continue
            end
            endV[step + length_vec[v]] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p.p_phys,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
            endt[step + length_vec[v]] = integrator.t
        end
    end
    return integrator.u
end

d = load("results/outputs0106_fullcyc/$(VAH)_HMC.jld2")
chain = d["chain"]

ω = mean(chain[:ω].data[:,1])
εₑ⁻ = mean(chain[:εₑ⁻].data[:,1])
εₑ⁺ = mean(chain[:εₑ⁺].data[:,1])
frac_sol_am_neg = mean(chain[:frac_sol_am_neg].data[:,1])
frac_sol_am_pos = mean(chain[:frac_sol_am_pos].data[:,1])


x⁻₀ = 0.6

εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

p_phys = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)

p_nn = initial_params(nn)

@model function fit_nn(cycle_array_vec, interpolated_vectors_vec, distribution_dict)
    p_nn ~ MvNormal(zeros(length(initial_params(nn))), 1e-9)

    ω = rand(distribution_dict[:ω]["Initial"])
    εₑ⁻ = rand(distribution_dict[:εₑ⁻]["Initial"])
    εₑ⁺ = rand(distribution_dict[:εₑ⁺]["Initial"])
    frac_sol_am_neg = rand(distribution_dict[:frac_sol_am_neg]["Initial"])
    frac_sol_am_pos = rand(distribution_dict[:frac_sol_am_pos]["Initial"])

    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p_phys = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)

    ics = ComponentVector(ω=ω, εₑ⁻=εₑ⁻, εₑ⁺=εₑ⁺, frac_sol_am_neg=frac_sol_am_neg, frac_sol_am_pos=frac_sol_am_pos)
    p = ComponentVector(p_nn = p_nn, p_phys = p_phys, ics=ics)
    u_out = evaluator(p, cycle_array_vec, interpolated_vectors_vec)

    ω_out, εₑ⁻_out, εₑ⁺_out, frac_sol_am_neg_out, frac_sol_am_pos_out  = u_out[8:12]


    probability_ω = Distributions.logpdf(distribution_dict[:ω]["Final"], ω_out)
    probability_εₑ⁻ = Distributions.logpdf(distribution_dict[:εₑ⁻]["Final"], εₑ⁻_out)
    probability_εₑ⁺ = Distributions.logpdf(distribution_dict[:εₑ⁺]["Final"], εₑ⁺_out)
    probability_frac_sol_am_neg = Distributions.logpdf(distribution_dict[:frac_sol_am_neg]["Final"], frac_sol_am_neg_out)
    probability_frac_sol_am_pos = Distributions.logpdf(distribution_dict[:frac_sol_am_pos]["Final"], frac_sol_am_pos_out)

    Turing.@addlogprob! probability_ω
    Turing.@addlogprob! probability_εₑ⁻
    Turing.@addlogprob! probability_εₑ⁺
    Turing.@addlogprob! probability_frac_sol_am_neg
    Turing.@addlogprob! probability_frac_sol_am_pos
    return nothing
end

model = fit_nn(cycle_array_vec, interpolated_vectors_vec, distribution_dict)

chain = sample(model, MH(), 1000)


