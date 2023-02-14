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
using ProgressMeter
using PythonPlot

cache = CellFitElectrolyte.initialize_cache(Float64)
num_rows = 6
num_cols = 4
fig, axes = subplots(6,4,figsize=(8,10))

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

filenames = readdir("results/outputs0117_elec")
vahs = similar(filenames)
for (i,file) in enumerate(filenames)
    vah = split(file,"_")[1]
    vahs[i] = vah
end
cells = unique(vahs)

for (k,cell) in enumerate(cells)
    println(cell)
    row = Int(floor((k-1)/num_cols)) + 1
    col = k%num_cols
    if col == 0
        col = 4
    end
cell_error_vec = Float64[]
cycle_vec = Int[]
strings = []
for file in readdir("results/errors_0117")
    inner_cell = split(file,"_")[1]
    if inner_cell != cell
        continue
    end
    f = split(file, "_errors.jld2")[1]
    try
        test1234234 = load("results/outputs0117_elec/$(f)_HMC.jld2")
        d = load("results/errors_0117/$file")
        rmse = d["rmse"]
        append!(cell_error_vec, rmse)
        cycle = parse(Int,split(file, "_")[2])
        append!(cycle_vec, cycle*ones(1000))
        for n in 1:1000
            push!(strings, string(f))
        end
    catch
        continue
    end
end

sorted_cycles_by_error = cycle_vec[sortperm(cell_error_vec)]
med = sorted_cycles_by_error[Int(floor(length(sorted_cycles_by_error)))]

VAH = "$(cell)_$(med)"
split1 = split(VAH,['H','_'])
cycle = parse(Int,split1[3])


df = CSV.read("/Users/abills/Datasets/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
#filter!(row->row.Ns>=4,df)

d = load("results/outputs0117_elec/$(VAH)_HMC.jld2")
chain = d["chain"]

vfull = 4.2
if cell==7
	vfull=4.0
elseif cell ==23
	vfull=4.1
end


initialcond = Dict("Starting Voltage[V]"=>vfull,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
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
adsf = parse(Int, cell[end-1:end])
other_cycle_array = CellFitElectrolyte.load_airbus_cyclearrays()["cycle_array_vector"][adsf][cycle]

new_num_steps = Int(other_cycle_array[1])
new_times = other_cycle_array[2:new_num_steps+1]
new_types = other_cycle_array[new_num_steps + 2:new_num_steps+1+new_num_steps]
new_values = other_cycle_array[new_num_steps+2+new_num_steps:new_num_steps+1+2*new_num_steps]


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
        
        idx = searchsortedlast(new_times, integrator.t)
        type = new_types[idx]
        value = new_values[idx]
        #calculate voltage on first step
        if step == 1
            if any(integrator.u .< 0)
                endV[step] = integrator.u[findfirst(integrator.u .< 0)]
                continue
            end
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
            Voltage = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        else
            Voltage = endV[step-1]
        end
    
        Temp = interpolated_temperature[step]
        @pack! p = Temp
        input_type = types[step]
        input_value = values[step]
        end_time::T = times[step+1]
        
        if type==0
            input_value = 0
        elseif type==1
            input_value = value/Voltage
        elseif type==2
            error("can't do this")
        elseif type==3
            input_value = value;
        elseif type==4
            input_value = value
        elseif type==5
            value = input_value
        end
        

        @pack! p = input_type,input_value
        while integrator.t < end_time
            step!(integrator)
        end
        if any(integrator.u .< 0)
            endV[step] = integrator.u[findfirst(integrator.u .< 0)]
            continue
        end
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
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        endt[step] = integrator.t
    end
    return endV
end

function chain_to_nt(chain,idx)
    return (ω = chain[:ω].data[idx,1], εₑ⁺=chain[:εₑ⁺].data[idx,1],εₑ⁻ = chain[:εₑ⁻].data[idx,1], frac_sol_am_pos = chain[:frac_sol_am_pos].data[idx,1],frac_sol_am_neg=chain[:frac_sol_am_neg].data[idx,1])
end

function fit_cfe(interpolated_voltage, nt)
    # Prior distributions.
    #n_li ~ truncated(Normal(0.2, 0.01), 0.16, 0.22)
    @unpack ω, εₑ⁺, εₑ⁻, frac_sol_am_pos, frac_sol_am_neg = nt

    x⁻₀ = 0.6

    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)
    
    
    predicted = evaluator(p)
    return predicted
end
PythonPlot.matplotlib.rcParams["font.size"] = 10
N = length(chain)

rmse = zeros(N)
i=1
c = chain_to_nt(chain,i)
predicted = fit_cfe(interpolated_voltage,c)
axes[row-1,col-1].plot(interpolated_time[1:end-1], predicted,color="xkcd:grey")
axes[row-1,col-1].plot(interpolated_time[1:end-1], predicted,color="xkcd:grey",label="Model")
@showprogress for i =2:N
    c = chain_to_nt(chain,i)
    predicted = fit_cfe(interpolated_voltage,c)
    axes[row-1,col-1].plot(interpolated_time[1:end-1], predicted,color="xkcd:grey")
end
axes[row-1,col-1].plot(interpolated_time, interpolated_voltage,color="tab:blue", label="Experiment")
axes[row-1,col-1].legend()
axes[row-1,col-1].set_xlabel("Time [sec]", fontsize=18)
axes[row-1,col-1].set_ylabel("Cell voltage [V]", fontsize=18)
axes[row-1,col-1].grid(alpha=0.5)
end
fig.tight_layout()
fig.savefig("figs/si/voltage_medians.png",bbox_inches="tight")
