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

cache = CellFitElectrolyte.initialize_cache(Float64)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

CELLS = ["VAH01", "VAH01", "VAH01"]
CYCLES = [[2, 200, 400, 600, 845], [2, 200, 400, 600, 845], [2, 200, 400, 600, 845]]

figure(1)
clf()
figure(2)
clf()


for (i,c) in enumerate(CELLS)
    for (j, cycle) in enumerate(CYCLES[i])
    VAH = c*"_$cycle"
    split1 = split(VAH,['H','_'])
    cell = parse(Int,split1[2])
    cycle = parse(Int,split1[3])


df = CSV.read("/Users/abills/Datasets/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
filter!(row->row.Ns>=4,df)

vfull = 4.2
if cell==7
	vfull=4.0
elseif cell ==23
	vfull=4.1
end


initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => df.TemperatureC[1].+273.15)
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
        Temp = interpolated_temperature[step]
        @pack! p = Temp
        input_type = types[step]
        input_value = values[step]
        end_time::T = times[step+1]
        @pack! p = input_type,input_value
        while integrator.t < times[step+1]
            step!(integrator)
        end
        if any(integrator.u .< 0)
            endV[step] = integrator.u[findfirst(integrator.u .< 0)]
            println("shid")
            continue
        end
        cₛˢ⁺ = integrator.u[7]
        cₛˢ⁻ = integrator.u[1]
        x⁺ = (cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min)
        x⁻ = (cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min)
        if (x⁺ >= 1)
            endV[step] = 50*x⁺
            println("shid")
            continue
        elseif (x⁻ >= 1)
            endV[step] = 50*x⁻
            println("shid")
            continue
        end
        endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        endt[step] = integrator.t
    end
    return endV
end

thing = []

chain = load("results/outputs0114_elec/$(VAH)_HMC.jld2")["chain"]
params = DataFrame(chain)
for n in 1:size(params)[1]
    frac_sol_am_neg = params.frac_sol_am_neg[n]
    frac_sol_am_pos = params.frac_sol_am_pos[n]
    εₑ⁺ = params.εₑ⁺[n]
    εₑ⁻ = params.εₑ⁻[n]
    x⁻₀ = 0.6
    ω = params.ω[n]

    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)

    V = evaluator(p)
    figure(1)
    subplot(length(CELLS), length(CYCLES[1]), length(CYCLES[1])*(i-1) + j)
    plot(interpolated_time[1:end-1], V, "xkcd:grey", alpha=0.1)
    err = sqrt(mean((interpolated_voltage[1:end-1] .- V).^2))
    append!(thing, err)
end
figure(1)
subplot(length(CELLS), length(CYCLES[1]), length(CYCLES[1])*(i-1) + j)
plot(interpolated_time, interpolated_voltage)
figure(2)
subplot(length(CELLS), length(CYCLES[1]), length(CYCLES[1])*(i-1) + j)
hist(thing)
end
end

figure(1)
savefig("voltage_trajectories.pdf")
figure(2)
savefig("error_histograms.pdf")




#=
εₛ⁻ = kde(chain[:εₛ⁻].data[:,1])
εₛ⁺ = kde(chain[:εₛ⁺].data[:,1])
εᵧ⁺ = kde(chain[:εᵧ⁺].data[:,1])
εᵧ⁻ = kde(chain[:εᵧ⁻].data[:,1])
ω = kde(chain[:ω].data[:,1])
n_li = kde(chain[:n_li].data[:,1])



figure(1);
clf()
subplot(321)
plot(εₛ⁻.x, εₛ⁻.density)
xlabel("εₛ⁻")
ylabel("Density")
grid()
subplot(322)
plot(ω.x, ω.density)
ylabel("Density")
xlabel("SEI Resistance")
grid()
subplot(323)
plot(n_li.x,n_li.density)
ylabel("Density")
xlabel("Moles Li")
grid()
subplot(324)
plot(εₛ⁺.x,εₛ⁺.density)
ylabel("Density")
xlabel("εₛ⁺")
grid()
subplot(325)
plot(εᵧ⁺.x,εᵧ⁺.density)
ylabel("Density")
xlabel("εᵧ⁺")
grid()
subplot(326)
plot(εᵧ⁻.x,εᵧ⁻.density)
ylabel("Density")
xlabel("εᵧ⁻")
grid()
=#
