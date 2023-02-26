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

VAH = "VAH01_2"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("data/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]
#filter!(row->row.Ns>=4,df)

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
other_cycle_array = CellFitElectrolyte.load_airbus_cyclearrays()["cycle_array_vector"][cell][cycle]

new_num_steps = Int(other_cycle_array[1])
new_times = other_cycle_array[2:new_num_steps+1]
new_types = other_cycle_array[new_num_steps + 2:new_num_steps+1+new_num_steps]
new_values = other_cycle_array[new_num_steps+2+new_num_steps:new_num_steps+1+2*new_num_steps]


num_steps = Int(cycle_array[1])
times = cycle_array[2:num_steps+1]
types = cycle_array[num_steps+2:num_steps+1+num_steps]
values = cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]



cellgeometry = CellFitElectrolyte.cell_geometry_new()
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
    func = ODEFunction((du, u, p, t)->CellFitElectrolyte.equations_electrolyte_allocating_new(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
    prob = ODEProblem(func,u,(0.0,times[end]),p)
    integrator = init(prob,QNDF(autodiff=false),save_everystep=false, tstops = times, verbose=false)

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
            Voltage = CellFitElectrolyte.calc_voltage_new(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
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
            if integrator.sol.retcode != :Default
                println("hi")
                endV[step:end] .= 50
                return endV
            end
        end
        if any(integrator.u .< 0)
            #println(integrator.u)
            endV[step:end] .= 50*integrator.u[findfirst(integrator.u .< 0)]
            return endV
            #continue
        end
        cₛˢ⁺ = integrator.u[7]
        cₛˢ⁻ = integrator.u[1]
        x⁺ = (cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min)
        x⁻ = (cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min)
        if (x⁺ >= 1)
            #println("yo")
            endV[step] = 50*x⁺
            continue
        elseif (x⁻ >= 1)
            #println("dog")
            endV[step] = 50*x⁻
            continue
        end
        endV[step] = CellFitElectrolyte.calc_voltage_new(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv,values[step])
        endt[step] = integrator.t
    end
    return endV
end





@model function fit_cfe(interpolated_voltage)
    # Prior distributions.
    #n_li ~ truncated(Normal(0.2, 0.01), 0.16, 0.22)
    ω ~ truncated(Normal(0.02, 0.001), 0.0, 0.05)
    x⁻₀ = 0.6

    c_e_0 ~ truncated(Normal(1000, 500),1000,4000)
    
    #coordinate transforms to stay in a nice area
    εₑ⁺ ~ Uniform(0.001, 0.75)
    εₑ⁻ ~ Uniform(0.001, 0.75)

    frac_sol_am_pos ~ Uniform(0.1, 1.0)
    frac_sol_am_neg ~ Uniform(0.1, 1.0)

    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p = ComponentVector(θₛ⁻ = 9.130391775012598e-10, θₑ = 1e-9, θₛ⁺ = 3.728671559985511e-8, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 1e-1, k₀⁻ = 1e-1, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = c_e_0, κ = 0.01, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)
    
    
    predicted = evaluator(p)

    # Observations.
    interpolated_voltage[1:end-1] ~ MvNormal(predicted, 0.05)

    return nothing
end


model = fit_cfe(interpolated_voltage)

θ = (
    ω = 0.02,
    c_e_0 = 1000,
    εₑ⁺ = 0.2,
    εₑ⁻ = 0.2,
    frac_sol_am_pos = 0.75,
    frac_sol_am_neg = 0.75,
)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model, PG(10), MCMCSerial(), 100, 3; progress=true, θ = θ)
d = Dict("chain" => chain)

save("$(VAH)_HMC.jld2", d)
sleep(5)

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
