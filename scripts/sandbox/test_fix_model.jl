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

VAH = "VAH01_20"
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("data/$(VAH).csv",DataFrame)
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



cellgeometry = CellFitElectrolyte.cell_geometry()
function evaluator(p::ComponentVector{T}) where {T}
    # Handle Initial Conditions
    u::Array{T,1}  = Array{T,1}(undef,7)
    CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    #return u
    Temp = df.TemperatureC[1].+273.15
    @pack! p = Temp
    du = similar(u)
    input_type::T = types[1]
    input_value::T = values[1]
    @pack! p = input_type,input_value

    #Create Function and Initialize Integrator
    func = ODEFunction((du, u, p, t)->CellFitElectrolyte.equations_electrolyte_allocating_new(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
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





function fit_cfe(vec)

    θₛⁿ = vec[1]
    θₛᵖ = vec[2]
    θₑ = vec[3]
    # Prior distributions.
    #n_li ~ truncated(Normal(0.2, 0.01), 0.16, 0.22)
    #ω ~ truncated(Normal(0.02, 0.001), 0.01, 0.05)
    ω = vec[4]
    x⁻₀ = 0.6
    
    #coordinate transforms to stay in a nice area
    #εₑ⁺ ~ Uniform(0.05, 0.5)
    εₑ⁺ = vec[5]
    #εₑ⁻ ~ Uniform(0.05, 0.5)
    εₑ⁻ = vec[6]
    #frac_sol_am_pos ~ Uniform(0.5, 1.0)
    #frac_sol_am_neg ~ Uniform(0.5, 1.0)
    frac_sol_am_pos = vec[7]
    frac_sol_am_neg = vec[8]

    κ = vec[9]


    εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
    εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
    εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
    εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

    p = ComponentVector(θₛ⁻ = θₛⁿ, θₑ = θₑ, θₛ⁺ = θₛᵖ, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 1000., κ = κ, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)
    
    
    predicted = evaluator(p)

    # Observations.
    #interpolated_voltage[1:end-1] ~ MvNormal(predicted, 0.1)

    return predicted
end


predicted = fit_cfe(interpolated_voltage)



[9.083011990824314e-8, 4.882065444878806e-6, 0.0005461023009076725, 0.03433620706589757, 0.9405241587940202, 0.9270719936285006, 8.806917876188582, 8.601714488599779, 0.1346959901335123]
parent = [2.6145309903155526e-10, 2.4126123473971843e-7, 0.00018203759659602146, 0.027302814250432884, 0.30579391025803365, 0.07328409664905508, 0.7695552367778169, 0.8598952621868902, 0.8464583490637971]
ub = similar(parent)
lb = similar(parent)
for (i,p) in enumerate(parent)
    ub[i] = 10*p
    lb[i] = 0.1*p
end

ub[5:8] .= 1
lb[5:8] .= 0.001

#CellFitElectrolyte.anneal(fit_cfe, parent, ub, lb)

# Sample 3 independent chains with forward-mode automatic differentiation (the default).
#chain = sample(model, NUTS(0.65), MCMCSerial(), 1000, 1; progress=true)
#d = Dict("chain" => chain)


#save("$(VAH)_HMC.jld2", d)
#sleep(5)


predicted = fit_cfe(parent)
fig, axes = subplots()
axes.plot(interpolated_time[1:end-1], predicted)
axes.plot(interpolated_time, interpolated_voltage)