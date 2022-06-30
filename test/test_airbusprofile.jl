using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CSV
using DataFrames
using Test
using ProgressMeter

type = Float64

cache = CellFitElectrolyte.initialize_cache(type)
p = CellFitElectrolyte.p_transport()
CellFitElectrolyte.fill_transport!(cache.A,p.θₛ⁻,p.θₑ,p.θₛ⁺)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()

cycle_array_vec = CellFitElectrolyte.load_airbus_cyclearrays()
cycle_array = cycle_array_vec["cycle_array_vector"][1][20]
df = CSV.read("data/VAH01_20.csv",DataFrame)
df.times = df.times.-df.times[1]
current = -df.ImA./1000
p.Tamb = df.TemperatureC[1].+273.15
cycle_array = CellFitElectrolyte.current_profile(current,df.times)

initialcond = Dict(
    "Starting Voltage[V]"=>4.2,
    "Ambient Temperature[K]"=>298.0,
    "Capacity[mAh]"=>3000.0)

cellgeometry = CellFitElectrolyte.cell_geometry()


u  = Array{type,1}(undef,14)
u = CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
u[13] = p.Tamb
du = similar(u)

mass_mat = Matrix{type}(I,14,14)
mass_mat[14,14] = 0.0



num_steps = Int(cycle_array[1])
times = @view cycle_array[2:num_steps+1]
types = @view cycle_array[num_steps+2:num_steps+1+num_steps]
values = @view cycle_array[num_steps+2+num_steps:num_steps+1+2*num_steps]

func = ODEFunction((du,u,p,t)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv),mass_matrix=mass_mat)
prob = ODEProblem(func,u,(0.0,times[end]),p)

integrator = init(prob,Rodas4(autodiff=false),dtmax=0.1,save_everystep=false)

@showprogress "simulating..." for step in 1:num_steps-1
    input_type = types[step]
    input_value = values[step]
    if input_type==5
        input_type = 3
        input_value = values[step]
        @pack! p = input_type,input_value
        V = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
        while V<=initialcond["Starting Voltage[V]"]
            step!(integrator)
            V = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
        end
        input_type = 2
        input_value = initialcond["Starting Voltage[V]"]
        @pack! p = input_type,input_value
        while integrator.u[14]<=-0.1
            step!(integrator)
        end
        input_type = 0
        input_value = 0
        @pack! p = input_type,input_value
        while integrator.t<=times[step+1]
            step!(integrator)
        end
    else
        end_time = times[step+1]
        @pack! p = input_type,input_value
        step!(integrator,end_time-integrator.t,true)
    end
    savevalues!(integrator,true)
end

sol = integrator.sol

concentrations = sol[1:7,:]
Vₛ⁻ = cellgeometry.Vₛ⁻
Vₛ⁺ = cellgeometry.Vₛ⁺
Vₑ⁻ = cellgeometry.Vₑ⁻
Vₑˢ = cellgeometry.Vₑˢ
Vₑ⁺ = cellgeometry.Vₑ⁺

εₛ⁻ = sol[8,1]
εₑ⁻ = sol[9,1]
εₑˢ = sol[10,1]
εₑ⁺ = sol[11,1]
εₛ⁺ = sol[12,1]

mass_vector = [Vₛ⁻*εₛ⁻,Vₛ⁻*εₛ⁻,Vₑ⁻*εₑ⁻,Vₑˢ*εₑˢ,Vₑ⁺*εₑ⁺,Vₛ⁺*εₛ⁺,Vₛ⁺*εₛ⁺]
n_li = concentrations'*mass_vector
n_li_change = maximum(abs.((n_li.-(n_li[1]))./n_li[1]))

V = CellFitElectrolyte.calc_voltage(sol,p,sol.t,cache,cellgeometry,cathodeocv,anodeocv)


x_s_n = sol[1,:]./anodeocv.c_s_max
x_b_n = sol[2,:]./anodeocv.c_s_max

x_s_p = sol[6,:]./cathodeocv.c_s_max
x_b_p = sol[7,:]./cathodeocv.c_s_max

myOCV = cathodeocv.(x_s_p,sol[13,:]).-anodeocv.(x_s_n,sol[13,:])

V_rmse = sqrt.(mean((V.-df.EcellV).^2))
println("Voltage RMSE is $V_rmse")

using MATLABPlots
figure(1)
clf()
subplot(3,1,1)
plot(sol.t,V,options=Dict("LineWidth"=>3))
hold_on()
plot(df.times,df.EcellV,options=Dict("LineWidth"=>3))
setgca(Dict(
    "FontSize"=>24,
    "FontName"=>"Open Sans",
))
#xlabel("Time[s]")
ylabel("Voltage[V]")
grid_on()


subplot(3,1,2)
plot(sol.t,sol[13,:],options=Dict("LineWidth"=>3))
hold_on()
plot(df.times,df.TemperatureC.+273.15,options=Dict("LineWidth"=>3))
setgca(Dict(
    "FontSize"=>24,
    "FontName"=>"Open Sans",
))
#xlabel("Time[s]")
ylabel("Temperature[K]")
grid_on()


subplot(3,1,3)
plot(sol.t,sol[1,:],options=Dict("LineWidth"=>3))
hold_on()
plot(sol.t,sol[2,:],options=Dict("LineWidth"=>3))
setgca(Dict(
    "FontSize"=>24,
    "FontName"=>"Open Sans",
))
xlabel("Time[s]")
ylabel("c^{-}[mol/m^3]")

@test n_li_change<1e-10




