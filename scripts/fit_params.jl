using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CSV
using DataFrames
using Test
using ProgressMeter
using LinearAlgebra
using Statistics

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

function evaluator(p)
T = eltype(p)
endV = Array{T,1}(undef,length(df.EcellV))
prob = ODEProblem(func,u,(0.0,times[end]),p)

integrator = init(prob,Rodas4(autodiff=false),save_everystep=false)

for step in 1:num_steps-1
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
    endV[step] = CellFitElectrolyte.calc_voltage(integrator.u,integrator.p,integrator.t,cache,cellgeometry,cathodeocv,anodeocv)
end

V_rmse::T = sqrt.(mean((endV.-df.EcellV).^2))

return V_rmse
end






