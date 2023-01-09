using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.DiffEqFlux
using CellFitElectrolyte.Turing
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
num_augmented_states = 2
predicted_states = [:ω]

#set up neural network (derived)
num_predicted_states = length(predicted_states)
input_dim = 8 + num_predicted_states + num_augmented_states
output_dim = num_predicted_states + num_augmented_states
nn = CellFitElectrolyte.build_nn(num_layers, activation, width, input_dim, output_dim)

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

p = ComponentVector{Any}(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.6, εₛ⁺ = 0.75, εᵧ⁺ = 0, εᵧ⁻ = 0, c = 50.0, h = Inf, Tamb = 320.0, Temp = 320.0, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = 0.01, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch=false, cccv_switch_2=false, vfull=vfull, ifull=-0.01)

p_nn_init = DiffEqFlux.initial_params(nn)

u = zeros(input_dim)
du = similar(u)

CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)

function f(du, u, p_nn, t)
    for (i, param) in enumerate(predicted_states)
        p[param] = u[8+i]
    end
    CellFitElectrolyte.equations_electrolyte_life_allocating((@view du[1:8]), (@view u[1:8]), p, t, cache, cellgeometry, cathodeocv, anodeocv)
    new_u = nn(u, p_nn)
    du[9:end] = new_u
    return nothing
end

input_type = 3
input_value = 4.2
@pack! p = input_type, input_value

vars = ones(length(u))
vars[8] = 0
mm = diagm(vars)

func = ODEFunction(f, mass_matrix=mm)
prob = ODEProblem(func,u,(0.0, 1.0),p_nn_init)

@model function fit_nn(prob)
    p_nn ~ MvNormal(zeros(length(p_nn_init)), 0.0001)
    sol = solve(prob, QNDF(autodiff=false), p = p_nn)
    ω = sol[9,end]
    ω ~ Normal(1.0, 0.000001)
end

model = fit_nn(prob)

chain = sample(model, NUTS(0.65), 100)
