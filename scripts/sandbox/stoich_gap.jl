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

initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 25 .+273.15)



cellgeometry = CellFitElectrolyte.cell_geometry()

ω = 0.02


#coordinate transforms to stay in a nice area
εₑ⁺ = 0.2
εₑ⁻ = 0.2

frac_sol_am_pos = 1.0
frac_sol_am_neg = 1.0

εₛ⁻ = (1 - εₑ⁻)*frac_sol_am_neg
εₛ⁺ = (1 - εₑ⁺)*frac_sol_am_pos
εᵧ⁻ = 1 - εₛ⁻ - εₑ⁻
εᵧ⁺ = 1 - εₛ⁺ - εₑ⁺

x⁻₀ = 0.4

val = []
x_vec = 0.4:0.01:0.7

for x⁻₀ in x_vec

p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = εₛ⁻, εₛ⁺ = εₛ⁺, εᵧ⁺ = εᵧ⁺, εᵧ⁻ = εᵧ⁻, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = x⁻₀, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0)

u = zeros(7)
CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)


moles_cyclable_anode = 2*cellgeometry.Vₛ⁻ * εₛ⁻ * u[1]
moles_cyclable_cathode = 2*(cathodeocv.c_s_max - u[6])*cellgeometry.Vₛ⁺ * εₛ⁺
q_pos = u[6] ./ cathodeocv.c_s_max

moles_cyclable = min(moles_cyclable_anode, moles_cyclable_cathode)


append!(val, moles_cyclable)

end

figure(1)
clf()
scatter(x_vec, val)