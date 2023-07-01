using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.Parameters
using CellFitElectrolyte.DataInterpolations
using CellFitElectrolyte.DiffEqFlux
using CellFitElectrolyte.Turing
using BenchmarkTools
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
using Turing
using PythonPlot
using MultiKDE

#set up simulation
cache = CellFitElectrolyte.initialize_cache(Float64)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()
initialcond = Dict("Starting Voltage[V]"=>4.2,"Ambient Temperature[K]" => 300.0)
vfull = initialcond["Starting Voltage[V]"]
cellgeometry = CellFitElectrolyte.cell_geometry()

#Problem is now built. Next step is to load data. This section roughly corresponds to batch settings.
cell = 1
first_cycle = 500
last_cycle = 100

params = (:ω, :εₑ⁻, :εₑ⁻, :frac_sol_am_neg, :frac_sol_am_pos)

function chainoperator(chain, param, operator)
    return operator(chain[param].data[:,1])
end

function scottoperator(data)
    σ = std(data)
    n = length(data)
    return 1.06*σ*(n^(-1/5))
end


#Load data
FOLDERNAME = "results/0302/"
vah = lpad(cell, 2, "0")
first_chain = load(FOLDERNAME*"VAH$(vah)_$(first_cycle)_HMC.jld2")["chain"]
last_chain = load(FOLDERNAME*"VAH$(vah)_$(last_cycle)_HMC.jld2")["chain"]

dims = [ContinuousDim() for param in params]
bws = [chainoperator(first_chain, param, scottoperator) for param in params]


N = 1000
observations = Vector{Vector{Float64}}(undef, N)
@showprogress "getting observations" for i in 1:N
    observations[i] = [first_chain[params[j]].data[i,1] for j in 1:length(params)]
end
    

multi_kde = KDEMulti(dims, bws, observations)



