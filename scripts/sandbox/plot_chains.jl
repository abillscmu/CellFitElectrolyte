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
using PythonPlot
using StatsPlots

using Random
Random.seed!(14)

cache = CellFitElectrolyte.initialize_cache(Float64)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
p = CellFitElectrolyte.p_transport()

VAHs = ["VAH01_2", "VAH01_200", "VAH01_400", "VAH01_601", "VAH01_845"]

fig, axes = subplots(1,5, figsize=(8,2))
for (j,VAH) in enumerate(VAHs)
split1 = split(VAH,['H','_'])
cell = parse(Int,split1[2])
cycle = parse(Int,split1[3])


df = CSV.read("/Users/abills/Datasets/cycle_individual_data/$(VAH).csv",DataFrame)
df.times = df.times.-df.times[1]

d = load("results/0302/$(VAH)_HMC.jld2")
chain = d["chain"]
StatsPlots.plot(chain, left_margin=4StatsPlots.Plots.mm, title=VAH)
savefig("figs/si/convergence_$VAH.pdf")
println(VAH)
println(summarize(chain).nt.rhat)

end

