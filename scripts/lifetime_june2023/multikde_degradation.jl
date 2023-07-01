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

axes3d = PythonPlot.pyimport("mpl_toolkits.mplot3d.axes3d")

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

params = (:εₑ⁺, :εₑ⁻, :ω)

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

dims = [ContinuousDim(), ContinuousDim(), ContinuousDim()]
bws = [chainoperator(first_chain, params[1], scottoperator), chainoperator(first_chain, params[2], scottoperator), chainoperator(first_chain, params[3], scottoperator)]


N = 1000
observations = Vector{Vector{Float64}}(undef, N)
@showprogress "getting observations" for i in 1:N
    observations[i] = [first_chain[params[1]].data[i,1], first_chain[params[2]].data[i,1], first_chain[params[3]].data[i,1]]
end
    

multi_kde = KDEMulti(dims, bws, observations)



min_x_1 = chainoperator(first_chain, params[1], minimum)
max_x_1 = chainoperator(first_chain, params[1], maximum)
mean_x_1 = chainoperator(first_chain, params[1], mean)
range_x_1 = max_x_1 - min_x_1


min_x_2 = chainoperator(first_chain, params[2], minimum)
max_x_2 = chainoperator(first_chain, params[2], maximum)
mean_x_2 = chainoperator(first_chain, params[2], mean)
range_x_2 = max_x_2 - min_x_2

min_x_3 = chainoperator(first_chain, params[3], minimum)
max_x_3 = chainoperator(first_chain, params[3], maximum)
mean_x_3 = chainoperator(first_chain, params[3], mean)
range_x_3 = max_x_3 - min_x_3

x_1r = range(min_x_1 - range_x_1/5, stop=max_x_1 + range_x_1/5, length=25)
x_2r = range(min_x_2 - range_x_2/5, stop=max_x_2 + range_x_2/5, length=25)
x_3r = range(min_x_3 - range_x_3/5, stop=max_x_3 + range_x_3/5, length=25)

x_g = [[x_1, x_2, x_3] for x_1 in x_1r for x_2 in x_2r for x_3 in x_3r]
z = zeros(length(x_g))
@showprogress "computing pd on grid" for (i,x) in enumerate(x_g)
    z[i] = MultiKDE.pdf(multi_kde, x)
end
x_1_g = [x[1] for x in x_g]
x_2_g = [x[2] for x in x_g]
x_3_g = [x[3] for x in x_g]

m = z.>1
#=
figure(1)
clf()
scatter3D(x_1_g[m], x_2_g[m], x_3_g[m], c=z[m]./maximum(z))
=#

df = DataFrame(x = x_1_g, y=x_2_g, z=x_3_g, p=z)