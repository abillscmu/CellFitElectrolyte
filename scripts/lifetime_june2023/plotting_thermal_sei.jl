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
using Optim
using MultiKDE

@load "thermal_sei.jld2"


#cell = "VAH09"
fig,ax = subplots(2,Int(length(keys(distribution_dict))/2))
for (i,cell) in enumerate(keys(distribution_dict))
    i_idx = (i+1)%3
    if i > length(keys(distribution_dict))/2
        j_idx = 1
    else
        j_idx = 0
    end
    old_i = i_idx
    i_idx = j_idx
    j_idx = old_i
    ax[i_idx,j_idx].plot(vcat(1, fitting_cycles[cell]), my_sol[cell][:frac_sol_am_neg][1])
    xs = collect(keys(distribution_dict[cell]))
    ys = zeros(1000,length(xs))
    for (i,x) in enumerate(xs)
        ys[:,i] .= distribution_dict[cell][x][:frac_sol_am_neg]
    end
    ax[i_idx,j_idx].violinplot(ys,positions=xs,widths=maximum(xs)/10)
    ax[i_idx,j_idx].set_xlabel("Cycle number")
    ax[i_idx,j_idx].set_ylabel("fₛ⁻")
    ax[i_idx,j_idx].set_title(cell)
end
fig.tight_layout()
fig.savefig("test.pdf",bbox_inches="tight")