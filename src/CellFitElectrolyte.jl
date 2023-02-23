module CellFitElectrolyte

const R = 8.314
const F = 96485.3365

using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEq
using Parameters
using StaticArrays
using OCV
using FileIO
using JLD2
using DataInterpolations
using PreallocationTools
using DiffEqFlux
using Turing
using DynamicHMC

include("parameters.jl")
include("cache.jl")
include("models.jl")
#include("simulator.jl")
include("initialize.jl")
include("inputbuilder.jl")
include("SimulatedAnnealing.jl")
include("models_allocating.jl")
include("degradation_nn.jl")
include("models_concentration_dependent.jl")
include("models_allocating_new.jl")

export simulate,initialize
end # module
