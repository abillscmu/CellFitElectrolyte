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

include("parameters.jl")
include("cache.jl")
include("models.jl")
#include("simulator.jl")
include("initialize.jl")
include("inputbuilder.jl")
include("SimulatedAnnealing.jl")

export simulate,initialize
end # module
