using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using BenchmarkTools
using Test

cache = CellFitElectrolyte.initialize_cache(Float64)
p = CellFitElectrolyte.p_transport()
CellFitElectrolyte.fill_transport!(cache.A,p.θₛ⁻,p.θₑ,p.θₛ⁺)


cellgeometry = CellFitElectrolyte.cell_geometry()
du = similar(u)
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
u = Array{Float64,1}(undef,14)
initialcond = Dict("Starting Voltage[V]"=>4.0,"Ambient Temperature[K]"=>298.0)
CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,anodeocv,cathodeocv)





benchmark = @benchmark CellFitElectrolyte.equations_electrolyte($du,$u,$p,1.0,$cache,$cellgeometry,$cathodeocv,$anodeocv,3,4.2)

@test benchmark.allocs<=0
@test benchmark.memory<=0
display(benchmark)
