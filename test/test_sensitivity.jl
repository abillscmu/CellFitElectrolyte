using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using ForwardDiff
using Test


function loss(p)
type = eltype(p)
cache = CellFitElectrolyte.initialize_cache(type)
CellFitElectrolyte.fill_transport!(cache.A,p.θₛ⁻,p.θₑ,p.θₛ⁺)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()

initialcond = Dict(
    "Starting Voltage[V]"=>4.2,
    "Ambient Temperature[K]"=>298.0,
    "Capacity[mAh]"=>3000.0)

cellgeometry = CellFitElectrolyte.cell_geometry()


u  = Array{type,1}(undef,13)

u = CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
du = similar(u)
end_dt = 3400.0

prob = ODEProblem((du,u,p,t)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv),u,(0.0,end_dt),p,saveat=1.0)

sol = solve(prob,Rodas4(autodiff=false))

concentrations = sol[1:7,:]
Vₛ⁻ = cellgeometry.Vₛ⁻
Vₛ⁺ = cellgeometry.Vₛ⁺
Vₑ⁻ = cellgeometry.Vₑ⁻
Vₑˢ = cellgeometry.Vₑˢ
Vₑ⁺ = cellgeometry.Vₑ⁺

εₛ⁻ = sol[8,1]
εₑ⁻ = sol[9,1]
εₑˢ = sol[10,1]
εₑ⁺ = sol[11,1]
εₛ⁺ = sol[12,1]

mass_vector = [Vₛ⁻*εₛ⁻,Vₛ⁻*εₛ⁻,Vₑ⁻*εₑ⁻,Vₑˢ*εₑˢ,Vₑ⁺*εₑ⁺,Vₛ⁺*εₛ⁺,Vₛ⁺*εₛ⁺]
n_li = concentrations'*mass_vector
n_li_change = maximum(abs.((n_li.-(n_li[1]))./n_li[1]))

V = CellFitElectrolyte.calc_voltage(sol,p,sol.t,cache,cellgeometry,cathodeocv,anodeocv)

V_ref = ones(3401)
sol = sum(abs.(V.-V_ref))
return sol

end

