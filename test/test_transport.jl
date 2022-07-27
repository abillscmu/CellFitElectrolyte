using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using CellFitElectrolyte.LinearAlgebra
using CellFitElectrolyte.Parameters
using Test

type = Float64

cache = CellFitElectrolyte.initialize_cache(type)
p = CellFitElectrolyte.p_transport()
CellFitElectrolyte.fill_transport!(cache.A,p.θₛ⁻,p.θₑ,p.θₛ⁺)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()

initialcond = Dict(
    "Starting Voltage[V]"=>4.2,
    "Ambient Temperature[K]"=>298.0,
    "Capacity[mAh]"=>3000.0)

cellgeometry = CellFitElectrolyte.cell_geometry()


u  = Array{type,1}(undef,7)

u = CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
du = similar(u)
end_dt = 2500.0

#mass_mat = Matrix{type}(I,14,14)
#mass_mat[14,14] = 0.0

func = ODEFunction((du,u,p,t)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv))
prob = ODEProblem(func,u,(0.0,end_dt),p)

sol = solve(prob,Tsit5(),save_everystep=false)

concentrations = sol[1:7,:]
Vₛ⁻ = cellgeometry.Vₛ⁻
Vₛ⁺ = cellgeometry.Vₛ⁺
Vₑ⁻ = cellgeometry.Vₑ⁻
Vₑˢ = cellgeometry.Vₑˢ
Vₑ⁺ = cellgeometry.Vₑ⁺

@unpack εₛ⁻,εₛ⁺,δ⁻,δ⁺,εₑˢ,R⁺,R⁻ = p
X⁺ = ((R⁺+δ⁺)^3-R⁺^3)/(R⁺^3)
X⁻ = ((R⁻+δ⁻)^3-R⁻^3)/(R⁻^3)
εₑ⁻ = 1-(1+X⁻)εₛ⁻
εₑ⁺ = 1-(1+X⁺)εₛ⁺

mass_vector = [Vₛ⁻*εₛ⁻,Vₛ⁻*εₛ⁻,Vₑ⁻*εₑ⁻,Vₑˢ*εₑˢ,Vₑ⁺*εₑ⁺,Vₛ⁺*εₛ⁺,Vₛ⁺*εₛ⁺]
n_li = concentrations'*mass_vector
n_li_change = maximum(abs.((n_li.-(n_li[1]))./n_li[1]))

V = CellFitElectrolyte.calc_voltage(sol,p,sol.t,cache,cellgeometry,cathodeocv,anodeocv,p.input_value.*ones(length(sol.t)),297 .*ones(length(sol.t)))


x_s_n = sol[1,:]./anodeocv.c_s_max
x_b_n = sol[2,:]./anodeocv.c_s_max

x_s_p = sol[6,:]./cathodeocv.c_s_max
x_b_p = sol[7,:]./cathodeocv.c_s_max

myOCV = cathodeocv.(x_s_p,297.15.*ones(length(sol.t))).-anodeocv.(x_s_n,297.15.*ones(length(sol.t)))


@test n_li_change<1e-10




