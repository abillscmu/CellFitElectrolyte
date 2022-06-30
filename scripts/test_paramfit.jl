using CellFitElectrolyte
using CellFitElectrolyte.ComponentArrays
using CellFitElectrolyte.OrdinaryDiffEq
using CellFitElectrolyte.OCV
using Test

type = Float64

cache = CellFitElectrolyte.initialize_cache(type)
p = CellFitElectrolyte.p_transport()
CellFitElectrolyte.fill_transport!(cache.A,p.θₛ⁻,p.θₑ,p.θₛ⁺)

cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()

initialcond = Dict(
    "Starting Voltage[V]"=>3.9,
    "Ambient Temperature[K]"=>298.0,
    "Capacity[mAh]"=>3000.0)

cellgeometry = CellFitElectrolyte.cell_geometry()



f(t) = 8.4


u  = Array{type,1}(undef,13)

u = CellFitElectrolyte.initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
du = similar(u)
end_dt = 100.0

prob = ODEProblem((du,u,p,t)->CellFitElectrolyte.equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv,f),u,(0.0,end_dt),p,dtmax=0.001)

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

V = CellFitElectrolyte.calc_voltage(sol,p,sol.t,cache,cellgeometry,cathodeocv,anodeocv,f)


x_s_n = sol[1,:]./anodeocv.c_s_max
x_b_n = sol[2,:]./anodeocv.c_s_max

x_s_p = sol[6,:]./cathodeocv.c_s_max
x_b_p = sol[7,:]./cathodeocv.c_s_max

myOCV = cathodeocv.(x_s_p,sol[13,:]).-anodeocv.(x_s_n,sol[13,:])

using MATLABPlots
figure(1)
clf()
subplot(4,1,1)
plot(sol.t,sol[1,:]./anodeocv.c_s_max)
hold_on()
plot(sol.t,sol[2,:]./anodeocv.c_s_max)
legend(["Surface Concentration","Bulk Concentration"])
ylabel("x^{-}")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))
grid_on()
subplot(4,1,2)
plot(sol.t,sol[3,:])
hold_on()
plot(sol.t,sol[4,:])
plot(sol.t,sol[5,:])
legend(["Negative","Separator","Positive"])
ylabel("c_{e}")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))
grid_on()
subplot(4,1,3)
plot(sol.t,sol[6,:]./cathodeocv.c_s_max)
hold_on()
plot(sol.t,sol[7,:]./cathodeocv.c_s_max)
legend(["Surface","Bulk"])
ylabel("x^{+}")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))
grid_on()
subplot(4,1,4)
plot(sol.t,V)
hold_on()
plot(sol.t,myOCV)
legend(["Voltage","OCV"])
ylabel("Voltage")
xlabel("Time[s]")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))
grid_on()

figure(2)
clf()
plot(sol.t,myOCV.-V)
xlabel("Time[s]")
ylabel("\\eta")
setgca(Dict("FontName"=>"Open Sans","FontSize"=>16))
grid_on()


@test n_li_change<1e-10




