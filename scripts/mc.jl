using KernelDensity

#expect to have a chain (with ω)
ω = 0.01

p = ComponentVector(θₛ⁻ = 3.238105814128935e-8, θₑ = 5.6464068552786306e-7, θₛ⁺ = 6.547741580032837e-5, R⁺ = 4.2902932816468984e-6, R⁻ = 1.7447548850488327e-6, β⁻ = 1.5, β⁺ = 1.5, βˢ = 1.5, εₛ⁻ = 0.75, εₛ⁺ = 0.6, εᵧ⁺ = 0, εᵧ⁻ = 0, c = 50.0, h = 0.1, Tamb = 298.15, Temp = 298.15, k₀⁺ = 0.002885522176210856, k₀⁻ = 1.7219544782420964, x⁻₀ = 0.6, εₑˢ = 0.8, cₑ₀ = 4175.451281358547, κ = 0.2025997972168558, t⁺ = 0.38, input_type = 3.0, input_value = 4.2, ω = ω, Eₑ = 50.0, Eₛ⁺ = 50.0, Eₛ⁻ = 50.0, cccv_switch_2 = false, cccv_switch = false, vfull = 4.2, ifull = -0.05)

predicted = cycle_array_evaluator(p)


figure(1);
pygui(true)
clf()
plot(df.times, df.EcellV)
plot(df.times, predicted)