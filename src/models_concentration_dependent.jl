
function equations_electrolyte_concentrationdependent(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv)
    cₛˢ⁻,cₛᵇ⁻,cₑ⁻,cₑˢ,cₑ⁺,cₛᵇ⁺,cₛˢ⁺ = @view u[1:7]
    @unpack Temp = p
    @unpack input_type,input_value = p
    Iapp = input_value


    #Transport Parameters
    @unpack θₛ⁻,θₑ,θₛ⁺ = p
    @unpack β⁻,βˢ,β⁺ = p
    @unpack Eₛ⁺, Eₛ⁻, Eₑ = p
    #Kinetic Parameters
    @unpack k₀⁺,k₀⁻ = p
    #Thermal Parameters
    @unpack c,h = p
    #Electrolyte Parameters
    @unpack κ,t⁺ = p

    #Geometry
    @unpack R⁺,R⁻ = p
    @unpack Vₛ⁻,Vₛ⁺,T⁺,T⁻  = cellgeometry
    @unpack εₛ⁻,εₛ⁺,εₑˢ = p
    @unpack εᵧ⁺,εᵧ⁻ = p
    εₑ⁻ = 1-εₛ⁻-εᵧ⁻
    εₑ⁺ = 1-εₛ⁺-εᵧ⁺


    #Geometry
    a⁻ = 3*εₛ⁻/R⁻
    a⁺ = 3*εₛ⁺/R⁺
    A⁻ = 2Vₛ⁻*a⁻
    A⁺ = 2Vₛ⁺*a⁺
    
    J⁻ = Iapp/A⁻
    J⁺ = Iapp/A⁺

    @fastmath pos = εₑ⁺^β⁺
    @fastmath neg = εₑ⁻^β⁻

    #Fill transport and apply bruggeman corrections (and temp later)
    A_NS = SMatrix{2,2}(-θₛ⁻, θₛ⁻, θₛ⁻, -θₛ⁻) * exp(Eₛ⁻/293-Eₛ⁻/Temp)
    A_PS = SMatrix{2,2}(-θₛ⁺, θₛ⁺, θₛ⁺, -θₛ⁺) * exp(Eₛ⁺/293-Eₛ⁺/Temp)

    θₑ⁻ = 2*electrolyte_diffusivity(cₑ⁻, Temp)/T⁻
    θₑ⁺ = 2*electrolyte_diffusivity(cₑ⁺, Temp)/T⁺
    θₑˢ = 2*electrolyte_diffusivity(cₑˢ, Temp)/T⁺

    A_E = SMatrix{3,3}(-θₑ⁻, θₑ⁻*neg, 0, θₑ⁻*neg, -θₑˢ*(neg+pos), θₑ⁺*pos, 0, θₑ⁺*pos, -θₑ⁺*pos)

    du_NS = A_NS * (@view u[1:2])
    du_E = A_E * (@view u[3:5])
    du_PS = A_PS * (@view u[6:7])
    du_transport = vcat(du_NS, du_E, du_PS)

    #Apply current
    cache_control = cache.controller*Iapp
    du_pre_mm = du_transport + cache_control

    mm_sn = 1/(cellgeometry.Vₛ⁻*εₛ⁻)
    mm_sp = 1/(cellgeometry.Vₛ⁺*εₛ⁺)

    mm = SVector(mm_sn, mm_sn, 1/(cellgeometry.Vₑ⁻*εₑ⁻), 1/(cellgeometry.Vₑˢ*εₑˢ), 1/(cellgeometry.Vₑ⁺*εₑ⁺), mm_sp, mm_sp)

    du[1:7] .= du_pre_mm .* mm
    return nothing
end


#Johannes Landesfeind and Hubert A. Gasteiger 2019 J. Electrochem. Soc. 166 A3079
function electrolyte_conductivity(c, T)
    #c is in mol/m^3
    c = c/1000
    p₁ = 7.98e-1
    p₂ = 2.28e2
    p₃ = -1.22e0
    p₄ = 5.09e-1
    p₅ = -4e-3
    p₆ = 3.79e-3
    κ = p₁*(1 + (T - p₂))*c*(1 + ((p₃*sqrt(c))+p₄*(1+p₅*exp(1000/T))*c))/(1 + ((c^4)*(p₆*exp(1000/T))))
    #κ is in mS/cm:
    κ = κ/10
    return κ
end

function electrolyte_diffusivity(c, T)
    #c is in mol/m^3
    c = c/1000
    p₁ = 1.47e3
    p₂ = 1.33e0
    p₃ = -1.69e3
    p₄ = -5.63e2
    D = (10^-6)*(p₁*exp(p₂*c)*exp(p₃/T)*exp(p₄*c/T))
    #D is in cm^2/s
    D = D*.0001
end

function electrolyte_TDF(c, T)
    c = c/1000

    p₁ = -5.58e0
    p₂ = 7.17e0
    p₃ = 3.8e-2
    p₄ = 1.91e0
    p₅ = -6.65e-2
    p₆ = -5.08e-5
    p₇ = 1.1e-1
    p₈ = -6.1e-3
    p₉ = 1.51e-4

    TDF = p₁ + (p₂*c) + (p₃*T) + (p₄*(c^2)) + (p₅*c*T) + (p₆*(T^2)) + (p₇*(c^3)) + (p₈*(c^2)*T) + (p₉*c*(T^2))
    return 2
end

function concentration_overpotential(cₑ⁺,cₑ⁻,t⁺,Temp,Δx, TDF)
    @fastmath ΔΦ = TDF .*(1 .-t⁺).*R.*Temp.*log.(cₑ⁺./cₑ⁻)./(F)
    return ΔΦ
end


function calc_voltage_concentrationdependent( u, p, t, cache, cellgeometry, cathodeocv::RKPolynomial, anodeocv, Iapp)

    cₛˢ⁻ = u[1]
    cₛᵇ⁻ = u[2]
    cₑ⁻ = u[3]
    cₑˢ = u[4]
    cₑ⁺ = u[5]
    cₛᵇ⁺ = u[6]
    cₛˢ⁺ = u[7]


    #Transport Parameters
    @unpack θₛ⁻,θₑ,θₛ⁺ = p
    #Kinetic Parameters
    @unpack k₀⁺,k₀⁻ = p
    #Thermal Parameters
    @unpack c,h = p
    #Electrolyte Parameters
    @unpack κ,t⁺,β⁻,β⁺ = p

    

    #Geometry
    @unpack R⁺,R⁻= p
    @unpack Vₛ⁻,Vₛ⁺,T⁺,T⁻  = cellgeometry
    @unpack εₛ⁻,εₛ⁺,εₑˢ,Temp = p
    @unpack ω = p
    @unpack εᵧ⁺,εᵧ⁻ = p

    εₑ⁻ = 1-εₛ⁻-εᵧ⁻
    εₑ⁺ = 1-εₛ⁺-εᵧ⁺

    #Current Density
    a⁻ = 3 *εₛ⁻/R⁻
    a⁺ = 3 *εₛ⁺/R⁺
    A⁻ = 2 *Vₛ⁻*a⁻
    A⁺ = 2 *Vₛ⁺*a⁺
    J⁻ = Iapp/A⁻
    J⁺ = Iapp/A⁺

    κ⁺ = electrolyte_conductivity(cₑ⁺, Temp)
    κ⁻ = electrolyte_conductivity(cₑ⁻, Temp)

    TDF⁺ = electrolyte_TDF(cₑ⁺, Temp)
    TDF⁻ = electrolyte_TDF(cₑ⁻, Temp)

    #Calculate Voltages
    U⁺ = calcocv(cathodeocv,(cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),Temp)
    U⁻ = calcocv(anodeocv,(cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),Temp)
    
    J₀⁻ = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺ = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊ = butler_volmer(J₀⁺,J⁺,Temp)
    η₋ = butler_volmer(J₀⁻,J⁻,Temp)
    
    ηc₋ = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,Temp,T⁻, TDF⁻)
    ηc₊ = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,Temp,T⁺, TDF⁺)
    
    ηₒ₋ = electrolyte_ohmic(εₑ⁻,β⁻,κ⁻,Iapp,T⁻)
    ηₒ₊ = electrolyte_ohmic(εₑ⁺,β⁺,κ⁺,Iapp,T⁺)

    #Calculate SEI Overpotential
    ηₛ = sei_ohmic(ω,Iapp)
    #Thermal Equations
    V = U⁺-U⁻-η₊-η₋-ηc₋-ηc₊-ηₒ₋-ηₒ₊-ηₛ
    return V
end