#u:
#Concentrations
# cₛˢ⁻   [1]
# cₛᵇ⁻   [2]
# cₑ⁻    [3]
# cₑˢ    [4]
# cₑ⁺    [5]
# cₛᵇ⁺    [6]
# cₛˢ⁺    [7]
# IF S.M.:
# Iapp [8]
#Volume Fractions
# εₛ⁻   [9]
# εₑ⁺   [10]
# εₑ⁻   [11]
# εₑˢ   [12]
# εₛ⁺    [13]
#Thermals
# T     [14]

function equations_electrolyte_allocating(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv)
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
    A_E = SMatrix{3,3}(-θₑ*neg, θₑ*neg, 0, θₑ*neg, -θₑ*(neg+pos), θₑ*pos, 0, θₑ*pos, -θₑ*pos) * exp(Eₑ/293-Eₑ/Temp)

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

function calc_voltage( u, p, t, cache, cellgeometry, cathodeocv::RKPolynomial, anodeocv, Iapp)

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

    #Calculate Voltages
    U⁺ = calcocv(cathodeocv,(cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),Temp)
    U⁻ = calcocv(anodeocv,(cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),Temp)
    
    J₀⁻ = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺ = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊ = butler_volmer(J₀⁺,J⁺,Temp)
    η₋ = butler_volmer(J₀⁻,J⁻,Temp)
    
    ηc₋ = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,Temp,T⁻)
    ηc₊ = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,Temp,T⁺)
    
    ηₒ₋ = electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,T⁻)
    ηₒ₊ = electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,T⁺)

    #Calculate SEI Overpotential
    ηₛ = sei_ohmic(ω,Iapp)
    #Thermal Equations
    V = U⁺-U⁻-η₊-η₋-ηc₋-ηc₊-ηₒ₋-ηₒ₊-ηₛ
    return V
end


function get_eps_gamma_from_δ(R, δ, εₛ)
    X = ((R+δ)^3-R^3)/(R^3)
    εₑ = 1-(1+X)εₛ
    εᵧ = 1 - εₑ - εₛ
    return εᵧ
end