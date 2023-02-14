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

function equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv)
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

    #Fill transport and apply bruggeman corrections (and temp later)
    fill_transport!(cache.A,θₛ⁻,θₑ,θₛ⁺,εₑ⁻,εₑ⁺,β⁻,β⁺)
    

    #Apply corrections for temperature to diffusion
    arrhenius!(view(cache.A, 1:2, 1:2),Eₛ⁻,Temp)
    arrhenius!(view(cache.A, 3:5, 3:5),Eₑ,Temp)
    arrhenius!(view(cache.A, 6:7, 6:7),Eₛ⁺,Temp)

    #Calculate Transport Matrix
    mul!(cache.du_transport,cache.A,@view u[1:7])

    #Apply Current
    mul!(cache.control,cache.controller,Iapp)
    du[1:7].=cache.du_transport.+cache.control

    #Apply Mass Matrix
    volume_correction!(cache.mm_cache,cellgeometry,εₛ⁻, εₑ⁻,εₑˢ,εₑ⁺,εₛ⁺)
    du[1:7] .= (@view du[1:7]).*cache.mm_cache
    return nothing
end

function fill_transport!(A,θₛ⁻,θₑ,θₛ⁺)
    #Positive Electrode Solid Concentration
    A[1,1] = -θₛ⁻
    A[1,2] = θₛ⁻
    A[2,1] = θₛ⁻
    A[2,2] = -θₛ⁻

    #Electrolyte Concentration
    A[3,3] = -θₑ
    A[3,4] = θₑ
    A[4,3] = θₑ
    A[4,4] = -2θₑ
    A[4,5] = θₑ
    A[5,4] = θₑ
    A[5,5] = -θₑ

    #Negative Electrode Solid Concentration
    A[6,6] = -θₛ⁺
    A[6,7] = θₛ⁺
    A[7,6] = θₛ⁺
    A[7,7] = -θₛ⁺
end

function fill_transport!(A,θₛ⁻,θₑ,θₛ⁺,εₑ⁻,εₑ⁺,β⁻,β⁺)
    #Positive Electrode Solid Concentration
    A[1,1] = -θₛ⁻
    A[1,2] = θₛ⁻
    A[2,1] = θₛ⁻
    A[2,2] = -θₛ⁻

    @fastmath neg = εₑ⁻^β⁻
    @fastmath pos = εₑ⁺^β⁺   

    #Electrolyte Concentration
    A[3,3] = -θₑ*neg
    A[3,4] = θₑ*neg
    A[4,3] = θₑ*neg
    A[4,4] = -θₑ*(neg+pos)
    A[4,5] = θₑ*pos
    A[5,4] = θₑ*pos
    A[5,5] = -θₑ*pos

    #Negative Electrode Solid Concentration
    A[6,6] = -θₛ⁺
    A[6,7] = θₛ⁺
    A[7,6] = θₛ⁺
    A[7,7] = -θₛ⁺
end


#ASSUMPTION: f=1
function volume_correction!(mm_cache,cellgeometry,εₛ⁻,εₑ⁻,εₑˢ,εₑ⁺,εₛ⁺)
    mm_cache[1:2] .= 1/(cellgeometry.Vₛ⁻*εₛ⁻)
    mm_cache[3] = 1/(cellgeometry.Vₑ⁻*εₑ⁻)
    mm_cache[4] = 1/(cellgeometry.Vₑˢ*εₑˢ)
    mm_cache[5] = 1/(cellgeometry.Vₑ⁺*εₑ⁺)
    mm_cache[6:7] .= 1/(cellgeometry.Vₛ⁺*εₛ⁺)
end

function arrhenius!(A,E,Temp)
    correction = exp(E/293-E/Temp)
    #note: can only get away with this when correction is a scalar
    mul!(A,A,correction)
end

function sei_ohmic(ω,Iapp)
    η_sei = ω*Iapp
    return η_sei
end

function sei_growth(c_s_bulk, k_sei, α_sei, Temp, η_sei, D_sei, δ)
    kinetic_growth = CellFitElectrolyte.F*k_sei*exp(-α_sei*CellFitElectrolyte.F*η_sei/(CellFitElectrolyte.R*Temp))
    diffusive_growth = δ / (CellFitElectrolyte.F * D_sei)
    i_sei = c_s_bulk / ((1 / kinetic_growth) + diffusive_growth)
    return i_sei
end

function sei_delta_growth(I_sei, δ, R, εₑ, Vₘ)
    δ̇ = ((R + δ)/(3*(1-εₑ)))*(I_sei/2*CellFitElectrolyte.F)*Vₘ
    return δ̇
end

function plating(j0_pl, a, α_pl, Temp, η_pl)
    i_pl = a*j0_pl*exp(-α_pl*CellFitElectrolyte.F*η_pl/(CellFitElectrolyte.R*Temp))
    return i_pl
end

function active_material_dissolution(i0_diss, a, α_diss, Temp, η_diss)
    i_diss = a*i0_diss*exp(-α_diss*CellFitElectrolyte.F*η_diss/(CellFitElectrolyte.R*Temp))
    return i_diss
end

function active_material_dissolution_volfrac(i_diss, q, n)
    dvolfrac_dt = i_diss/(CellFitElectrolyte.F*q*n)
    return dvolfrac_dt
end


function exchange_current_density(cₛ,cₑ,c_s_max,k₀)
    J₀ = k₀.*sqrt.(cₑ.*(c_s_max.-cₛ).*cₛ)
    return J₀
end

function butler_volmer(J₀,J,T)
    return @fastmath (R.*T./(0.5.*F)).*asinh.(J./(2J₀))
end

function concentration_overpotential(cₑ⁺::Array,cₑ⁻::Array,t⁺,Temp,Δx)
    @fastmath ΔΦ = 2 .*(1 .-t⁺).*R.*Temp.*log.(cₑ⁺./cₑ⁻)./(F)
    return ΔΦ
end

function concentration_overpotential(cₑ⁺,cₑ⁻,t⁺,Temp,Δx)
    @fastmath ΔΦ = 2 *(1 -t⁺)*R*Temp*log(cₑ⁺/cₑ⁻)/(F)
    return ΔΦ
end

function electrolyte_ohmic(ε,β,κ,Iapp::Array,Δx)
    @fastmath η = Δx.*Iapp./(ε.^β.*κ)
    return η
end

function electrolyte_ohmic(ε,β,κ,Iapp,Δx)
    @fastmath η = Δx*Iapp/(ε^β*κ)
    return η
end

function calc_voltage(sol,p,t::Array,cache,cellgeometry,cathodeocv,anodeocv,Iapp::Array,Temp::Array)

    cₛˢ⁻ = @view sol[1,:]
    cₛᵇ⁻ = @view sol[2,:]
    cₑ⁻ = @view sol[3,:]
    cₑˢ = @view sol[4,:]
    cₑ⁺ = @view sol[5,:]
    cₛᵇ⁺ = @view sol[6,:]
    cₛˢ⁺ = @view sol[7,:]


    #Transport Parameters
    @unpack θₛ⁻,θₑ,θₛ⁺ = p
    #Kinetic Parameters
    @unpack k₀⁺,k₀⁻ = p
    #Thermal Parameters
    @unpack c,h = p
    #Electrolyte Parameters
    @unpack κ,t⁺,β⁻,β⁺ = p

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


    J⁻ = Iapp./A⁻
    J⁺ = Iapp./A⁺

    #Calculate Voltages
    U⁺ = cathodeocv.((cₛˢ⁺.-cathodeocv.c_s_min)./(cathodeocv.c_s_max-cathodeocv.c_s_min),Temp)
    U⁻ = anodeocv.((cₛˢ⁻.-anodeocv.c_s_min)./(anodeocv.c_s_max-anodeocv.c_s_min),Temp)
    
    J₀⁻ = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺ = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊ = butler_volmer(J₀⁺,J⁺,Temp)
    η₋ = butler_volmer(J₀⁻,J⁻,Temp)
    
    ηc₋ = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,Temp,T⁻)
    ηc₊ = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,Temp,T⁺)
    
    ηₒ₋ = electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,T⁻)
    ηₒ₊ = electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,T⁺)

    ηₛ = sei_ohmic(ω,Iapp)
    #Thermal Equations
    V = U⁺.-U⁻.-η₊.-η₋.-ηc₋.-ηc₊.-ηₒ₋.-ηₒ₊.-ηₛ
    return V
end


function calc_voltage(u::Array{T,1},p::ComponentVector{T},t::T,cache::cache{T},cellgeometry::ComponentVector{T},cathodeocv::OCV.OpenCircuitVoltage,anodeocv::OCV.OpenCircuitVoltage,Iapp) where {T}

    cₛˢ⁻::T = u[1]
    cₛᵇ⁻::T = u[2]
    cₑ⁻::T = u[3]
    cₑˢ::T = u[4]
    cₑ⁺::T = u[5]
    cₛᵇ⁺::T = u[6]
    cₛˢ⁺::T = u[7]


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
    a⁻::T = 3 *εₛ⁻/R⁻
    a⁺::T = 3 *εₛ⁺/R⁺
    A⁻::T = 2 *Vₛ⁻*a⁻
    A⁺::T = 2 *Vₛ⁺*a⁺
    J⁻::T = Iapp/A⁻
    J⁺::T = Iapp/A⁺

    #Calculate Voltages
    U⁺::T = calcocv(cathodeocv,(cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),Temp)
    U⁻::T = calcocv(anodeocv,(cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),Temp)
    
    J₀⁻::T = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺::T = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊::T = butler_volmer(J₀⁺,J⁺,Temp)
    η₋::T = butler_volmer(J₀⁻,J⁻,Temp)
    
    ηc₋::T = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,Temp,T⁻)
    ηc₊::T = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,Temp,T⁺)
    
    ηₒ₋::T = electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,T⁻)
    ηₒ₊::T = electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,T⁺)

    #Calculate SEI Overpotential
    ηₛ = sei_ohmic(ω,Iapp)
    #Thermal Equations
    V::T = U⁺-U⁻-η₊-η₋-ηc₋-ηc₊-ηₒ₋-ηₒ₊-ηₛ
    return V
end


function equations_electrolyte_life(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv)
    cₛˢ⁻,cₛᵇ⁻,cₑ⁻,cₑˢ,cₑ⁺,cₛᵇ⁺,cₛˢ⁺ = @view u[1:7]
    @unpack Temp = p
    @unpack input_type,input_value = p
    Iapp = u[8]


    #Transport Parameters
    @unpack θₛ⁻,θₑ,θₛ⁺ = p
    @unpack β⁻,βˢ,β⁺ = p
    @unpack Eₛ⁺, Eₛ⁻, Eₑ = p
    #Kinetic Parameters
    @unpack k₀⁺,k₀⁻,ω = p
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

    #Fill transport and apply bruggeman corrections (and temp later)
    fill_transport!(cache.A,θₛ⁻,θₑ,θₛ⁺,εₑ⁻,εₑ⁺,β⁻,β⁺)

    #Calculate Transport Matrix
    mul!(cache.du_transport,cache.A,@view u[1:7])

    #Apply Current
    mul!(cache.control,cache.controller,Iapp)
    du[1:7].=cache.du_transport.+cache.control

    #Apply corrections for temperature to diffusion
    arrhenius!(view(cache.A, 1:2, 1:2),Eₛ⁻,Temp)
    arrhenius!(view(cache.A, 3:5, 3:5),Eₑ,Temp)
    arrhenius!(view(cache.A, 6:7, 6:7),Eₛ⁺,Temp)


    #Apply Mass Matrix
    volume_correction!(cache.mm_cache,cellgeometry,εₛ⁻, εₑ⁻,εₑˢ,εₑ⁺,εₛ⁺)
    du[1:7] .= (@view du[1:7]).*cache.mm_cache

 

    #Calculate Voltages
    U⁺ = cathodeocv((cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),Temp)
    U⁻ = anodeocv((cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),Temp)
    
    
    J₀⁻ = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺ = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊ = butler_volmer(J₀⁺,J⁺,Temp)
    η₋ = butler_volmer(J₀⁻,J⁻,Temp)
    
    ηc₋ = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,Temp,T⁻)
    ηc₊ = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,Temp,T⁺)
    
    ηₒ₋ = electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,T⁻)
    ηₒ₊ = electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,T⁺)

    ηₛ = sei_ohmic(ω, Iapp)

    η = η₊+η₋+ηc₋+ηc₊+ηₒ₋+ηₒ₊-ηₛ


        #Calculate Current
    if input_type==0
        du[8] = Iapp-0
    elseif input_type==1
        Voltage = U⁺-U⁻-η
        du[8] = input_value-(Iapp*Voltage)
    elseif input_type==2
        Voltage = U⁺-U⁻-η
        du[8] = input_value-Voltage 
    elseif input_type==3
        du[8] = Iapp-input_value;
    elseif input_type==4
        du[8] = Iapp-0
    elseif input_type==5
        Voltage = U⁺-U⁻-η
        if p.cccv_switch == true
            if p.cccv_switch_2 == true
                du[8] = Iapp - 0
            else
                du[8] = Voltage - p.vfull
            end
        else
            du[8] = Iapp - input_value
        end
    else
        @warn "condition not recognized"
        du[8] = Iapp-0
    end
    #Thermal Equations
    return nothing
end

