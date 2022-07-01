#u:
#Concentrations
# cₛˢ⁻   [1]
# cₛᵇ⁻   [2]
# cₑ⁻    [3]
# cₑˢ    [4]
# cₑ⁺    [5]
# cₛᵇ⁺    [6]
# cₛˢ⁺    [7]
#Volume Fractions
# εₛ⁻   [8]
# εₑ⁺   [9]
# εₑ⁻   [10]
# εₑˢ   [11]
# εₛ⁺    [12]
#Thermals
# T     [13]

function equations_electrolyte(du,u,p,t,cache,cellgeometry,cathodeocv,anodeocv)
    Iapp = u[14]
    cₛˢ⁻,cₛᵇ⁻,cₑ⁻,cₑˢ,cₑ⁺,cₛᵇ⁺,cₛˢ⁺ = @view u[1:7]
    @unpack Tamb = p


    #Transport Parameters
    @unpack θₛ⁻,θₑ,θₛ⁺ = p
    @unpack β⁻,βˢ,β⁺ = p
    #Kinetic Parameters
    @unpack k₀⁺,k₀⁻ = p
    #Thermal Parameters
    @unpack c,h = p
    #Electrolyte Parameters
    @unpack κ,t⁺ = p

    #Geometry
    @unpack R⁺,R⁻ = p
    @unpack Vₛ⁻,Vₛ⁺,T⁺,T⁻  = cellgeometry

    εₛ⁻,εₑ⁻,εₑˢ,εₑ⁺,εₛ⁺ = @view u[8:12]
    T = u[13]

    #Fill transport and apply bruggeman corrections (and temp later)
    fill_transport!(cache.A,θₛ⁻,θₑ,θₛ⁺,εₑ⁻,εₑ⁺,β⁻,β⁺)

    #Calculate Transport Matrix
    mul!(cache.du_transport,cache.A,@view u[1:7])

    #Apply Current
    mul!(cache.control,cache.controller,Iapp)
    du[1:7].=cache.du_transport.+cache.control

    #Apply corrections for temperature to diffusion



    #Apply Mass Matrix
    volume_correction!(cache.mm_cache,cellgeometry,εₛ⁻, εₑ⁻,εₑˢ,εₑ⁺,εₛ⁺)
    du[1:7] .= (@view du[1:7]).*cache.mm_cache

    #Current Density
    a⁻ = 3*εₛ⁻/R⁻
    a⁺ = 3*εₛ⁺/R⁺
    A⁻ = 2Vₛ⁻*a⁻
    A⁺ = 2Vₛ⁺*a⁺
    J⁻ = Iapp/A⁻
    J⁺ = Iapp/A⁺

    #Calculate Voltages
    #U⁺ = cathodeocv((cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),T)
    #U⁻ = anodeocv((cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),T)
    J₀⁻ = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺ = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊ = butler_volmer(J₀⁺,J⁺,T)
    η₋ = butler_volmer(J₀⁻,J⁻,T)
    
    ηc₋ = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,T,T⁻)
    ηc₊ = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,T,T⁺)
    
    ηₒ₋ = electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,T⁻)
    ηₒ₊ = electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,T⁺)

    η = η₊+η₋+ηc₋+ηc₊+ηₒ₋+ηₒ₊
    du[13] = (η*Iapp-h*(T-Tamb))/c





    du[8]=0
    du[9]=0
    du[10]=0
    du[11]=0
    du[12]=0
    @unpack input_type,input_value = p
        #Calculate Current
        if input_type==0
            du[14] = Iapp-0
            # println(Iapp)
        elseif input_type==1
            U⁺ = calcocv(cathodeocv,(cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),T)
            U⁻ = calcocv(anodeocv,(cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),T)
            Voltage = U⁺-U⁻-η
            du[14] = input_value-(Iapp*Voltage)
        elseif input_type==2
            U⁺ = calcocv(cathodeocv,(cₛˢ⁺-cathodeocv.c_s_min)/(cathodeocv.c_s_max-cathodeocv.c_s_min),T)
            U⁻ = calcocv(anodeocv,(cₛˢ⁻-anodeocv.c_s_min)/(anodeocv.c_s_max-anodeocv.c_s_min),T)
            Voltage = U⁺-U⁻-η
            du[14] = input_value-Voltage 
        elseif input_type==3
            du[14] = Iapp-input_value;
        elseif input_type==4
            du[14] = Iapp-0
        else
            @warn "condition not recognized"
            du[14] = Iapp-0
        end
    #Thermal Equations
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




function exchange_current_density(cₛ,cₑ,c_s_max,k₀)
    J₀ = k₀.*sqrt.(cₑ.*(c_s_max.-cₛ).*cₛ)
    return J₀
end

function butler_volmer(J₀,J,T)
    return @fastmath (R.*T./(0.5.*F)).*asinh.(J./(2J₀))
end

function concentration_overpotential(cₑ⁺,cₑ⁻,t⁺,T,Δx)
    @fastmath ΔΦ = 2 .*(1 .-t⁺).*R.*T.*log.(cₑ⁺./cₑ⁻)./(F)
    return ΔΦ
end

function electrolyte_ohmic(ε,β,κ,Iapp,Δx)
    return Δx.*Iapp./(ε.^β.*κ)
end

function calc_voltage(sol,p,t::Array,cache,cellgeometry,cathodeocv,anodeocv)
    Iapp = @view sol[14,:]

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
    @unpack R⁺,R⁻= p
    @unpack Vₛ⁻,Vₛ⁺,T⁺,T⁻  = cellgeometry

    εₛ⁻ = @view sol[8,:]
    εₑ⁻ = @view sol[9,:]
    εₑˢ = @view sol[10,:]
    εₑ⁺ = @view sol[11,:]
    εₛ⁺ = @view sol[12,:]
    T = @view sol[13,:]

    #Current Density
    a⁻ = 3 .*εₛ⁺./R⁻
    a⁺ = 3 .*εₛ⁺/R⁺
    A⁻ = 2 .*Vₛ⁻.*a⁻
    A⁺ = 2 .*Vₛ⁺.*a⁺
    J⁻ = Iapp./A⁻
    J⁺ = Iapp./A⁺

    #Calculate Voltages
    U⁺ = cathodeocv.((cₛˢ⁺.-cathodeocv.c_s_min)./(cathodeocv.c_s_max-cathodeocv.c_s_min),T)
    U⁻ = anodeocv.((cₛˢ⁻.-anodeocv.c_s_min)./(anodeocv.c_s_max-anodeocv.c_s_min),T)
    
    J₀⁻ = exchange_current_density(cₛˢ⁻,cₑ⁻,anodeocv.c_s_max,p.k₀⁻)
    J₀⁺ = exchange_current_density(cₛˢ⁺,cₑ⁺,cathodeocv.c_s_max,p.k₀⁺)
    
    η₊ = butler_volmer(J₀⁺,J⁺,T)
    η₋ = butler_volmer(J₀⁻,J⁻,T)
    
    ηc₋ = concentration_overpotential(cₑ⁻,cₑˢ,t⁺,T,T⁻)
    ηc₊ = concentration_overpotential(cₑˢ,cₑ⁺,t⁺,T,T⁺)
    
    ηₒ₋ = electrolyte_ohmic(εₑ⁻,β⁻,κ,Iapp,T⁻)
    ηₒ₊ = electrolyte_ohmic(εₑ⁺,β⁺,κ,Iapp,T⁺)
    #Thermal Equations
    V = U⁺.-U⁻.-η₊.-η₋.-ηc₋.-ηc₊.-ηₒ₋.-ηₒ₊
    return V
end


function calc_voltage(u::Array{T,1},p::ComponentVector{T},t::T,cache::cache{T},cellgeometry::ComponentVector{T},cathodeocv::RKPolynomial{Vector{T},T},anodeocv::RKPolynomial{Vector{T},T}) where {T}
    Iapp::T = u[14]

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

    εₛ⁻::T = u[8]
    εₑ⁻::T = u[9]
    εₑˢ::T = u[10]
    εₑ⁺::T = u[11]
    εₛ⁺::T = u[12]
    Temp::T = u[13]

    #Current Density
    a⁻::T = 3 *εₛ⁺/R⁻
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
    #Thermal Equations
    V::T = U⁺-U⁻-η₊-η₋-ηc₋-ηc₊-ηₒ₋-ηₒ₊
    return V
end

