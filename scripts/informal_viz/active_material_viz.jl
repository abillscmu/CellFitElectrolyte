function get_n_li_from_u_and_p(u, p)
    cₛˢ⁻ = u[1]
    cₛᵇ⁻ = u[2]
    cₑ⁻ = u[3]
    cₑˢ = u[4]
    cₑ⁺ = u[5]
    cₛᵇ⁺ = u[6]
    cₛˢ⁺ = u[7]

    @unpack Vₛ⁻,Vₛ⁺,Vₑ⁻,Vₑˢ,Vₑ⁺  = cellgeometry
    @unpack εₛ⁻,εₛ⁺,δ⁻,δ⁺,εₑˢ,Temp = p
    @unpack R⁺,R⁻= p
    X⁺ = ((R⁺+δ⁺)^3-R⁺^3)/(R⁺^3)
    X⁻ = ((R⁻+δ⁻)^3-R⁻^3)/(R⁻^3)
    εₑ⁻ = 1-(1+X⁻)εₛ⁻ 
    εₑ⁺ = 1-(1+X⁺)εₛ⁺

    Veffₛ⁻ = εₛ⁻*Vₛ⁻
    Veffₛ⁺ = εₛ⁺*Vₛ⁺
    Veffₑ⁻ = εₑ⁻*Vₑ⁻
    Veffₑˢ = εₑˢ*Vₑˢ
    Veffₑ⁺ = εₑ⁺*Vₑ⁺
    
    n_li_m = Veffₛ⁻*cₛˢ⁻ + Veffₛ⁻*cₛᵇ⁻ + Veffₑ⁻*cₑ⁻ + Veffₑˢ*cₑˢ + Veffₑ⁺*cₑ⁺ + Veffₛ⁺*cₛˢ⁺ + Veffₛ⁺*cₛᵇ⁺

    return n_li_m
end


