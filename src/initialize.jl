function p_transport()
    p = ComponentArray(
        #Transport
        θₛ⁻ = 6.5e-9,
        θₑ = 4e-6,
        θₛ⁺ = 5e-8,
        R⁺ = 1e-6,
        R⁻ = 1e-6,
        β⁻ = 1.5,
        β⁺ = 1.5,
        βˢ = 1.5,
        εₛ⁻ = 0.6,
        εₛ⁺ = 0.6,
        δ⁻ = 0.0,
        δ⁺ = 0.0,

        #Thermal
        c = 50.0,
        h = 0.1,
        Tamb = 298.15,
        Temp = 298.15,

        #Kinetic
        k₀⁺ = 1e-4,
        k₀⁻ = 1e-1,
        #Initial
        x⁻₀ = 0.6,


        εₑˢ = 0.8,

        cₑ₀ = 1000.0,
        κ = 8.5e-2,
        t⁺ = 0.3,
        #Load
        input_type = 3,
        input_value = 4.2,
        E = 0.0,

        #SEI
        ω_+ = 0.0,
        ω_- = 0.0,
    )
end

function initial_conditions!(u,p,cellgeometry,initialcond,cathodeocv,anodeocv)
    V = initialcond["Starting Voltage[V]"]
    T = initialcond["Ambient Temperature[K]"]
    x⁻₀ = p.x⁻₀
    V⁺₀ = V + calcocv(anodeocv,x⁻₀,297.0)
    x⁺₀ = OCV.get_x_from_voltage(cathodeocv,V⁺₀,297.0)
    c_n_init = x⁻₀*(anodeocv.c_s_max-anodeocv.c_s_min)
    c_p_init = x⁺₀*(cathodeocv.c_s_max-cathodeocv.c_s_min)

    
    u[1:2] .= c_n_init
    u[3:5] .= p.cₑ₀
    u[6:7] .= c_p_init
    
    #=
    u[8] = p.εₛ⁻
    u[9] = p.εₑ⁻
    u[10] = p.εₑˢ
    u[11] = p.εₑ⁺
    u[12] = p.εₛ⁺

    u[13] = T
    =#

    return u
end

function cell_geometry()
    #Calculate Cell Volumes Based on Capacity
    T⁻ = 126e-6
    T⁺ = 125e-6
    Tˢ = 16e-6

    W⁻ = 800e-3
    W⁺ = 800e-3

    L⁺ = 70e-3
    L⁻ = 70e-3

    Vₛ⁻ = W⁻*L⁻*T⁻
    Vₛ⁺ = W⁺*L⁺*T⁺
    Vₑ⁻ = W⁻*L⁻*T⁻
    Vₑˢ = W⁺*L⁺*Tˢ
    Vₑ⁺ = W⁺*L⁺*T⁺

    

    cellgeometry = ComponentArray(
        Vₛ⁻ = W⁻*L⁻*T⁻/2,
        Vₛ⁺ = W⁺*L⁺*T⁺/2,
        Vₑ⁻ = W⁻*L⁻*T⁻,
        Vₑˢ = W⁺*L⁺*Tˢ,
        Vₑ⁺ = W⁺*L⁺*T⁺,
        T⁺ = T⁺,
        T⁻ = T⁻,
)

end

function initialize_airbus_ocv()
    A_p = [-123466.90208905171,
-65215.58255076725,
-77739.69488015345,
-132750.8136541972,
-109183.62667589705]

A_n = [-9.04729951400691e6,
-8.911003506471587e6,
-9.04657963355355e6,
-8.904669509837592e6,
-9.0363250622869e6,
-9.05878665345401e6,
-9.606335422964232e6,
-8.023042975317075e6,
-2.3190474522951595e6,
 1.4303914546788693e6]

U_0_p = 3.4271972387993173      # U_0_p
U_0_n = -46.09780594535385    # U_0_n


c_s_max⁻ = 50000.0
c_s_max⁺ = c_s_max⁻/1.15
c_s_min⁺ = 0.0
c_s_min⁻ = 0.0

cathode_ocv = OCV.RKPolynomial(A_p,U_0_p,c_s_max⁺,c_s_min⁺,length(A_p))
anode_ocv   = OCV.RKPolynomial(A_n,U_0_n,c_s_max⁻,c_s_min⁻,length(A_n))
return cathode_ocv,anode_ocv
end
