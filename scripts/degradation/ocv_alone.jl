using CellFitElectrolyte, JLD2, Turing, OCV, NLsolve


cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
cellgeometry = CellFitElectrolyte.cell_geometry()

capacity_dict = Dict()



for file in ["VAH01_840_HMC.jld2"]
    println(file)
    cell = split(file,"_")[1]
    cycle = parse(Int,split(file,"_")[2])
    if cell != "VAH01"
        continue
    end
    chain = try
        d = load("results/outputs0117_elec/$file")
        chain = d["chain"]
    catch
        @warn "problem with file $file"
        continue
    end

    if !(cell in keys(capacity_dict))
        capacity_dict[cell] = Dict()
    end
    εₑ⁺ = chain[:εₑ⁺].data[:,1]
    εₑ⁻ = chain[:εₑ⁻].data[:,1]
    fₛ⁺ = chain[:frac_sol_am_pos].data[:,1]
    fₛ⁻ = chain[:frac_sol_am_neg].data[:,1]
    εₛ⁺ = (1 .- εₑ⁺).*fₛ⁺
    εₛ⁻ = (1 .- εₑ⁻).*fₛ⁻

    Vₛ⁺ = 2 .* cellgeometry.Vₛ⁺ .* εₛ⁺
    Vₛ⁻ = 2 .* cellgeometry.Vₛ⁻ .* εₛ⁻

    x⁻₀ = 0.6
    V⁺₀ = 4.2 + calcocv(anodeocv,x⁻₀,297.0)
    x⁺₀ = OCV.get_x_from_voltage(cathodeocv,V⁺₀,297.0)
    c_n_init = x⁻₀*(anodeocv.c_s_max-anodeocv.c_s_min)
    c_p_init = x⁺₀*(cathodeocv.c_s_max-cathodeocv.c_s_min)



    function calc_ocv(capacity, Vₛ⁻, Vₛ⁺)
        capacity = capacity[1]
        n_li_n_begin = c_n_init*Vₛ⁻
        n_li_p_begin = c_p_init*Vₛ⁺

        n_li_n_end = n_li_n_begin - capacity
        n_li_p_end = n_li_p_begin + capacity
    
        x⁺ = (n_li_p_end/Vₛ⁺)/cathodeocv.c_s_max
        x⁻ = (n_li_n_end/Vₛ⁻)/anodeocv.c_s_max
        try
            V = calcocv(cathodeocv, x⁺, 297.0) - calcocv(anodeocv, x⁻, 297.0)
            return [V - 3.0]
        catch
            return [1e12]
        end
    end

    capacities = zeros(length(chain))

    for i in 1:length(chain)
        Vsn = Vₛ⁻[i]
        Vsp = Vₛ⁺[i]
        f = (capacity) -> calc_ocv(capacity, Vsn, Vsp)
        cap = nlsolve(f, [3/26.8])
        capacities[i] = cap.zero[1]
    end
    capacity_dict[cell][cycle] = capacities
end
