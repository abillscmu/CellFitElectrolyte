using JLD2, Turing, CellFitElectrolyte, KDEDistributions, KernelDensity, CellFitElectrolyte.OCV

FOLDERNAME = "results/0302/"

cellgeometry = CellFitElectrolyte.cell_geometry()
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()


predicted_states = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]
derived_states = [:n_li, :εₛ⁻, :εₛ⁺]
predicted_states = vcat(predicted_states, derived_states)
data_dict = Dict()
for file in readdir(FOLDERNAME)
    #println(file)
    if file == "VAH25_353_HMC.jld2"
        continue
    end
    #if occursin(vah, file)
        vah = split(file, "_")[1]
        filename = FOLDERNAME*file
        chain = try
            d = load(filename)
            chain = d["chain"]
        catch
            @warn "problem with file $filename"
            continue
        end
        cycle_num = parse(Int, split(file, "_")[2])
        if vah in keys(data_dict)
            append!(data_dict[vah]["cycles"], cycle_num)
        else
            data_dict[vah] = Dict(
                "cycles" => [cycle_num,],
                "distributions" => Dict(),
                "mean" => Dict(),
                "time" => Dict()
            )
        end
        data_dict[vah]["distributions"][cycle_num] = Dict()
        data_dict[vah]["mean"][cycle_num] = Dict()
        data_dict[vah]["time"][cycle_num] = (chain.info.stop_time - chain.info.start_time)./3600
        for sym in predicted_states
            if !(sym in derived_states)
                distribution_data = chain[sym].data[:,1]
            elseif sym == :n_li
                εₑ⁻ = chain[:εₑ⁻].data[:,1]
                εₑ⁺ = chain[:εₑ⁺].data[:,1]
                frac_sol_am_pos = chain[:frac_sol_am_pos].data[:,1]
                frac_sol_am_neg = chain[:frac_sol_am_neg].data[:,1]
                c_e_0 = 2000
            
                V = 4.2
                x⁻₀ = 0.6
                V⁺₀ = V + calcocv(anodeocv,x⁻₀,297.0)
                x⁺₀ = OCV.get_x_from_voltage(cathodeocv,V⁺₀,297.0)
                c_n_init = x⁻₀*(anodeocv.c_s_max-anodeocv.c_s_min)
                c_p_init = x⁺₀*(cathodeocv.c_s_max-cathodeocv.c_s_min)
            
                #electrolyte
                n_li_elec = (εₑ⁻.*cellgeometry.Vₑ⁻ .+ εₑ⁺*cellgeometry.Vₑ⁺ .+ 0.8.*cellgeometry.Vₑˢ).*c_e_0
                n_li_pos = ((1 .- εₑ⁺) .* frac_sol_am_pos) .* cellgeometry.Vₑ⁺ .* c_p_init
                n_li_neg = ((1 .- εₑ⁻) .* frac_sol_am_neg) .* cellgeometry.Vₑ⁻  .* c_n_init
                n_li = n_li_pos .+ n_li_neg .+ n_li_elec
                distribution_data = n_li
            elseif sym == :εₛ⁻
                distribution_data = (1 .- chain[:εₑ⁻].data[:,1]).*chain[:frac_sol_am_neg].data[:,1]
            elseif sym == :εₛ⁺
                distribution_data = (1 .- chain[:εₑ⁺].data[:,1]).*chain[:frac_sol_am_pos].data[:,1]
            else
                error("symbol $sym not recognized")
            end
            distribution_kde = kde(distribution_data)
            distribution_bw = KernelDensity.default_bandwidth(distribution_data)
            distribution = KDEDistribution(distribution_data, distribution_kde, distribution_bw)
            data_dict[vah]["distributions"][cycle_num][sym] = distribution
            data_dict[vah]["mean"][cycle_num][sym] = mean(distribution_data)
        end
end