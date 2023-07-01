using PythonCall
using PythonPlot
using ProgressMeter
#pygui(true)
pymannkendall = pyimport("pymannkendall")
#T(z,N) = sqrt(4N+10)*z/(3*sqrt(N-1))

sym = :ω
cells = unique([split(f,"_")[1] for f in readdir("results/outputs0117_elec/")])
zscoredict = Dict()
cells=["VAH01"]
for cell in cells
max_cycle = maximum(data_dict[cell]["cycles"])
num_cycs = length(data_dict[cell]["cycles"])
println(cell)
total_cycs = []
    for sym in [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos, :n_li, :εₛ⁺, :εₛ⁻]
        cycs = []
        println(sym)
        
        for last_cycle in 5:max_cycle
            cycles = sort(data_dict[cell]["cycles"][data_dict[cell]["cycles"].<last_cycle])
            medians = [median(data_dict[cell]["distributions"][cyc][sym].data) for cyc in cycles]
            pmk = pymannkendall.original_test(medians)
            j_result = pyconvert(String, pmk.trend)
            if j_result == "increasing"
                b = 1
            elseif j_result == "decreasing"
                b = -1
            elseif j_result == "no trend"
                b = 0
            end
            append!(cycs, pyconvert(Float64,pmk.z))
            append!(cycles, last_cycle)
        end
        if sym == :ω
            total_cycs = cycs
        else
            total_cycs = hcat(total_cycs,cycs)
        end
    end
    cycles = 5:max_cycle
    sig_tau_vec = [1.96 for N in 5:max_cycle]
    fig, ax = subplots(figsize=(5,2.5))
    total_cycs = Array{Float64}(total_cycs)
    tot_min = minimum(total_cycs)
    tot_max = maximum(total_cycs)
    ps = ax.plot(cycles, total_cycs)
    ax.plot(cycles,sig_tau_vec,"k--")
    ax.plot(cycles, -sig_tau_vec,"k--")
    ax.legend(String.([:ω, :εₑ⁻, :εₑ⁺, "fₛ⁻", "fₛ⁺","nₗᵢ", :εₛ⁻, :εₛ⁺,"No Trend p=0.05"]),loc="lower center",bbox_to_anchor=(0.5, 1.05), fontsize=14, ncol=3)
    for i in 1:(size(total_cycs)[1] - 1)
        for (j,sym) in enumerate([:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos,:n_li, :εₛ⁻, :εₛ⁺])
            this = abs(total_cycs[i,j])
            next = abs(total_cycs[i+1,j])
            a = this>1.96
            c = next>1.96
            if xor(a,c)
                plot([cycles[i],cycles[i]],[tot_min,tot_max],ps[j-1].get_color(),alpha=0.5)
            end
        end
    end
    ylim([tot_min,tot_max])
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.tick_params(labelsize=14)
    zscoredict[cell] = Dict()
    zscoredict[cell]["cycles"] = cycles
    zscoredict[cell]["mk"] = total_cycs
    zscoredict[cell]["sig_tau"] = sig_tau_vec
    #title(cell)
    xlabel("Cycles Seen",fontsize=15)
    ylabel("Mann-Kendall Z Score",fontsize=15)
    grid()
    savefig("figs/mk/$(cell)_zscore.png",bbox_inches="tight")
    savefig("figs/mk/$(cell)_zscore.pdf",bbox_inches="tight")
    #fig.savefig("$cell.png")
end