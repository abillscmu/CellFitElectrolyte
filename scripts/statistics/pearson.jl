using PythonPlot
using ProgressMeter
using Statistics
pygui(true)


sym = :ω
cells = unique([split(f,"_")[1] for f in readdir("results/outputs0117_elec/")])
cells = ["VAH01",]
for cell in cells
max_cycle = maximum(data_dict[cell]["cycles"])
num_cycs = length(data_dict[cell]["cycles"])
println(cell)
total_cycs = []
    for sym in [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]
        cycs = []
        for last_cycle in 5:max_cycle
            cycles = sort(data_dict[cell]["cycles"][data_dict[cell]["cycles"].<last_cycle])
            medians = [median(data_dict[cell]["distributions"][cyc][sym].data) for cyc in cycles]
            r = Statistics.cor(cycles,medians)
            append!(cycs, r)
        end
        if sym == :ω
            total_cycs = cycs
        else
            total_cycs = hcat(total_cycs,cycs)
        end
    end

    total_cycs = Array{Float64}(total_cycs)
    fig, ax = subplots()
    im = ax.plot(total_cycs)
    fig.tight_layout()
    ax.set_title(cell)
    #fig.savefig("$cell.png")
end