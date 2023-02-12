using PythonCall
using PythonPlot
using ProgressMeter
pygui(true)
pymannkendall = pyimport("pymannkendall")


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
            pmk = pymannkendall.original_test(medians)
            j_result = pyconvert(String, pmk.trend)
            if j_result == "increasing"
                b = 1
            elseif j_result == "decreasing"
                b = -1
            elseif j_result == "no trend"
                b = 0
            end
            append!(cycs, pyconvert(Float64,pmk.Tau))
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
    ax.set_yticks([0,1,2,3,4],labels=String.([:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]))
    fig.tight_layout()
    cbar = fig.colorbar(im,ticks=[-1, 0, 1],ax=ax)
    cbar.ax.set_yticks([-1, 0, 1],labels=["Decreasing","No Trend","Increasing"])
    ax.set_title(cell)
    #fig.savefig("$cell.png")
end