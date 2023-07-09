using JLD2, PythonPlot

@load "files_to_use.jld2"

first_cells = ["VAH01", "VAH02","VAH05","VAH06","VAH07","VAH09","VAH10","VAH11","VAH12","VAH13","VAH15","VAH16","VAH17"]
extremes = ["VAH02", "VAH12", "VAH06", "VAH16", "VAH01", "VAH11", "VAH30", "VAH09"]
baselines = ["VAH01", "VAH17", "VAH27"]

params = (:ω, :εₑ⁺, :εₑ⁻, :frac_sol_am_neg, :frac_sol_am_pos)
@showprogress for cell in baselines
    for p in params
        fig, ax = subplots()
        xs = files_to_use[cell]
        ys = zeros(1000,length(xs))
        for (i,x) in enumerate(xs)
            d = load("results/0302/$(cell)_$(x)_HMC.jld2")
            chain = d["chain"]
            ys[:,i] .= chain[p].data[:,1]
        end
        ax.violinplot(ys, positions=xs,widths=50)
        fig.savefig("new_plots/$(cell)_$(p).pdf", bbox_inches = "tight")
    end
end