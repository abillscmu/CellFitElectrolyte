using DataFrames, PythonPlot, CSV, PythonCall, PythonPlot, ProgressMeter, LinearAlgebra
deg_df = CSV.read("data/degstuff.csv",DataFrame)
#pygui(true)
pymannkendall = pyimport("pymannkendall")
mpl = PythonPlot.matplotlib
TF(z,N) = sqrt(4N+10)*z/(3*sqrt(N-1))

sym = :ω
cells = unique([split(f,"_")[1] for f in readdir("results/outputs0117_elec/")])
#cells = ["VAH01"]
celldict = Dict()
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
    increasing_mask = total_cycs .> 1.96
    decreasing_mask = total_cycs .< -1.96
    total_mask = increasing_mask .- decreasing_mask

    diagnosis_dict = Dict(
        name => Dict(
            "diagnosis" => zeros(5:max_cycle),
            "symptoms" => Array(filter(:Name => (Name) -> Name == name,deg_df))[2:end]
        ) for name in deg_df[!,"Name"]
    )
    diagnoses = zeros(max_cycle, length(deg_df[!,"Name"]))
    for j in 1:length(deg_df[!,"Name"])
        additive = 0.015*(j-4)
        diagnoses[:,j] .+= additive
    end
    for cycle in 5:max_cycle
        ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos, n_li, εₛ⁺, εₛ⁻ = total_mask[cycle-4,:]
        comparison_vec = [εₑ⁻,εₑ⁺,ω,n_li,εₛ⁺,εₛ⁻]
        for (j,disease) in enumerate(deg_df[!,"Name"])
            distance = dot(comparison_vec,diagnosis_dict[disease]["symptoms"])./(norm(comparison_vec)*norm(diagnosis_dict[disease]["symptoms"]))
            if isnan(distance)
                continue
            end
            diagnoses[cycle, j] += distance
        end
    end
    cycles = 5:max_cycle
    sig_tau_vec = [1.96 for N in 5:max_cycle]
    nrows=8
    figh = 0.35 + 0.15 + (nrows + (nrows - 1) * 0.1) * 0.22
    fig, ax = subplots(nrows=nrows + 1, figsize=(8.5, figh))
    fig.subplots_adjust(top=1 - 0.35 / figh, bottom=0.15 / figh,
                        left=0.2, right=0.99)
    total_cycs = Array{Float64}(total_cycs)
    for j in 1:size(diagnoses)[2]
        ax[j-1].imshow(Array(diagnoses[5:end,j]'),aspect="auto",interpolation="none",cmap=mpl.colormaps["Greys"], vmin=0, vmax=1)
        ax[j-1].text(-0.01, 0.0, deg_df.Name[j]*"  ", va="center", ha="right", fontsize=16)
        ax[j-1].axes.get_yaxis().set_visible(false)
        ax[j-1].axes.get_xaxis().set_visible(false)
    end
    #ax.legend(deg_df[!,"Name"],loc="center left",bbox_to_anchor=(1,0.5))
    #box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    #ax[0].set_title(cell)
    ax[-1].set_xlabel("Cycles Seen", fontsize=16)
    ax[-1].axes.get_xaxis().set_visible(true)
    ax[-1].tick_params(labelsize=16)
    #fig.tight_layout()
    fig.savefig("figs/diagnosis/$(cell)_diagnosis.png",bbox_inches="tight")
    fig.savefig("figs/diagnosis/$(cell)_diagnosis.pdf",bbox_inches="tight")
    #fig.savefig("$cell.png")
    celldict[cell] = diagnoses
end