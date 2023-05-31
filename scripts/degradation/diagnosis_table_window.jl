using DataFrames, PythonPlot, CSV, PythonCall, PythonPlot, ProgressMeter
deg_df = CSV.read("data/degstuff.csv",DataFrame)
pygui(true)
pymannkendall = pyimport("pymannkendall")
T(z,N) = sqrt(4N+10)*z/(3*sqrt(N-1))

sym = :ω
cells = unique([split(f,"_")[1] for f in readdir("results/outputs0117_elec/")])
cells = ["VAH01"]
for cell in cells
max_cycle = maximum(data_dict[cell]["cycles"])
num_cycs = length(data_dict[cell]["cycles"])
println(cell)
total_cycs = []
    for sym in [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos, :n_li, :εₛ⁺, :εₛ⁻]
        cycs = []
        println(sym)
        
        for last_cycle in 45:max_cycle
            first_cycle = last_cycle - 40
            cycles = sort((data_dict[cell]["cycles"][(data_dict[cell]["cycles"].<last_cycle).&(data_dict[cell]["cycles"].>first_cycle)]))
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
    println(size(total_mask))

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
    for cycle in 45:max_cycle
        ω, εₑ⁻, εₑ⁺, frac_sol_am_neg, frac_sol_am_pos, n_li, εₛ⁺, εₛ⁻ = total_mask[cycle-44,:]
        comparison_vec = [εₑ⁻,εₑ⁺,ω,n_li,εₛ⁺,εₛ⁻]
        for (j,disease) in enumerate(deg_df[!,"Name"])
            fraction = sum((comparison_vec .== diagnosis_dict[disease]["symptoms"]).*abs.(comparison_vec))./sum(abs.(diagnosis_dict[disease]["symptoms"]))
            diagnoses[cycle, j] += fraction
        end
    end

    cycles = 5:max_cycle
    sig_tau_vec = [1.96 for N in 5:max_cycle]
    fig = figure(1, (8.5,2));
    fig.clf()
    ax = subplot(111)
    total_cycs = Array{Float64}(total_cycs)
    ax.imshow(Array(diagnoses[5:end,:]'), aspect=max_cycle/(8*1.5))
    #ax.legend(deg_df[!,"Name"],loc="center left",bbox_to_anchor=(1,0.5))
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    title(cell)
    xlabel("Cycles Seen")
    ylabel("Fraction of Indications")
    grid()
    savefig("figs/diagnosis/$(cell)_diagnosis.png",bbox_inches="tight")
    savefig("figs/diagnosis/$(cell)_diagnosis.pdf",bbox_inches="tight")
    #fig.savefig("$cell.png")
end