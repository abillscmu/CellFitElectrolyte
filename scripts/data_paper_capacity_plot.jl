using CSV, DataFrames, Statistics, ProgressMeter

FOLDERNAME = "/Users/abills/Datasets/cycle_individual_data/"
markers = ["o","v","^","<",">","1","2","3","4","s","p","P","*","h","H","+","x","X","D","d","|","_"]
cell = "VAH01"
cap_data_dict = Dict()
fig,axes = subplots()
cells = [cell for cell in keys(ex_data_dict["cells"])]
for (i,cell) in enumerate(cells)
    println(cell)

cmap = PythonPlot.get_cmap("Blues_r")
cap_data_dict[cell] = Dict(
    "capacity"=>Float64[],
    "cycle" => Int[]
)
for file in readdir(FOLDERNAME)
    if file == ".DS_Store"
        continue
    end
    arr = split(file, ['_','.'])
    vah = arr[1]
    if vah != cell
        continue
    end
    cycle = parse(Int,arr[2])
    df = CSV.read(FOLDERNAME*file, DataFrame)
    if df.Ns[1] != 3
        continue
    end
    cycle_life = maximum(ex_data_dict["cells"][cell]["cycles"]) + 200
    normalized_color = (cycle+200)/(cycle_life+200)
    df = filter(row-> row.Ns in [3,5,7], df)
    df = filter(row-> row.ImA < -300, df)
    append!(cap_data_dict[cell]["capacity"],df.QDischargemAh[end]./1000)
    append!(cap_data_dict[cell]["cycle"],cycle)

end
axes.scatter(cap_data_dict[cell]["cycle"],cap_data_dict[cell]["capacity"],marker=markers[i])
axes.set_xlabel("Cycle")
axes.set_ylabel("Capacity [Ah]")
end
axes.legend(cells,loc="center left",ncols=2,bbox_to_anchor=(1.02, 0.5))
fig.savefig("figs/capacity_plot.pdf",bbox_inches="tight")
fig.savefig("figs/capacity_plot.png",bbox_inches="tight")