using CSV, DataFrames, Statistics, ProgressMeter, PythonPlot

FOLDERNAME = "/Users/abills/Datasets/cycle_individual_data/"
cells = [cell for cell in keys(ex_data_dict["cells"])]
@showprogress for cell in cells
fig,axes = subplots(figsize=(2.133,1.6))
s = []
cmap = PythonPlot.get_cmap("Blues_r")
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
    axes.plot(df.QDischargemAh./1000, df.EcellV,c=cmap(normalized_color))

end
norm = PythonPlot.matplotlib.colors.Normalize(vmin=0, vmax=maximum(ex_data_dict["cells"][cell]["cycles"]))
cbar = fig.colorbar(PythonPlot.matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes)
cbar.set_label("Cycle Number")
axes.set_xlabel("Capacity [Ah]")
axes.set_ylabel("Voltage [V]")
fig.savefig("figs/$(cell)_capplot.pdf",bbox_inches="tight")
fig.savefig("figs/$(cell)_capplot.png",bbox_inches="tight")
end

