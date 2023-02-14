using JLD2, PythonPlot, PythonCall, LinearAlgebra, Statistics
np = pyimport("numpy")


cell_errors = Dict()
cell_error_vec = Float64[]
strings = []
for file in readdir("results/errors_0117")
    cell = split(file,"_")[1]
    f = split(file, "_errors.jld2")[1]
    if !(cell in keys(cell_errors))
        cell_errors[cell] = []
    end
    try
        d = load("results/errors_0117/$file")
        rmse = d["rmse"]
        append!(cell_errors[cell], rmse)
        append!(cell_error_vec, rmse)
        for n in 1:1000
            push!(strings, string(f))
        end
    catch
        continue
    end
end

num_rows = 6
num_cols = 4
fig, axes = subplots(num_rows,num_cols, figsize=(8,10))
fig2, axes2 = subplots(num_rows,num_cols, figsize=(8,10))
for (k,cell) in enumerate(keys(cell_errors))
    i = Int(floor((k-1)/num_cols)) + 1
    j = k%num_cols
    if j == 0
        j = 4
    end
    cell_error_vec = cell_errors[cell]
    top_1 = quantile(cell_error_vec, 0.99)
    cell_error_cut = cell_error_vec[cell_error_vec .< top_1]
    med = median(cell_error_vec)
    PythonPlot.matplotlib.rcParams["font.size"] = 10
    #fig, axes = subplots()
    bins = 10 .^range(-3,stop=1,length=1000)
    axes[i-1,j-1].hist(cell_error_cut, bins=bins,density=true)
    axes[i-1,j-1].axvline(x=med, color="k", linestyle="--")
    axes[i-1,j-1].set_xscale("log")
    axes[i-1,j-1].grid(which="both",alpha=0.3)
    axes[i-1,j-1].set_xlabel("Voltage RMSE [V]")
    axes[i-1,j-1].set_ylabel("Density")
    axes[i-1,j-1].set_title(cell)
    #fig.savefig("figs/si/error_distribution/$(cell)_voltage_error.pdf",bbox_inches="tight")
    #fig.savefig("figs/si/error_distribution/$(cell)_voltage_error.png",bbox_inches="tight")

    #fig2, axes2 = subplots()
    bins = 10 .^range(-3,stop=1,length=1000)
    axes2[i-1,j-1].hist(cell_error_cut, bins=bins,cumulative=true, density=true)
    axes2[i-1,j-1].set_xscale("log")
    axes2[i-1,j-1].grid(which="both",alpha=0.3)
    axes2[i-1,j-1].set_xlabel("Voltage RMSE [V]")
    axes2[i-1,j-1].set_ylabel("Cumulative Density")
    axes2[i-1,j-1].axvline(x=med, color="k", linestyle="--")
    axes2[i-1,j-1].set_title(cell)
    #fig2.savefig("figs/si/error_distribution/$(cell)_voltage_error_cdf.pdf",bbox_inches="tight")
    #fig2.savefig("figs/si/error_distribution/$(cell)_voltage_error_cdf.png",bbox_inches="tight")
end
axes[num_rows-1,2].axis("off")
axes[num_rows-1,3].axis("off")
axes2[num_rows-1,2].axis("off")
axes2[num_rows-1,3].axis("off")
fig.tight_layout()
fig2.tight_layout()
fig.savefig("figs/si/error_distribution/hist.pdf",bbox_inches="tight")
fig.savefig("figs/si/error_distribution/hist.png",bbox_inches="tight")
fig2.savefig("figs/si/error_distribution/cdf.pdf",bbox_inches="tight")
fig2.savefig("figs/si/error_distribution/cdf.png",bbox_inches="tight")

