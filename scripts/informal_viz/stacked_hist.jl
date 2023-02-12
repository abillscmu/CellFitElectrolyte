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

top_1 = quantile(cell_error_vec, 0.99)


cell_error_cut = cell_error_vec[cell_error_vec .< top_1]
med = median(cell_error_vec)
PythonPlot.matplotlib.rcParams["font.size"] = 16
fig, axes = subplots()
bins = 10 .^range(-3,stop=1,length=1000)
axes.hist(cell_error_cut, bins=bins,density=true)
axes.axvline(x=med, color="k", linestyle="--")
axes.set_xscale("log")
axes.grid(which="both",alpha=0.3)
axes.set_xlabel("Voltage RMSE [V]")
axes.set_ylabel("Density")
fig.savefig("figs/voltage_error.pdf",bbox_inches="tight")
fig.savefig("figs/voltage_error.png",bbox_inches="tight")

fig2, axes2 = subplots()
bins = 10 .^range(-3,stop=1,length=1000)
axes2.hist(cell_error_cut, bins=bins,cumulative=true, density=true)
axes2.set_xscale("log")
axes2.grid(which="both",alpha=0.3)
axes2.set_xlabel("Voltage RMSE [V]")
axes2.set_ylabel("Cumulative Density")
axes2.set_yticks(0:0.1:1)
axes2.axvline(x=med, color="k", linestyle="--")
fig2.savefig("figs/voltage_error_cdf.pdf",bbox_inches="tight")
fig2.savefig("figs/voltage_error_cdf.png",bbox_inches="tight")

