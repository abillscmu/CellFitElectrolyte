using JLD2, PythonPlot, PythonCall, LinearAlgebra, Statistics
np = pyimport("numpy")

individual = false
recompute = true
plot_med_cycle = true

if recompute
    cell_errors = Dict()
    cell_error_vec = Float64[]
    cells = []
    strings = []
    for file in readdir("results/0302_errors")
        cell = split(file,"_")[1]
        f = split(file, "_HMC.jld2")[1]
        if !(cell in keys(cell_errors))
            cell_errors[cell] = []
        end
        try
            d = load("results/0302_errors/$file")
            rmse = d["rmse"]
            vah = split(file, "_errors")[1]
            append!(cell_errors[cell], rmse)
            append!(cell_error_vec, rmse)
            for r in rmse
                push!(cells, vah)
            end
            for n in 1:1000
                push!(strings, string(f))
            end
        catch
            continue
        end
    end
end




if individual
    num_rows = 6
    num_cols = 4
    fig, ax = subplots(num_rows,num_cols, figsize=(8,10))
    fig2, ax2 = subplots(num_rows,num_cols, figsize=(8,10))
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
        #fig, ax = subplots()
        bins = 10 .^range(-3,stop=1,length=1000)
        ax[i-1,j-1].hist(cell_error_cut, bins=bins,density=true)
        ax[i-1,j-1].axvline(x=med, color="k", linestyle="--")
        ax[i-1,j-1].set_xscale("log")
        ax[i-1,j-1].grid(which="both",alpha=0.3)
        ax[i-1,j-1].set_xlabel("Voltage RMSE [V]")
        ax[i-1,j-1].set_ylabel("Density")
        ax[i-1,j-1].set_title(cell)
        #fig.savefig("figs/si/error_distribution/$(cell)_voltage_error.pdf",bbox_inches="tight")
        #fig.savefig("figs/si/error_distribution/$(cell)_voltage_error.png",bbox_inches="tight")

        #fig2, ax2 = subplots()
        bins = 10 .^range(-3,stop=1,length=1000)
        ax2[i-1,j-1].hist(cell_error_cut, bins=bins,cumulative=true, density=true)
        ax2[i-1,j-1].set_xscale("log")
        ax2[i-1,j-1].grid(which="both",alpha=0.3)
        ax2[i-1,j-1].set_xlabel("Voltage RMSE [V]")
        ax2[i-1,j-1].set_ylabel("Cumulative Density")
        ax2[i-1,j-1].axvline(x=med, color="k", linestyle="--")
        ax2[i-1,j-1].set_title(cell)
        #fig2.savefig("figs/si/error_distribution/$(cell)_voltage_error_cdf.pdf",bbox_inches="tight")
        #fig2.savefig("figs/si/error_distribution/$(cell)_voltage_error_cdf.png",bbox_inches="tight")
    end
    ax[num_rows-1,2].axis("off")
    ax[num_rows-1,3].axis("off")
    ax2[num_rows-1,2].axis("off")
    ax2[num_rows-1,3].axis("off")
    fig.tight_layout()
    fig2.tight_layout()
    fig.savefig("figs/si/error_distribution/hist.pdf",bbox_inches="tight")
    fig.savefig("figs/si/error_distribution/hist.png",bbox_inches="tight")
    fig2.savefig("figs/si/error_distribution/cdf.pdf",bbox_inches="tight")
    fig2.savefig("figs/si/error_distribution/cdf.png",bbox_inches="tight")
else
    fig, ax = subplots()
    fig2, ax2 = subplots()
    top_1 = quantile(cell_error_vec, 0.99)
    cell_error_cut = cell_error_vec#[cell_error_vec .< top_1]
    med = median(cell_error_vec)
    PythonPlot.matplotlib.rcParams["font.size"] = 16
    #fig, ax = subplots()
    bins = 10 .^range(-3,stop=1,length=1000)
    ax.hist(cell_error_cut, bins=bins,density=true)
    ax.axvline(x=med, color="k", linestyle="--")
    ax.set_xscale("log")
    ax.grid(which="both",alpha=0.3)
    ax.set_xlabel("Voltage RMSE [V]")
    ax.set_ylabel("Density")
    #fig.savefig("figs/si/error_distribution/$(cell)_voltage_error.pdf",bbox_inches="tight")
    #fig.savefig("figs/si/error_distribution/$(cell)_voltage_error.png",bbox_inches="tight")

    fig2, ax2 = subplots()
    bins = 10 .^range(-3,stop=1,length=1000)
    ax2.hist(cell_error_cut, bins=bins,cumulative=true, density=true)
    ax2.set_xscale("log")
    ax2.grid(which="both",alpha=0.3)
    ax2.set_xlabel("Voltage RMSE [V]")
    ax2.set_ylabel("Cumulative Density")
    ax2.axvline(x=med, color="k", linestyle="--")
    #fig2.savefig("figs/si/error_distribution/$(cell)_voltage_error_cdf.pdf",bbox_inches="tight")
    #fig2.savefig("figs/si/error_distribution/$(cell)_voltage_error_cdf.png",bbox_inches="tight")
end
