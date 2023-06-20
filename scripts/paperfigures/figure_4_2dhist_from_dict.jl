using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall
np = pyimport("numpy")
#pygui(true)

PythonPlot.matplotlib.rcParams["font.size"] = 10


FOLDERNAME = "results/outputs0117_elec/"
CELLS = ["VAH01", "VAH02", "VAH05", "VAH06", "VAH12", "VAH30"]
#CELLS = ["VAH01", "VAH02"]
SYMBOLS = SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg, :n_li, :εₛ⁻, :εₛ⁺]
num_rows = length(SYMBOLS)
num_cols = length(CELLS)

ylimits = Dict(
    :ω => pylist([0.01,0.05]),
    :εₑ⁻ => pylist([0.05,0.5]),
    :εₑ⁺ => pylist([0.05, 0.5]),
    :frac_sol_am_pos => pylist([0.5,1.0]),
    :frac_sol_am_neg => pylist([0.5,1.0]),
    :n_li => pylist([0.05, 0.2  ]),
    :εₛ⁺ => pylist([0.25, 0.95]),
    :εₛ⁻ => pylist([0.25, 0.95])
)
ylabels = Dict(
    :ω => "ω [Ω]",
    :εₑ⁻ => "εₑ⁻",
    :εₑ⁺ => "εₑ⁺",
    :frac_sol_am_pos => "fₛ⁺",
    :frac_sol_am_neg => "fₛ⁻",
    :n_li => "nₗᵢ [mol]",
    :εₛ⁺ => "εₛ⁺",
    :εₛ⁻ => "εₛ⁻"
)


fig, axes = subplots(nrows=num_rows, ncols=num_cols, figsize=(8.5,6))
fig.subplots_adjust(hspace=0.05, wspace=0.05)

  
for (j,CELL) in enumerate(CELLS)
    max_cyc = maximum(data_dict[CELL]["cycles"])
    min_cyc = minimum(data_dict[CELL]["cycles"])
    xlimit_thiscell = pylist([min_cyc, max_cyc])
    for (i, SYMBOL) in enumerate(SYMBOLS)
        ylimit_thissym = ylimits[SYMBOL]
        yname = ylabels[SYMBOL]
        val = []
        cyc = []
        for cell in 0:2500
        try
            data = data_dict[CELL]["distributions"][cell][SYMBOL].data
            if length(data) != 1000
                continue
            else
                append!(val, data)
                append!(cyc, cell.*ones(1000))
            end  
        catch
            @warn "problem with $cell"
            continue
        end
    end

    axes[i-1,j-1].hexbin(cyc, val, gridsize=100, extent=pytuple((xlimit_thiscell[0], xlimit_thiscell[1], ylimit_thissym[0], ylimit_thissym[1])), mincnt = 10)
    #colorbar(label="Density")
    #grid(0.2)
    if i == num_rows
        axes[i-1,j-1].set_xlabel("Cycle")
        max_x = pyconvert(Int64,xlimit_thiscell[1])
        min_x = pyconvert(Int64,xlimit_thiscell[0])
        x_range = max_x - min_x
        top_mark = round(0.75*x_range + min_x)
        bot_mark = round(0.25*x_range + min_x)
        axes[i-1,j-1].set_xlim(xlimit_thiscell)
        axes[i-1,j-1].set_xticks(pylist([bot_mark, top_mark]))
        axes[i-1,j-1].set_xticklabels(pylist([bot_mark, top_mark]),rotation=45, ha="right",rotation_mode="anchor")
    else
        max_x = pyconvert(Int64,xlimit_thiscell[1])
        min_x = pyconvert(Int64,xlimit_thiscell[0])
        x_range = max_x - min_x
        top_mark = round(0.75*x_range + min_x)
        bot_mark = round(0.25*x_range + min_x)
        axes[i-1,j-1].set_xlim(xlimit_thiscell)
        axes[i-1,j-1].set_xticks(pylist([bot_mark, top_mark]))
        axes[i-1,j-1].set_xticklabels([])
        axes[i-1,j-1].xaxis.set_tick_params(bottom=false)
    end
    if i > 1
        axes[i-1,j-1].set_xlim(axes[i-2,j-1].get_xlim())
    end
    if i == 1
        axes[i-1,j-1].set_title(CELL)#,rotation=-45,rotation_mode="anchor",ha="right")
    end
    if j == 1
        axes[i-1,j-1].set_ylabel(yname)
        max_y = pyconvert(Float64,ylimit_thissym[1])
        min_y = pyconvert(Float64,ylimit_thissym[0])
        y_range = max_y - min_y
        top_mark = round(0.75*y_range + min_y,digits=3)
        bot_mark = round(0.25*y_range + min_y,digits=3)
        axes[i-1,j-1].set_ylim(ylimit_thissym)
        axes[i-1,j-1].set_yticks(pylist([bot_mark, top_mark]))
    else
        max_y = pyconvert(Float64,ylimit_thissym[1])
        min_y = pyconvert(Float64,ylimit_thissym[0])
        y_range = max_y - min_y
        top_mark = round(0.75*y_range + min_y,digits=3)
        bot_mark = round(0.25*y_range + min_y,digits=3)
        axes[i-1,j-1].set_yticks(pylist([bot_mark, top_mark]))
        axes[i-1,j-1].set_yticklabels([])
        axes[i-1,j-1].yaxis.set_tick_params(bottom=false)
        axes[i-1,j-1].set_ylim(axes[i-1,j-2].get_ylim())
    end
    axes[i-1,j-1].grid(alpha=0.2, which="major")
    #title("Series: $CELL")
end
end 
#tight_layout()
norm = PythonPlot.matplotlib.colors.Normalize(vmin=0,vmax=2)
cmap = PythonPlot.matplotlib.pyplot.get_cmap("viridis")
sm = PythonPlot.matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array(np.array([0.0]))
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])

cb = fig.colorbar(sm, cax=cbar_ax,ticks=[0,1,2],boundaries=np.arange(-0.05,2.1,.1))
cb.ax.set_yticklabels(["Low\nDensity", "Medium\nDensity", "High\nDensity"], fontsize=10)
savefig("figs/real/2dhist_supp.png", bbox_inches="tight")
savefig("figs/real/2dhist_supp.pdf", bbox_inches="tight")