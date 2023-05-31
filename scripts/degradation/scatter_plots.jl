using PythonPlot, Statistics, PythonCall


xs = ["Max Temperature", "Max DOD", "Max Current", "Cycle Life","Mean Temperature", "Min Current"]
ys = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg, :n_li, :εₛ⁻, :εₛ⁺]

labeldict = Dict(
    "Max Temperature" => "Max Temperature [K]",
    "Max DOD" => "Max DOD [mAh]",
    "Max Current" => "Max Current [A]",
    "Cycle Life" => "Cycle Life",
    "Mean Temperature" => "Mean Temperature [K]",
    "Min Current" => "Min Current [A]"
)

nrow = length(ys)
ncol = length(xs)


x_type = "absolute"
y_type = "delta_N"

cells = [k for k in keys(data_dict)]


fig, axes = subplots(nrows=nrow,ncols=ncol);
fig.subplots_adjust(hspace=0.05, wspace=0.05)
for (j,x) in enumerate(xs)
    if x in keys(data_dict["VAH01"]["distributions"][2])
        first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][x] for cell in cells]
        last_val = [data_dict[cell]["mean"][maximum(data_dict[cell]["cycles"])][x] for cell in cells]
        if x_type == "first"
            x_axis = first_val
            xlabeltext = "Initial $x"
        elseif x_type == "last"
            x_axis = last_val
            xlabeltext = "Final $x"
        elseif x_type == "delta"
            x_axis = last_val .- first_val
            xlabeltext = "Δ$x"
        elseif x_type == "delta_N"
            cycle_life = [maximum(data_dict[cell]["cycles"]) for cell in cells]
            x_axis = (last_val .- first_val)./cycle_life
            xlabeltext = "Δ$x/Ν"
        end
    elseif x in keys(ex_data_dict["summary"]["VAH01"])
        x_axis = [ex_data_dict["summary"][cell][x] for cell in cells]
        xlabeltext = labeldict[x]
    elseif x == "Cycle Life"
        x_axis = [maximum(data_dict[cell]["cycles"]) for cell in cells]
        xlabeltext = x
    else
        error("x not found")
    end
    for (i,y) in enumerate(ys)
        if y in keys(data_dict["VAH01"]["distributions"][2])
            first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][y] for cell in cells]
            last_val = [data_dict[cell]["mean"][maximum(data_dict[cell]["cycles"])][y] for cell in cells]
            if y_type == "first"
                y_axis = first_val
                ylabeltext = "Initial $y"
            elseif y_type == "last"
                y_axis = last_val
                ylabeltext = "Final $y"
            elseif y_type == "delta"
                y_axis = last_val .- first_val
                ylabeltext = "Δ$y"
            elseif y_type == "delta_N"
                cycle_life = [maximum(data_dict[cell]["cycles"]) for cell in cells]
                y_axis = (last_val .- first_val)./cycle_life
                ylabeltext = "Δ$y/N"
            else
                error("y not found")
            end
        elseif y in keys(ex_data_dict["summary"]["VAH01"])
            y_axis = [ex_data_dict["summary"][cell][y] for cell in cells]
            ylabeltext = labeldict[y]
        elseif y == "Cycle Life"
            y_axis = [maximum(data_dict[cell]["cycles"]) for cell in cells]
            ylabeltext = y
        else
            error("y not found")
        end
        if false
            axes[i-1,j-1].hist(x_axis)
            if i != 1
                axes[i-1,j-1].set_xticks(pylist())
                axes[i-1,j-1].set_xlim(axes[i-2,j-1].get_xlim())
            end
            axes[i-1,j-1].set_yticks(pylist())
        elseif true
            axes[i-1,j-1].scatter(x_axis,y_axis,s=10)
            axes[i-1,j-1].xaxis.set_visible(false)
            axes[i-1,j-1].yaxis.set_visible(false)
            #axes[i-1,j-1].annotate("ρ = $C", (0.05, 0.95), xycoords="axes fraction", ha="left", va="top")
        else
            cmap = PythonPlot.matplotlib.cm.get_cmap("viridis")
            C = round(cor(x_axis, y_axis),digits=3)
            # make xaxis invisibel
            axes[i-1,j-1].xaxis.set_visible(false)
            # make spines (the box) invisible
            PythonPlot.setp(axes[i-1,j-1].spines.values(), visible=false)
            # remove ticks and labels for the left axis
            axes[i-1,j-1].tick_params(left=false, labelleft=false)
            #remove background patch (only needed for non-white background)
            axes[i-1,j-1].patch.set_visible(false)
            axes[i-1,j-1].scatter([0],[0],s=abs(C)*700, c=[C], cmap=cmap, vmin=-1,vmax=1)
            axes[i-1,j-1].set_yticks(pylist())
            axes[i-1,j-1].set_xticks(pylist())
            axes[i-1,j-1].annotate("ρ = $C", (0.05, 0.95), xycoords="axes fraction", ha="left", va="top")
        end

        # Set up ticks only on one side for the "edge" subplots...
        if j == 1
            axes[i-1,j-1].yaxis.set_ticks_position("left")
            axes[i-1,j-1].yaxis.set_visible(true)
            axes[i-1,j-1].set_ylabel(ylabeltext,rotation=0,ha="right")
        end
        if j == ncol
            axes[i-1,j-1].yaxis.set_ticks_position("right")
            axes[i-1,j-1].yaxis.set_visible(true)
            
        end
        if i == 1
            axes[i-1,j-1].xaxis.set_ticks_position("top")
            axes[i-1,j-1].xaxis.set_visible(true)
        end
        if i == nrow
            axes[i-1,j-1].xaxis.set_ticks_position("bottom")
            axes[i-1,j-1].xaxis.set_visible(true)
            axes[i-1,j-1].set_xlabel(xlabeltext,rotation=20, rotation_mode="anchor",ha="right",va="top")
        end
    end
    #fig.tight_layout()
end
fig.set_figwidth(8)
#fig.savefig("figs/correlation_xalone.pdf",bbox_inches="tight")
#fig.savefig("figs/correlation_xalone.png",bbox_inches="tight")
