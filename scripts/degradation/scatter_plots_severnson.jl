using PythonPlot, Statistics, PythonCall




cycs = [50,100,150,200,250]
ys = repeat(["Cycle Life"],length(cycs))
xs = [k for k in keys(data_dict["VAH01"]["distributions"][2])]

nrow = length(ys)
ncol = length(xs)


x_type = "delta_N"
y_type = "delta_N"

cells = [k for k in keys(data_dict)]


fig, axes = subplots(nrows=nrow,ncols=ncol,figsize=(16,8));
fig.subplots_adjust(hspace=0.05, wspace=0.05)
for (j,x) in enumerate(xs)
    y_axis = [maximum(data_dict[cell]["cycles"]) for cell in cells]
    for (i,y) in enumerate(ys)
    ylabeltext = "Cycle Life"
    if x in keys(data_dict["VAH01"]["distributions"][2])
        first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][x] for cell in cells]
        last_val = [data_dict[cell]["mean"][sort(data_dict[cell]["cycles"])[searchsortedfirst(sort(data_dict[cell]["cycles"]),cycs[i]-5)]][x] for cell in cells]
        for lv in cycs[i]-9:cycs[i]+10
            last_val .+= [data_dict[cell]["mean"][sort(data_dict[cell]["cycles"])[searchsortedfirst(sort(data_dict[cell]["cycles"]),lv)]][x] for cell in cells]
        end
        last_val = last_val./21
        
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
        elseif x_type == "variance"
            x_axis = zeros(length(cells))
            for (k,cell) in enumerate(cells)
                cycles_we_have = data_dict[cell]["cycles"][data_dict[cell]["cycles"].<=cycs[i]]
                x_axis[k] = log10(var([data_dict[cell]["mean"][cyc][x] for cyc in cycles_we_have]))
            end
            xlabeltext = "var($x)"
        else
            error("x not found")
        end
    elseif x in keys(ex_data_dict["summary"]["VAH01"])
        x_axis = [ex_data_dict["summary"][cell][x] for cell in cells]
        xlabeltext = x
    elseif x == "Cycle Life"
        x_axis = [maximum(data_dict[cell]["cycles"]) for cell in cells]
        xlabeltext = x
    else
        error("x not found")
    end
        ylabeltext = ylabeltext*"(1-$(cycs[i]))"
        println(size(x_axis))
        println(size(y_axis))
        axes[i-1,j-1].scatter(x_axis,y_axis,s=10)
        C = round(cor(x_axis, y_axis),digits=3)
        axes[i-1,j-1].xaxis.set_visible(false)
        axes[i-1,j-1].yaxis.set_visible(false)

        # Set up ticks only on one side for the "edge" subplots...
        if j == 1
            axes[i-1,j-1].yaxis.set_ticks_position("left")
            axes[i-1,j-1].yaxis.set_visible(true)
            axes[i-1,j-1].set_ylabel(ylabeltext,rotation=45,ha="right",fontdict=pydict(Dict("size"=>14)))
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
            axes[i-1,j-1].set_xlabel(xlabeltext,fontdict=pydict(Dict("size"=>14)))
        end
        axes[i-1,j-1].annotate("ρ = $C", (0.05, 0.95), xycoords="axes fraction", ha="left", va="top")
    end
    #fig.tight_layout()
end

fig.savefig("figs/correlation_cycle_life.pdf")
fig.savefig("figs/correlation_cycle_life.png")
