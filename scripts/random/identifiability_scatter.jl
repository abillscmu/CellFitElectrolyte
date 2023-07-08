using PythonPlot, Statistics, PythonCall, KernelDensity
np = pyimport("numpy")


SYMBOLS = SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg, :n_li, :εₛ⁻, :εₛ⁺]
xs = ys = SYMBOLS

mycolors = ["blue","orange","green"]

nrow = length(ys)
ncol = length(xs)


cells = keys(data_dict)
for CELL in cells
println(CELL)

low = minimum(data_dict[CELL]["cycles"])
mid = median(data_dict[CELL]["cycles"])
if mid%1 != 0
    mid = floor(mid)
end
high = maximum(data_dict[CELL]["cycles"])

cycles = [low, mid, high]
range  = high - low
min_cell = -range
max_cell = range*1.25
norm = PythonPlot.matplotlib.colors.Normalize(vmin = min_cell, vmax=max_cell)
cmap = PythonPlot.matplotlib.cm.get_cmap("Blues")
mycolors = [cmap(norm(low)), cmap(norm(mid)), cmap(norm(high))]


x_type = "absolute"
y_type = "delta_N"

cells = [k for k in keys(data_dict)]


fig, ax = subplots(nrows=nrow,ncols=ncol,figsize=(15,9));
fig.subplots_adjust(hspace=0.05, wspace=0.05)

for (j,x) in enumerate(xs)
    if x == :frac_sol_am_neg
        xlabeltext = "fₛ⁻"
    elseif x == :frac_sol_am_pos
        xlabeltext = "fₛ⁺"
    elseif x == :n_li
        xlabeltext = "nₗᵢ"
    else
        xlabeltext = String(x)
    end
    for (i,y) in enumerate(ys)
        if y == :frac_sol_am_neg
            ylabeltext = "fₛ⁻"
        elseif y == :frac_sol_am_pos
            ylabeltext = "fₛ⁺"
        elseif y == :n_li
            ylabeltext = "nₗᵢ"
        else
            ylabeltext = String(y)
        end
        if i == j
            for (k,CYCLE) in enumerate(cycles)
                y_axis = data_dict[CELL]["distributions"][CYCLE][y].data
                x_axis = data_dict[CELL]["distributions"][CYCLE][x].data
                U = kde(x_axis)
                ax[i-1,j-1].fill_between(x=U.x,y1=U.density,alpha=0.5, color=cmap(norm(CYCLE)))
                ax[i-1,j-1].plot(U.x,U.density, color=cmap(norm(CYCLE)))
            end
        elseif i>j
            cs = zeros(3)
            for (k,CYCLE) in enumerate(cycles)
                y_axis = data_dict[CELL]["distributions"][CYCLE][y].data
                x_axis = data_dict[CELL]["distributions"][CYCLE][x].data
                cs[k] = round(cor(x_axis, y_axis),digits=3)
            end
            textArr = [pylist(string.([cycles[k],cs[k]])) for k in 1:length(cycles)]
            textColour = pylist([pylist(repeat([mycolors[k]],2)) for k in 1:length(cycles)])
            textArr = pylist(textArr)
            colLabels = pylist(["Cycle","ρ"])
            mytab = ax[i-1,j-1].table(cellText=textArr,cellColours=textColour, loc="center",colLabels=colLabels)
            for cell in mytab._cells
                mytab._cells[cell].set_alpha(.5)
            end
            ax[i-1,j-1].set_xticklabels([])
            ax[i-1,j-1].set_yticklabels([])
            ax[i-1,j-1].xaxis.set_ticks_position("none")
            ax[i-1,j-1].yaxis.set_ticks_position("none")
            ax[i-1,j-1].spines["top"].set_visible(false)
            ax[i-1,j-1].spines["right"].set_visible(false)
            ax[i-1,j-1].spines["bottom"].set_visible(false)
            ax[i-1,j-1].spines["left"].set_visible(false)

        elseif i<j
            for (k,CYCLE) in enumerate(cycles)
                y_axis = data_dict[CELL]["distributions"][CYCLE][y].data
                x_axis = data_dict[CELL]["distributions"][CYCLE][x].data
                ax[i-1,j-1].scatter(x_axis,y_axis,s=0.5, c=CYCLE*ones(length(x_axis)),cmap="Blues", vmax=max_cell, vmin=min_cell)
            end
        end

        y_axis = Float64[]
        for CYCLE in cycles
            y_axis = vcat(y_axis, data_dict[CELL]["distributions"][CYCLE][y].data)
        end
        max_y = maximum(y_axis)
        min_y = minimum(y_axis)


        δy = (max_y - min_y)
        t1 = round(δy*0.2 + min_y,digits=3)
        t2 = round(δy*0.8 + min_y,digits=3)
        ax[i-1,j-1].set_yticks(pylist([t1,t2]), labels=[])
        ax[i-1,j-1].grid(axis="x")

        if i != j
            ax[i-1,j-1].grid(axis="y")
            x_axis = Float64[]
            for CYCLE in cycles
                x_axis = vcat(x_axis, data_dict[CELL]["distributions"][CYCLE][x].data)
            end
            max_x = maximum(x_axis)
            min_x = minimum(x_axis)
    
            δx = (max_x - min_x)
            t1 = round(δx*0.2 + min_x,digits=3)
            t2 = round(δx*0.8 + min_x,digits=3)
            ax[i-1,j-1].set_xticks(pylist([t1,t2]))
            ax[i-1,j-1].xaxis.set_ticks_position("top")
            ax[i-1,j-1].xaxis.set_visible(true)
        end

        if i == j
            if i > 1
                lims = ax[i-2,j-1].get_xlim()
                ax[i-1,j-1].set_xlim(lims)
            end
            x_axis = Float64[]
            for CYCLE in cycles
                x_axis = vcat(x_axis, data_dict[CELL]["distributions"][CYCLE][x].data)
            end
            max_x = maximum(x_axis)
            min_x = minimum(x_axis)
    
            δx = (max_x - min_x)
            t1 = round(δx*0.2 + min_x,digits=3)
            t2 = round(δx*0.8 + min_x,digits=3)
            ax[i-1,j-1].set_xticks(pylist([t1,t2]))
            ax[i-1,j-1].xaxis.set_ticks_position("top")
        end



        ax[i-1,j-1].xaxis.set_visible(true)
        ax[i-1,j-1].yaxis.set_visible(true)
        ax[i-1,j-1].xaxis.set_tick_params(length=0)
        ax[i-1,j-1].yaxis.set_tick_params(length=0)
        ax[i-1,j-1].set_yticklabels([])
        ax[i-1,j-1].set_xticklabels([])
        

        # Set up ticks only on one side for the "edge" subplots...
        if (j == 1)
            ax[i-1,j-1].set_ylabel(ylabeltext,fontdict=pydict(Dict("size"=>22)))
            ax[i-1,j-1].yaxis.set_visible(true)
            ax[i-1,j-1].set_yticklabels([])
            ax[i-1,j-1].yaxis.set_tick_params(length=0)
        end
        if (j == ncol) & (i != nrow)
            y_axis = Float64[]
            for CYCLE in cycles
                y_axis = vcat(y_axis, data_dict[CELL]["distributions"][CYCLE][y].data)
            end

            max_y = maximum(y_axis)
            min_y = minimum(y_axis)


            δy = (max_y - min_y)
            t1 = round(δy*0.2 + min_y,digits=3)
            t2 = round(δy*0.8 + min_y,digits=3)
            ax[i-1,j-1].set_yticks(pylist([t1,t2]), labels=pylist([t1,t2]))
            ax[i-1,j-1].yaxis.set_tick_params(labelsize=18)
                
            ax[i-1,j-1].yaxis.set_ticks_position("right")
            ax[i-1,j-1].yaxis.set_tick_params(labelsize=18)
            ax[i-1,j-1].yaxis.set_visible(true)
        end
        if i == 1
            x_axis = Float64[]
            for CYCLE in cycles
                x_axis = vcat(x_axis, data_dict[CELL]["distributions"][CYCLE][x].data)
            end
            max_x = maximum(x_axis)
            min_x = minimum(x_axis)

            δx = (max_x - min_x)
            t1 = round(δx*0.2 + min_x,digits=3)
            t2 = round(δx*0.8 + min_x,digits=3)
            ax[i-1,j-1].set_xticks(pylist([t1,t2]))
            ax[i-1,j-1].xaxis.set_ticks_position("top")
            ax[i-1,j-1].xaxis.set_visible(true)
            ax[i-1,j-1].set_xticklabels(pylist([t1,t2]), rotation=15, ha="left",fontdict=pydict(Dict("size"=>18)), rotation_mode="anchor")
            ax[i-1,j-1].yaxis.set_tick_params(length=0)
            ax[i-1,j-1].yaxis.set_tick_params(length=0)
            #ax[i-1,j-1].grid(axis="y")
            #ax[i-1,j-1].grid(axis="x")
        end
        if i == nrow
            ax[i-1,j-1].xaxis.set_visible(true)
            ax[i-1,j-1].set_xticklabels([])
            ax[i-1,j-1].set_xlabel(xlabeltext,fontdict=pydict(Dict("size"=>22)))
            ax[i-1,j-1].tick_params(bottom=false, axis="x")
            
        end
        if i > j
            ax[i-1,j-1].grid()
        end
    end
    #fig.tight_layout()
end

norm = PythonPlot.matplotlib.colors.Normalize(vmin=min_cell,vmax=max_cell)
cmap = PythonPlot.matplotlib.pyplot.get_cmap("Blues")
sm = PythonPlot.matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array(np.array([0.0]))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.925, 0.15, 0.025, 0.7])

cb = fig.colorbar(sm, cax=cbar_ax,ticks=[cycles[1],cycles[3]-10],boundaries=np.arange(0,cycles[3]))
cb.ax.set_yticklabels(["Early\nLife ", "Late\nLife"], fontsize=14)
fig.savefig("figs/identifiability/$(CELL)_identifiability.png",bbox_inches="tight")
fig.savefig("figs/identifiability/$(CELL)_identifiability.pdf",bbox_inches="tight")
end
