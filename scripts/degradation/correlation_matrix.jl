using PythonPlot, Statistics, PythonCall


ys = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg, :n_li, :εₛ⁻, :εₛ⁺]
xs = ["Max Temperature", "Max DOD", "Max Current", "Cycle Life","Mean Temperature", "Min Current"]


nrow = length(ys)
ncol = length(xs)

xlabeltextarr = Array{String}(undef, ncol)
ylabeltextarr = Array{String}(undef, nrow)


x_type = "absolute"
y_type = "delta_N"

cells = [k for k in keys(data_dict)]

Z = zeros(nrow, ncol)

for (j,x) in enumerate(xs)
    if x in keys(data_dict["VAH01"]["distributions"][2])
        first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][x] for cell in cells]
        last_val = [data_dict[cell]["mean"][maximum(data_dict[cell]["cycles"])][x] for cell in cells]
        if x == :frac_sol_am_neg
            x_sym = "fₛ⁻"
        elseif x == :frac_sol_am_pos
            x_sym = "fₛ⁺"
        elseif x == :n_li
            x_sym = "nₗᵢ"
        else
            x_sym = string(x)
        end
        if x_type == "first"
            x_axis = first_val
            xlabeltext = "Initial $x_sym"
        elseif x_type == "last"
            x_axis = last_val
            xlabeltext = "Final $x_sym"
        elseif x_type == "delta"
            x_axis = last_val .- first_val
            xlabeltext = "Δ$x_sym"
        elseif x_type == "delta_N"
            cycle_life = [maximum(data_dict[cell]["cycles"]) for cell in cells]
            x_axis = (last_val .- first_val)./cycle_life
            xlabeltext = "Δ$x_sym/Ν"
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
    for (i,y) in enumerate(ys)
        if y in keys(data_dict["VAH01"]["distributions"][2])
            first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][y] for cell in cells]
            last_val = [data_dict[cell]["mean"][maximum(data_dict[cell]["cycles"])][y] for cell in cells]
            if y == :frac_sol_am_neg
                y_sym = "fₛ⁻"
            elseif y == :frac_sol_am_pos
                y_sym = "fₛ⁺"
            elseif y == :n_li
                y_sym = "nₗᵢ"
            else
                y_sym = string(y)
            end
            if y_type == "first"
                y_axis = first_val
                ylabeltext = "Initial $y_sym"
            elseif y_type == "last"
                y_axis = last_val
                ylabeltext = "Final $y_sym"
            elseif y_type == "delta"
                y_axis = last_val .- first_val
                ylabeltext = "Δ$y_sym"
            elseif y_type == "delta_N"
                cycle_life = [maximum(data_dict[cell]["cycles"]) for cell in cells]
                y_axis = (last_val .- first_val)./cycle_life
                ylabeltext = "Δ$y_sym/N"
            else
                error("y not found")
            end
        elseif y in keys(ex_data_dict["summary"]["VAH01"])
            y_axis = [ex_data_dict["summary"][cell][y] for cell in cells]
            ylabeltext = y
        elseif y == "Cycle Life"
            y_axis = [maximum(data_dict[cell]["cycles"]) for cell in cells]
            ylabeltext = y
        else
            error("y not found")
        end
        ylabeltextarr[i] = ylabeltext
        xlabeltextarr[j] = xlabeltext
        C = round(cor(x_axis, y_axis),digits=3)
        Z[i,j] = C
    end
    #fig.tight_layout()
end

fig, ax = subplots(figsize=(8,6.4))
cmap = PythonPlot.matplotlib.colormaps["PRGn"]
pcm = ax.imshow(Z,aspect="auto", cmap=cmap)
colorbar(pcm, ax=ax)
#=
for i in 1:nrow
    for j in 1:ncol
        text(j-1, i-1, Z[i, j], ha="center", va="center", color="w",fontsize=13)
    end
end
=#
ax.set_xticks(0:length(xs)-1,labels=xlabeltextarr,rotation=45,ha="right",fontsize=15)
ax.set_yticks(0:length(ys)-1,labels=ylabeltextarr,fontsize=15)
fig.align_ylabels()
fig.tight_layout()
fig.savefig("figs/correlations/correlation_heatmap_$(y_type).pdf",bbox_inches="tight")
fig.savefig("figs/correlations/correlation_heatmap_$(y_type).png",bbox_inches="tight")
