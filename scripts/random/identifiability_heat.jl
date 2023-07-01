using PythonPlot, Statistics, PythonCall


SYMBOLS = SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg, :n_li, :εₛ⁻, :εₛ⁺]
xs = ys = SYMBOLS


nrow = length(ys)
ncol = length(xs)
CELL = "VAH01"
CYCLE = 2
Z = zeros(length(xs),length(ys))

xlabeltextarr = Array{String}(undef, ncol)
ylabeltextarr = Array{String}(undef, nrow)


cells = [k for k in keys(data_dict)]



for (j,x) in enumerate(xs)
    x_axis = data_dict[CELL]["distributions"][CYCLE][x].data
    xlabeltextarr[j] = String(x)
    for (i,y) in enumerate(ys)
        y_axis = data_dict[CELL]["distributions"][CYCLE][y].data
        C = round(cor(x_axis, y_axis),digits=3)
        Z[i,j] = C
        ylabeltextarr[i] = String(y)
    end
    #fig.tight_layout()
end

fig = figure(1)
clf()
imshow(Z, aspect="auto")

for i in 1:nrow
    for j in 1:ncol
        text(j-1, i-1, Z[i, j], ha="center", va="center", color="w")
    end
end
xticks(ticks = 0:length(xs)-1,labels=xlabeltextarr,rotation=45,ha="right")
yticks(ticks=0:length(ys)-1,labels=ylabeltextarr, rotation=45, ha="right")
fig.align_ylabels()
tight_layout()

