using PythonPlot

num_cells = length(keys(celldict))
cells = [k for k in keys(celldict)]
num_cycs = 3
num_mechanisms = length(deg_df.Name)
heatmap_matrix = zeros(num_mechanisms, num_cells, num_cycs)


for k in 1:num_cycs
    for (j,cell) in enumerate(cells)
        num_cycs_this_cell = size(celldict[cell])[1]
        cyc = Int(round(k*num_cycs_this_cell/num_cycs))
        for (i,mechanism) in enumerate(deg_df.Name)
            heatmap_matrix[i,j,k] = celldict[cell][cyc, i]
        end
    end
end

cmap = PythonPlot.matplotlib.colormaps["PRGn"]

fig,ax = subplots(nrows=3, figsize = (5, 9))
pos = ax[0].pcolormesh(heatmap_matrix[:,:,1],vmin=-1,vmax=1,cmap=cmap)
ax[0].set_yticks((1:num_mechanisms).-0.5, deg_df.Name, fontsize=14)
ax[0].get_xaxis().set_visible(false)
ax[1].pcolormesh(heatmap_matrix[:,:,2],vmin=-1,vmax=1,cmap=cmap)
ax[1].set_yticks((1:num_mechanisms).-0.5, deg_df.Name, fontsize=14)
ax[1].get_xaxis().set_visible(false)
ax[2].pcolormesh(heatmap_matrix[:,:,3],vmin=-1,vmax=1,cmap=cmap)
ax[2].set_yticks((1:num_mechanisms).-0.5, deg_df.Name, fontsize=14)
ax[2].set_xticks((1:length(cells)).-0.5,cells,rotation=90, ha="right",rotation_mode="anchor", fontsize=14)
cbar = fig.colorbar(pos, ax=ax, ticks = [-1 , 1], orientation="horizontal", location="top")
cbar.ax.set_xticklabels(["Least Similar" , "Most Similar"], fontsize=14)
fig.savefig("figs/heatmap_cycles.pdf",bbox_inches="tight")
fig.savefig("figs/heatmap_cycles.png",bbox_inches="tight")