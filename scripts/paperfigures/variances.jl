using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall, KDEDistributions
np = pyimport("numpy")
#pygui(true)

FOLDERNAME = "results/outputs0117_elec/"
files = readdir(FOLDERNAME)
cells = similar(files)
for (i,f) in enumerate(files)
    cells[i] = split(f, "_")[1]
end
cells = unique(cells)
files = 1
SYMBOL = :frac_sol_am_pos
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]
SYMBOL_LABELS = ["ω [Ω]","εₑ⁻","εₑ⁺","fₛ⁻","fₛ⁺"]


#cells = ["VAH01"]
for CELL in cells
    fig, ax = subplots(length(SYMBOLS),1,figsize=(12.5,8))
thing = Dict()


cyc = []
for sym in SYMBOLS
minval = []
maxval = []

val = []


for cell in 0:2500
    cell_to_load = "$(FOLDERNAME)$(CELL)_$(cell)_HMC.jld2"
    chain = try
        d = load(cell_to_load)
        chain = d["chain"]
        data = chain[sym].data[:,1]
        U = kde(chain[sym].data[:,1])
        fisher = std(data)
        append!(val, fisher)
        append!(minval, minimum(data))
        append!(maxval, maximum(data))
        if sym == :ω
            append!(cyc, cell)
        end
    catch
        @warn "problem with $cell"
        continue
    end
end
thing[sym] = Dict(
    "min" => minval,
    "max" => maxval,
    "val" => val
)
end

num_rows = length(SYMBOLS)
num_cols = length(thing[:ω]["val"])

arr = zeros(num_rows, num_cols)
cmap = PythonPlot.get_cmap("Blues_r")

for (i,sym) in enumerate(SYMBOLS)

    indices = sortperm(cyc)
    cyc_sorted = cyc[indices]
    val_sorted = thing[sym]["val"][indices]
    
    cell_num = parse(Int,CELL[end-1:end])/35
    ax[i-1].plot(cyc_sorted,val_sorted)
    ax[i-1].set_ylabel("σ($(SYMBOL_LABELS[i]))")
end
ax[-1].set_xlabel("Cycle Number")
fig.savefig("figs/si/variances/$CELL.png",bbox_inches="tight")
fig.savefig("figs/si/variances/$CELL.pdf",bbox_inches="tight")
end


