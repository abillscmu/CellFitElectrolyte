using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0106_fullcyc/"
CELLS = ["VAH27", "VAH28", "VAH30"]
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg]



figure(1)
clf()

num_rows = length(CELLS)
num_cols = length(SYMBOLS)

for (i,CELL) in enumerate(CELLS)
    for (j, SYMBOL) in enumerate(SYMBOLS)
        val = []
        cyc = []
        for cell in 0:2500
        cell_to_load = "$(FOLDERNAME)$(CELL)_$(cell)_HMC.jld2"
        chain = try
            d = load(cell_to_load)
            chain = d["chain"]
            append!(val, chain[SYMBOL].data[:,1])
            append!(cyc, cell*ones(1000))
        catch
            @warn "problem with $cell"
            continue
        end
    end

    subplot(num_rows, num_cols, num_cols*(i-1) + j)
    hist2D(cyc, val, bins=100)
    #colorbar(label="Density")
    grid(alpha=0.8)
    ylabel(string(SYMBOL))
    if i == 3
        xlabel("Cycle Number")
    end
    #title("Series: $CELL")
end
end 

