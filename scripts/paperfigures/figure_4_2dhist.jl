using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall
np = pyimport("numpy")
#pygui(true)


FOLDERNAME = "results/outputs0117_elec/"
CELLS = ["VAH06","VAH07","VAH09"]
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg]



figure(1)
clf()
PythonPlot.matplotlib.rcParams["font.size"] = 16
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
    #grid(0.2)
    if j == 1
        ylabel(CELL*": "*string(SYMBOL))
    else
        ylabel(string(SYMBOL))
    end
    if i == length(CELLS)
        xlabel("Cycle Number")
    end
    #title("Series: $CELL")
end
end 
tight_layout()
savefig("figs/2dhist$(CELLS[1])$(CELLS[2])$(CELLS[3]).png")
savefig("figs/2dhist$(CELLS[1])$(CELLS[2])$(CELLS[3]).eps")