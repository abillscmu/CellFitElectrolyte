using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0106_fullcyc/"
CELL = "VAH01"
SYMBOL = :ω

val = []
cyc = []

figure(1)
clf()



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

figure(1)
clf()
hist2D(cyc, val, bins=100)
colorbar(label="Density")
grid(alpha=0.8)
xlabel("Cycle Number")
ylabel(string(SYMBOL))
title("Series: $CELL")

