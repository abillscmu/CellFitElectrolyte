using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0106_fullcyc/"
CELL = "VAH01"
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg]
minmax = Dict(sym => Dict("min" => Inf, "max" => -Inf) for sym in SYMBOLS)
CELLS = ["VAH01", "VAH01", "VAH01"]
CYCLES = [2, 420, 840]

val = []
cyc = []

figure(1)
pygui(true)
subplot(3, 5, 1)
clf()
for (i,cycleasdf) in enumerate(CYCLES)
    for (j,SYMBOL) in enumerate(SYMBOLS)
        subplot(3, 5, 5*(i-1) + j)
        cycle = 2
        cell_to_load = "$(FOLDERNAME)$(CELL)_$(cycleasdf)_HMC.jld2"
        println(cycle)
        chain = load(cell_to_load)["chain"]
        data = chain[SYMBOL].data[:,1]
        if minimum(data) < minmax[SYMBOL]["min"]
            minmax[SYMBOL]["min"] = minimum(data)
        end
        if maximum(data) > minmax[SYMBOL]["max"]
            minmax[SYMBOL]["max"] = maximum(data)
        end
        U = kde(data)
        plot(U.x, U.density);
        hist(data)
        grid()
        if i == 3
            xlabel(string(SYMBOL))
        end
        if j == 1
            ylabel("Density")
        end
    end
end

for i in 1:length(CYCLES)
    for (j,sym) in enumerate(SYMBOLS)
        subplot(3,5,5*(i-1) + j)
        xlim([minmax[sym]["min"], minmax[sym]["max"]])
    end
end

