using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0112_elec/"
CELL = "VAH01"
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg, :c_e_0]
minmax = Dict(sym => Dict("min" => Inf, "max" => -Inf) for sym in SYMBOLS)
CELLS = ["VAH01"]
CYCLES = [200, 300, 400, 500, 600, 700, 800]

val = []
cyc = []

figure(1)
pygui(true)
subplot(length(CYCLES), length(SYMBOLS), 1)
clf()
for (i,cycleasdf) in enumerate(CYCLES)
    for (j,SYMBOL) in enumerate(SYMBOLS)
        subplot(length(CYCLES), length(SYMBOLS), length(SYMBOLS)*(i-1) + j)
        cell_to_load = "$(FOLDERNAME)$(CELL)_$(cycleasdf)_HMC.jld2"
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
        if i == length(CYCLES)
            xlabel(string(SYMBOL))
        end
        if j == 1
            ylabel("Density")
        end
    end
end

for i in 1:length(CYCLES)
    for (j,sym) in enumerate(SYMBOLS)
        subplot(length(CYCLES),length(SYMBOLS),length(SYMBOLS)*(i-1) + j)
        xlim([minmax[sym]["min"], minmax[sym]["max"]])
    end
end

