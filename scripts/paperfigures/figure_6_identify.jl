using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall, KDEDistributions
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0106_fullcyc/"
CELL = "VAH01"
SYMBOL = :frac_sol_am_pos
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_neg, :frac_sol_am_pos]

cyc = []

thing = Dict()



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
        append!(cyc, cell)
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

for (i,sym) in enumerate(SYMBOLS)
    println(sym)
    totalmin = minimum(thing[sym]["min"])
    totalmax = maximum(thing[sym]["max"])
    totalrange = totalmax-totalmin
    normalizedvar = thing[sym]["val"] ./ totalrange
    arr[i,:] .= normalizedvar
end

figure(1)
clf()
pcolor(arr)
xlabel("Cycle")
colorbar()
yticks(0.5:4.5, SYMBOLS)


