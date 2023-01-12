using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity, DataFrames, PythonCall
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs1214_nuts_vah01/"

figure(1)
clf()


sym = [:εₑ⁺, :εₑ⁻, :n_li, :frac_sol_am_neg, :frac_sol_am_pos, :ω]
cells = []
lower_bounds = Dict()
upper_bounds = Dict()
median = Dict()
cycle = Dict()

for s in sym
    lower_bounds[s] = Float64[]
    upper_bounds[s] = Float64[]
    median[s] = Float64[]
    cycle[s] = Int[]
end


for file in readdir(FOLDERNAME)
    filename = FOLDERNAME * file
    arr = parse(Int, split(split(file, "VAH")[2],"_")[1])
    if arr == 2
        continue
    end
    cell = parse(Int,split(file,"_")[2])
    chain = try
        d = load(filename)
        chain = d["chain"]
    catch
        @warn "problem with $file"
        continue
    end
    quantiles = quantile(chain, q=[0.025, 0.5, 0.975])
    for s in sym
        push!(lower_bounds[s], quantiles[s].nt[:var"2.5%"][1])
        push!(upper_bounds[s], quantiles[s].nt[:var"97.5%"][1])
        push!(median[s], quantiles[s].nt[:var"50.0%"][1])
        push!(cycle[s], cell)
    end
end

df_dict = Dict()

for s in sym
    df = DataFrame(cycle=cycle[s], lower_bounds=lower_bounds[s], upper_bounds=upper_bounds[s], median=median[s])
    sorted_df = sort(df, :cycle)
    df_dict[s] = sorted_df
end

function my_plot(s)
    figure(1)
    clf()
    fill_between(np.array(pyfloat.(df_dict[s].cycle),dtype=np.float),np.array(pyfloat.(df_dict[s].lower_bounds),dtype=np.float),np.array(pyfloat.(df_dict[s].upper_bounds),dtype=np.float), facecolor="grey", alpha = 1)
    plot(df_dict[s].cycle, df_dict[s].median,"b")
    legend(["95%","Median"])
    xlabel("Cycle")
    ylabel(string(s))
    return nothing
end

#=
figure(1);
#clf()
fill_between(np.array(pyfloat.(sorted_df.cells),dtype=np.float),np.array(pyfloat.(sorted_df.lower_bound),dtype=np.float),np.array(pyfloat.(sorted_df.upper_bound),dtype=np.float), facecolor="grey", alpha = 1)
plot(sorted_df.cells, sorted_df.mid,"b")
#plot(sorted_df.cells, sorted_df.mle, "g")
legend(["95%","Median"])
xlabel("Cycle")
ylabel(string(sym))
PythonPlot.matplotlib.rcParams["font.size"] = 12
grid()
=#