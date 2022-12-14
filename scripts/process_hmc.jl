using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity, DataFrames, PythonCall
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs1212_vah02_hmc/"

figure(1)
clf()

cells = []
n_li_lower_bound = []
n_li_upper_bound = []
n_li_50 = []
mle_vec = []

sym = :ω

for file in readdir(FOLDERNAME)
    filename = FOLDERNAME * file
    cell = parse(Int,split(file,"_")[2])
    chain = try
        d = load(filename)
        chain = d["chain"]
    catch
        @warn "problem with $file"
        continue
    end
    quantiles = quantile(chain, q=[0.025, 0.5, 0.975])
    low_n_li = quantiles[sym].nt[:var"2.5%"][1]
    high_n_li = quantiles[sym].nt[:var"97.5%"][1]
    mid_n_li = quantiles[sym].nt[:var"50.0%"][1]
    U = kde(chain[sym].data[:,1])
    mle = U.x[argmax(U.density)]
    push!(cells, cell)
    push!(n_li_lower_bound, low_n_li)
    push!(n_li_upper_bound, high_n_li)
    push!(n_li_50, mid_n_li)
    push!(mle_vec, mle)
end

df = DataFrame(cells=cells, lower_bound=n_li_lower_bound, upper_bound=n_li_upper_bound, mid=n_li_50, mle=mle_vec)

sorted_df = sort(df, :cells)

figure(1);
clf()
fill_between(np.array(pyfloat.(sorted_df.cells),dtype=np.float),np.array(pyfloat.(sorted_df.lower_bound),dtype=np.float),np.array(pyfloat.(sorted_df.upper_bound),dtype=np.float), facecolor="grey")
plot(sorted_df.cells, sorted_df.mid,"b")
plot(sorted_df.cells, sorted_df.mle, "g")
legend(["95%","50% Quantile", "MAP"])
xlabel("Cycle")
ylabel(string(sym))
PythonPlot.matplotlib.rcParams["font.size"] = 12
grid()