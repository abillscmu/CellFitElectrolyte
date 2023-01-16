using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity, DataFrames, PythonCall, CellFitElectrolyte, CellFitElectrolyte.OCV
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0114_elec/"

figure(1)
clf()

cells = []
n_li_lower_bound = []
n_li_upper_bound = []
n_li_50 = []
mle_vec = []

sym = :n_li
vah = "VAH01"
cell = for file in readdir(FOLDERNAME)
    if !occursin(vah, file)
        continue
    end
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

    εₑ⁻ = chain[:εₑ⁻]
    εₑ⁺ = chain[:εₑ⁺]
    frac_sol_am_pos = chain[:frac_sol_am_pos]
    frac_sol_am_neg = chain[:frac_sol_am_neg]
    c_e_0 = 4000

    V = 4.2
    x⁻₀ = 0.6
    V⁺₀ = V + calcocv(anodeocv,x⁻₀,297.0)
    x⁺₀ = OCV.get_x_from_voltage(cathodeocv,V⁺₀,297.0)
    c_n_init = x⁻₀*(anodeocv.c_s_max-anodeocv.c_s_min)
    c_p_init = x⁺₀*(cathodeocv.c_s_max-cathodeocv.c_s_min)

    #electrolyte
    n_li_elec = (εₑ⁻.*cellgeometry.Vₑ⁻ .+ εₑ⁺*cellgeometry.Vₑ⁺ .+ 0.8.*cellgeometry.Vₑˢ).*c_e_0
    n_li_pos = ((1 .- εₑ⁺) .* frac_sol_am_pos) .* cellgeometry.Vₑ⁺ .* c_p_init
    n_li_neg = ((1 .- εₑ⁻) .* frac_sol_am_neg) .* cellgeometry.Vₑ⁻  .* c_n_init
    n_li = n_li_pos .+ n_li_neg .+ n_li_elec



    
    U = kde(n_li[:,1])

    q = quantile(n_li[:,1])
    low_n_li = q[1]
    mid_n_li = q[2]
    high_n_li = q[3]

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
#clf()
fill_between(np.array(pyfloat.(sorted_df.cells),dtype=np.float),np.array(pyfloat.(sorted_df.lower_bound),dtype=np.float),np.array(pyfloat.(sorted_df.upper_bound),dtype=np.float), facecolor="grey", alpha = 1)
plot(sorted_df.cells, sorted_df.mid,"b")
plot(sorted_df.cells, sorted_df.mle, "g")
legend(["95%","Median"])
xlabel("Cycle")
ylabel(string(sym))
PythonPlot.matplotlib.rcParams["font.size"] = 12
grid()