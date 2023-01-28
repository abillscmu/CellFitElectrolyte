using CellFitElectrolyte, JLD2, PythonPlot, Turing, KernelDensity, PythonCall, CellFitElectrolyte.OCV
np = pyimport("numpy")
pygui(true)

FOLDERNAME = "results/outputs0117_elec/"
CELLS = ["VAH01","VAH02","VAH05"]
SYMBOLS = [:ω, :εₑ⁻, :εₑ⁺, :frac_sol_am_pos, :frac_sol_am_neg]



figure(1)
clf()
cellgeometry = CellFitElectrolyte.cell_geometry()
cathodeocv,anodeocv = CellFitElectrolyte.initialize_airbus_ocv()
sym = :n_li
num_rows = length(CELLS)
num_cols = length(SYMBOLS)
  
#for (i,CELL) in enumerate(CELLS)
    #for (j, SYMBOL) in enumerate(SYMBOLS)
    begin
        CELL = "VAH05"
        val = []
        cyc = []
        for cell in 2:2500
        cell_to_load = "$(FOLDERNAME)$(CELL)_$(cell)_HMC.jld2"
        chain = try
            d = load(cell_to_load)
            chain = d["chain"]
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

            append!(val, n_li)
            append!(cyc, cell*ones(1000))
        catch
            @warn "problem with $cell"
            continue
        end
    end

    hist2D(cyc, val, bins=100)
    #colorbar(label="Density")
    grid(alpha=0.8)
    ylabel("n_li")
    xlabel("cycle number")
    #title("Series: $CELL")
end

