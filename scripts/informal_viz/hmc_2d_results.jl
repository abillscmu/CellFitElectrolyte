using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity
pygui(true)

FOLDERNAME = "results/newnuts_5/"

figure(1)
clf()



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
    #data_LI = chain[:n_li].data[:,1]
    data_ω = chain[:ω].data[:,1]
    data_εₛ⁻ = chain[:εₑ⁻].data[:,1]
    data_εₛ⁺ = chain[:εₑ⁺].data[:,1]
    data_εᵧ⁺ = chain[:frac_sol_am_neg].data[:,1]
    data_εᵧ⁻ = chain[:frac_sol_am_pos].data[:,1]
    figure(1)
    subplot(232)
    scatter(cell, std(data_ω))
    subplot(233)
    scatter(cell, std(data_εₛ⁻))
    subplot(234)
    scatter(cell, std(data_εₛ⁺))
    subplot(235)
    scatter(cell, std(data_εᵧ⁻))
    subplot(236)
    scatter(cell, std(data_εᵧ⁺))
end
#=
figure(1)
subplot(231)
xlim([0,650])
subplot(232)
xlim([0,650])
subplot(233)
xlim([0,650])
subplot(234)
xlim([0,650])
subplot(235)
xlim([0,650])
subplot(236)
xlim([0,650])
=#
