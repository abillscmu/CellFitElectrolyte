using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity
pygui(true)

FOLDERNAME = "results/outputs1212_vah02_hmc/"

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
    data_LI = chain[:n_li].data[:,1]
    data_ω = chain[:ω].data[:,1]
    data_εₛ⁻ = chain[:εₛ⁻].data[:,1]
    data_εₛ⁺ = chain[:εₛ⁺].data[:,1]
    data_εᵧ⁺ = chain[:εᵧ⁺].data[:,1]
    data_εᵧ⁻ = chain[:εᵧ⁻].data[:,1]
    figure(1)
    subplot(231)
    scatter(cell, mean(data_LI))
    subplot(232)
    scatter(cell, mean(data_ω))
    subplot(233)
    scatter(cell, mean(data_εₛ⁻))
    subplot(234)
    scatter(cell, mean(data_εₛ⁺))
    subplot(235)
    scatter(cell, mean(data_εᵧ⁻))
    subplot(236)
    scatter(cell, mean(data_εᵧ⁺))
end

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
