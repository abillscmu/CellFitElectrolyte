using Turing, DynamicHMC, PythonPlot, JLD2, KernelDensity
pygui(true)

FOLDERNAME = "results/outputs1212_vah02_hmc/"

VAH = "VAH02_11"
clear = false

["57","250","500","622","11"]


d = JLD2.load(FOLDERNAME*VAH*"_HMC.jld2")


chain = d["chain"]


U_LI = kde(chain[:n_li].data[:,1])
U_ω = kde(chain[:ω].data[:,1])
U_εₛ⁻ = kde(chain[:εₛ⁻].data[:,1])
U_εₛ⁺ = kde(chain[:εₛ⁺].data[:,1])
U_εᵧ⁺ = kde(chain[:εᵧ⁺].data[:,1])
U_εᵧ⁻ = kde(chain[:εᵧ⁻].data[:,1])

figure(1)
if clear
    clf()
end
subplot(231)
plot(U_LI.x,U_LI.density)
xlabel("n_li [mol]")
ylabel("Density")
grid()
subplot(232)
plot(U_ω.x,U_ω.density)
xlabel("ω [Ω]")
ylabel("Density")
grid()
subplot(233)
plot(U_εₛ⁻.x,U_εₛ⁻.density)
xlabel("εₛ⁻")
ylabel("Density")
grid()
subplot(234)
plot(U_εₛ⁺.x,U_εₛ⁺.density)
xlabel("εₛ⁺")
ylabel("Density")
grid()
subplot(235)
plot(U_εᵧ⁻.x,U_εᵧ⁻.density)
xlabel("εᵧ⁻")
ylabel("Density")
grid()
subplot(236)
plot(U_εᵧ⁺.x,U_εᵧ⁺.density)
xlabel("εᵧ⁺")
ylabel("Density")
grid()


