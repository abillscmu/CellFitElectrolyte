#assumes data dict is loaded
using PythonPlot

sym = :frac_sol_am_pos
cells = [k for k in keys(data_dict)]
first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][sym] for cell in cells]
last_val = [data_dict[cell]["mean"][maximum(data_dict[cell]["cycles"])][sym] for cell in cells]
cycle_life = [maximum(data_dict[cell]["cycles"]) for cell in cells]

key = "Max DOD"

x_axis = [ex_data_dict["summary"][cell][key] for cell in cells]
Δ = (last_val .- first_val)
figure(1)
clf()
scatter(cycle_life, Δ)
ylabel("Average $sym change per cycle [%]")
xlabel(key)
C = cor(cycle_life, Δ)
title("Correlation: $C")
grid()

