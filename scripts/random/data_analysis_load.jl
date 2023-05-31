#assumes data dict is loaded
using PythonPlot

sym = :εₛ⁺
cells = [k for k in keys(data_dict)]
first_val = [data_dict[cell]["mean"][minimum(data_dict[cell]["cycles"])][sym] for cell in cells]
last_val = [data_dict[cell]["mean"][maximum(data_dict[cell]["cycles"])][sym] for cell in cells]
cycle_life = [maximum(data_dict[cell]["cycles"]) for cell in cells]

key = "Max Current"

x_axis = [ex_data_dict["summary"][cell][key] for cell in cells]
Δ = (last_val .- first_val)./cycle_life
figure(1)
clf()
PythonPlot.scatter(x_axis, Δ)
ylabel("Average $sym change per cycle [%]")
xlabel(key)
C = cor(x_axis, Δ)
title("Correlation: $C")
PythonPlot.grid()

