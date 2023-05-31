t = []
for cell in keys(data_dict)
    for cycle_num in keys(data_dict[cell]["time"])
        append!(t, data_dict[cell]["time"][cycle_num])
    end
end
t2 = sort(t)
i = 1:length(t2)
mt2 = [maximum(t2[1:n]) for n in i]
using PythonPlot
figure(1)
clf()
PythonPlot.step(collect(i), mt2)
PythonPlot.xlabel("Number of Cycles (Sorted by Walltime)")
PythonPlot.ylabel("Time to Analyze [h]")

figure(2)
clf()
PythonPlot.step(collect(i), cumsum(mt2))
PythonPlot.xlabel("Number of Cycles (Sorted by Walltime)")
PythonPlot.ylabel("Time to Analyze [h]")