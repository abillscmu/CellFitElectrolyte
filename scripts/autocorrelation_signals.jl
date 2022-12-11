using StatsBase

dfi = 6

sym = :δ⁻


cell = unique(gdf[dfi].cell)

df_sorted = sort(gdf[dfi],order(:cycle))

ac = autocor(df_sorted[!,sym])

upper_99_bound = 2.576/sqrt(nrow(df_sorted))
lower_99_bound = -2.576/sqrt(nrow(df_sorted))

mk_results = mk.original_test(df_sorted[!,sym])
println("$(mk_results.trend)")

figure(1)
clf()
PythonPlot.matplotlib.rcParams["font.size"] = 18
plot(ac, marker="o")
plot([-1, length(ac)],[upper_99_bound, upper_99_bound],"k-")
plot([-1, length(ac)],[lower_99_bound, lower_99_bound],"k-")
grid()
xlim([0, length(ac)])
xlabel("Lag")
ylabel("ACF")
title("Series : VAH$(lpad(Int(cell[1]), 2, "0"))")