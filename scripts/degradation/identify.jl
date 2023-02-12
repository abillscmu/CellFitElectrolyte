using PythonPlot
pygui(true)


cell = "VAH15"
sym = :frac_sol_am_neg

variances = [std(data_dict[cell]["distributions"][cycle][sym].data) for cycle in data_dict[cell]["cycles"]]
min_cyc = minimum(data_dict[cell]["cycles"])
max_cyc = maximum(data_dict[cell]["cycles"])
figure(1)
clf()
scatter(data_dict[cell]["cycles"],log.(variances))
std_var = std(Uniform(0.5,1.0))
plot([min_cyc,max_cyc],log.(std_var*ones(2)),"k--")
ylim([-10,0])