using PythonPlot #assume data dict is loaded

CELL = "VAH11"
CYCLE = 502
figure(1)
clf()
d_f = data_dict[CELL]["distributions"][CYCLE][:n_li].data
d_e = data_dict[CELL]["distributions"][CYCLE][:ω].data
scatter(d_f,d_e)