df = CSV.read("/Users/abills/Datasets/cycle_individual_data/VAH01_20.csv",DataFrame)
df.times = df.times .- df.times[1]
fig,axes = subplots(2,1,sharex=true)
axes[0].plot(df.times,df.EcellV)
axes[1].plot(df.times,df.TemperatureC)


df = CSV.read("/Users/abills/Datasets/cycle_individual_data/VAH01_400.csv",DataFrame)
df.times = df.times .- df.times[1]
axes[0].plot(df.times,df.EcellV)
axes[1].plot(df.times,df.TemperatureC)

df = CSV.read("/Users/abills/Datasets/cycle_individual_data/VAH01_845.csv",DataFrame)
df.times = df.times .- df.times[1]
axes[0].plot(df.times,df.EcellV)
axes[1].plot(df.times,df.TemperatureC)


axes[0].legend(["VAH01 Cycle 20","VAH01 Cycle 400","VAH01 Cycle 845"])
axes[1].set_xlabel("Time [s]")
axes[0].set_ylabel("Voltage [V]")
axes[1].set_ylabel("Temperature [\$^{\\mathrm{o}}\$C]")
fig.set_figwidth(8)
fig.savefig("figs/data_normalcycle.pdf",bbox_inches="tight")
fig.savefig("figs/data_normalcycle.png",bbox_inches="tight")