using CSV, DataFrames

df = CSV.read("data/outputs/outputs0802_full.csv",DataFrame)

p = CellFitElectrolyte.p_transport()

newp = deepcopy(p)


for param in keys(p)
    try
        newp[param] = mean(df[:,param])
    catch
        @warn "$param not found"
    end
end