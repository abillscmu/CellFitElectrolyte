FOLDERNAME = "results/outputs1208_vah02/"



using DataFrames,CSV

arr = readdir(FOLDERNAME)
df = CSV.read("$FOLDERNAME/$(arr[1])",DataFrame)
for a in arr[2:end]
    ndf = CSV.read("$FOLDERNAME/$(a)",DataFrame)
    global df = vcat(df,ndf)
end

CSV.write("$FOLDERNAME.csv",df)
