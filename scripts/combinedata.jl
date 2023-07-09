using DataFrames,CSV,ProgressMeter
dir = ARGS[1]
name=ARGS[2]
arr = readdir(dir)
df = CSV.read("$(dir)/$(arr[1])",DataFrame)
@showprogress for i in 2:(length(arr)-1)
	ndf = CSV.read("$(dir)/$(arr[i])",DataFrame)
	global df = vcat(df,ndf)
end

CSV.write(name,df)
