using CSV, DataFrames, Statistics, ProgressMeter

FOLDERNAME = "/Users/abills/Datasets/cycle_individual_data/"


ex_data_dict = Dict(
    "cells" => Dict(),
    "summary" => Dict()
)

@showprogress for file in readdir(FOLDERNAME)
    if file == ".DS_Store"
        continue
    end
    arr = split(file, ['_','.'])
    vah = arr[1]
    cycle = parse(Int,arr[2])
    df = CSV.read(FOLDERNAME*file, DataFrame)
    
    mean_temperature = mean(df.TemperatureC.+273)
    max_temperature = maximum(df.TemperatureC.+273)
    min_temperature = minimum(df.TemperatureC.+273)

    mean_current = mean(-df.ImA./1000)
    max_current = maximum(-df.ImA./1000)
    min_current = minimum(-df.ImA./1000)

    mean_voltage = mean(-df.EcellV./1000)
    max_voltage = maximum(df.EcellV)
    min_voltage = minimum(df.EcellV)

    mean_dod = mean(df.QDischargemAh)
    max_dod = maximum(df.QDischargemAh)
    min_dod = minimum(df.QDischargemAh)

    if !(vah in keys(ex_data_dict["cells"]))
        ex_data_dict["cells"][vah] = Dict(
            "cycles" => Int[],
            "data" => Dict()
        )
    end
    append!(ex_data_dict["cells"][vah]["cycles"], cycle)
    data = Dict(
        "Mean Temperature" => mean_temperature,
        "Max Temperature" => max_temperature,
        "Min Temperature" => min_temperature,
        "Mean Current" => mean_current,
        "Max Current" => max_current,
        "Min Current" => min_current,
        "Mean Voltage" => mean_voltage,
        "Max Voltage" => max_voltage,
        "Min Voltage" => min_voltage,
        "Mean DOD" => mean_dod,
        "Max DOD" => max_dod,
        "Min DOD" => min_dod
    )
    ex_data_dict["cells"][vah]["data"][cycle] = data
end

function get_mean(ex_data_dict, key, cell)
    return mean(ex_data_dict["cells"][cell]["data"][cycle][key] for cycle in ex_data_dict["cells"][cell]["cycles"])
end


for cell in keys(ex_data_dict["cells"])
    ex_data_dict["summary"][cell] = Dict()
    for key in keys(ex_data_dict["cells"]["VAH01"]["data"][1])
        val = get_mean(ex_data_dict, key, cell)
        ex_data_dict["summary"][cell][key] = val
    end
end

