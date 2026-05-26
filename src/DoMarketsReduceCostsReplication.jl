module DoMarketsReduceCostsReplication

using DataFrames
using StatFiles

export raw_data_dir, load_frw_data, count_plant_epochs

"""
    raw_data_dir()

Return the path to the openICPSR raw data folder used by the replication scripts.
"""
function raw_data_dir()
    return joinpath(@__DIR__, "..", "data", "raw", "openicpsr", "116286-V1")
end

"""
    load_frw_data(filename)

Load one of the Stata datasets from the openICPSR replication package.
"""
function load_frw_data(filename)
    path = joinpath(raw_data_dir(), filename)
    return DataFrame(load(path))
end

"""
    count_plant_epochs(df)

Count unique plant-epoch identifiers.
"""
function count_plant_epochs(df)
    return length(unique(df.plant_num2))
end

end