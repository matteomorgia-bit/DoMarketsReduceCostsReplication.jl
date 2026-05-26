module DoMarketsReduceCostsReplication

using DataFrames
using StatFiles

export raw_data_dir, load_frw_data, count_plant_epochs, prepare_input_data

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
"""
    prepare_input_data(df)

Reproduce the data-preparation steps at the start of `inputregs.do`.
The function sorts by plant epoch and year, creates `diff`, and creates
`iounug`.
"""
function prepare_input_data(df)
    sorted = sort(df, [:plant_num2, :yr_data])
    out = copy(sorted)

    n = nrow(out)
    out.diff = falses(n)

    for i in 2:n
        same_plant = out.plant_num2[i] == out.plant_num2[i - 1]
        consecutive_year = out.yr_data[i] - out.yr_data[i - 1] == 1
        out.diff[i] = same_plant && consecutive_year
    end

    out.iounug = out.ioudum .* out.hi_nug
    out.iounug = ifelse.(out.yr_data .< 1993, 0, out.iounug)

    return out
end
end