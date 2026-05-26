using DataFrames
using StatFiles

const RAW_DIR = joinpath("data", "raw", "openicpsr", "116286-V1")

function load_stata_file(filename)
    path = joinpath(RAW_DIR, filename)
    return DataFrame(load(path))
end

function count_plant_epochs(df)
    return length(unique(df.plant_num2))
end

function print_dataset_summary(filename)
    df = load_stata_file(filename)

    println("Dataset: ", filename)
    println("Observations: ", nrow(df))
    println("Plant-epochs: ", count_plant_epochs(df))
    println("Years: ", minimum(df.yr_data), " to ", maximum(df.yr_data))
    println()
end

print_dataset_summary("frw1extract_enf.dta")
print_dataset_summary("frw1extract_f.dta")