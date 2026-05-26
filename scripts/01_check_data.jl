using DoMarketsReduceCostsReplication

function print_dataset_summary(filename)
    df = load_frw_data(filename)

    println("Dataset: ", filename)
    println("Observations: ", nrow(df))
    println("Plant-epochs: ", count_plant_epochs(df))
    println("Years: ", minimum(df.yr_data), " to ", maximum(df.yr_data))
    println()
end

print_dataset_summary("frw1extract_enf.dta")
print_dataset_summary("frw1extract_f.dta")