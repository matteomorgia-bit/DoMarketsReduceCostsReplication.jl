using DoMarketsReduceCostsReplication

const TABLE34_VARS = [
    :ln_emp,
    :ln_nfexp,
    :ioudrg,
    :ioulaw,
    :iouretail,
    :mnc_post92,
    :mnc_post87,
    :anwage_util,
    :ln_mwhs,
    :fgddum,
    :lstsales,
    :stateyr,
    :plant_num2,
    :yr_data,
    :ioudum,
    :hi_nug,
    :iounug,
    :diff,
]

const TABLE5_VARS = [
    :ln_btu,
    :ioudrg,
    :ioulaw,
    :iouretail,
    :mnc_post92,
    :mnc_post87,
    :ln_mwhs,
    :fgddum,
    :lstsales,
    :stateyr,
    :plant_num2,
    :yr_data,
    :ioudum,
    :hi_nug,
    :iounug,
    :diff,
]

function print_variable_check(df, vars)
    for var in vars
        missing_count = count(ismissing, df[!, var])
        println(rpad(String(var), 14), " missing: ", missing_count)
    end
end

println("Tables 3 and 4 data")
enf = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
println("Observations: ", size(enf, 1))
println("Plant-epochs: ", count_plant_epochs(enf))
print_variable_check(enf, TABLE34_VARS)

println()
println("Table 5 data")
fuel = prepare_input_data(load_frw_data("frw1extract_f.dta"))
println("Observations: ", size(fuel, 1))
println("Plant-epochs: ", count_plant_epochs(fuel))
print_variable_check(fuel, TABLE5_VARS)