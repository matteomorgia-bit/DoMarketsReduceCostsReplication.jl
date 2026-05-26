using FixedEffectModels
using DoMarketsReduceCostsReplication

const BASE_VARS = [
    :ln_emp,
    :ioudrg,
    :mnc_post92,
    :mnc_post87,
    :anwage_util,
    :ln_mwhs,
    :fgddum,
]

df = prepare_input_data(load_frw_data("frw1extract_enf.dta"))

println("Table 3, column 1: labor input demand, GLS basic")
println("Observations: ", size(df, 1))
println("Plant-epochs: ", count_plant_epochs(df))

model = reg(
    df,
    @formula(ln_emp ~ ioudrg + mnc_post92 + mnc_post87 + anwage_util + ln_mwhs + fgddum + fe(yr_data) + fe(plant_num2)),
    save = :residuals,
)

rho = estimate_rho(residuals(model), df.diff)
println("Initial rho: ", rho)

for iter in 1:20
    add_prais_columns!(df, BASE_VARS, rho)

    global model = reg(
        df,
        @formula(prais_ln_emp ~ prais_ioudrg + prais_mnc_post92 + prais_mnc_post87 + prais_anwage_util + prais_ln_mwhs + prais_fgddum + fe(yr_data) + fe(plant_num2)),
        save = :residuals,
    )

    new_rho = estimate_rho(residuals(model), df.diff)
    println("Iteration ", iter, " rho: ", new_rho)

    if abs(new_rho - rho) < 0.005
        global rho = new_rho
        break
    end

    global rho = new_rho
end

println()
println("Final rho: ", rho)
println(model)