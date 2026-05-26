using FixedEffectModels
using DoMarketsReduceCostsReplication

df = prepare_input_data(load_frw_data("frw1extract_enf.dta"))

println("Table 3, column 1: labor input demand, GLS basic")
println("Observations: ", size(df, 1))
println("Plant-epochs: ", count_plant_epochs(df))

ols_model = reg(
    df,
    @formula(ln_emp ~ ioudrg + mnc_post92 + mnc_post87 + anwage_util + ln_mwhs + fgddum + fe(yr_data) + fe(plant_num2)),
    save = :residuals,
)

residuals_ols = residuals(ols_model)
rho = estimate_rho(residuals_ols, df.diff)

println("Initial rho estimate: ", rho)