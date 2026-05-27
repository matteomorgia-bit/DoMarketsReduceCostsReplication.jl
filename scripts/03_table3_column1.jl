using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra

const YEAR_VARS = Symbol.("year", 82:99)

function dummy_matrix(values)
    levels = sort(unique(values))
    return hcat([Float64.(values .== level) for level in levels]...)
end

function build_x(df, vars)
    main = Matrix{Float64}(df[:, vars])
    years = Matrix{Float64}(df[:, YEAR_VARS])
    plants = dummy_matrix(df.plant_num2)
    constant = ones(size(df, 1), 1)

    return hcat(main, years, plants, constant)
end

function transform_matrix(X, diff, rho)
    out = similar(X, Float64)

    for j in axes(X, 2)
        out[:, j] = prais_transform(X[:, j], diff, rho)
    end

    return out
end

function ols_fit(y, X)
    beta = pinv(X) * y
    residuals = y - X * beta
    return beta, residuals
end

function prais_gls(df, depvar, vars; maxiter = 20, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    X = build_x(df, vars)

    beta, residuals_raw = ols_fit(y, X)
    rho = estimate_rho(residuals_raw, df.diff)

    for iter in 1:maxiter
        y_prais = prais_transform(y, df.diff, rho)
        X_prais = transform_matrix(X, df.diff, rho)

        beta, _ = ols_fit(y_prais, X_prais)
        residuals_original = y - X * beta
        new_rho = estimate_rho(residuals_original, df.diff)

        if abs(new_rho - rho) < tolerance
            rho = new_rho
            break
        end

        rho = new_rho
    end

    return (; beta, rho)
end

function print_result(title, vars, result)
    println()
    println(title)
    println("rho: ", result.rho)
    for (name, coef) in zip(String.(vars), result.beta[1:length(vars)])
        println(rpad(name, 14), coef)
    end
end

vars_emp = [:ioudrg, :mnc_post92, :mnc_post87, :anwage_util, :ln_mwhs, :fgddum]
vars_no_wage = [:ioudrg, :mnc_post92, :mnc_post87, :ln_mwhs, :fgddum]

enf = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
fuel = prepare_input_data(load_frw_data("frw1extract_f.dta"))

println("GLS basic columns for Tables 3, 4, and 5")
println("ENF observations: ", size(enf, 1), ", plant-epochs: ", count_plant_epochs(enf))
println("Fuel observations: ", size(fuel, 1), ", plant-epochs: ", count_plant_epochs(fuel))

result_t3 = prais_gls(enf, :ln_emp, vars_emp)
result_t4 = prais_gls(enf, :ln_nfexp, vars_no_wage)
result_t5 = prais_gls(fuel, :ln_btu, vars_no_wage)

print_result("Table 3 column 1: ln_emp", vars_emp, result_t3)
print_result("Table 4 column 1: ln_nfexp", vars_no_wage, result_t4)
print_result("Table 5 column 1: ln_btu", vars_no_wage, result_t5)