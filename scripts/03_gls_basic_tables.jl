using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra

const YEAR_VARS = Symbol.("year", 82:99)

function dummy_matrix(values)
    levels = sort(unique(values))
    levels = levels[2:end]
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
    beta = X \ y
    residuals = y - X * beta
    return beta, residuals
end

function cluster_vcov(X, residuals, cluster)
    bread = inv(X' * X)
    meat = zeros(size(X, 2), size(X, 2))

    for g in unique(cluster)
        idx = findall(cluster .== g)
        xg = X[idx, :]
        ug = residuals[idx]
        score = xg' * ug
        meat .+= score * score'
    end

    n = size(X, 1)
    k = size(X, 2)
    g = length(unique(cluster))
    correction = (g / (g - 1)) * ((n - 1) / (n - k))

    return correction * bread * meat * bread
end

function prais_gls(df, depvar, vars; maxiter = 20, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    X = build_x(df, vars)

    beta, residuals_raw = ols_fit(y, X)
    rho = estimate_rho(residuals_raw, df.diff)

    final_x = nothing
    final_residuals = nothing

    for iter in 1:maxiter
        y_p = prais_transform(y, df.diff, rho)
        X_p = transform_matrix(X, df.diff, rho)

        beta, residuals_p = ols_fit(y_p, X_p)
        residuals_original = y - X * beta
        new_rho = estimate_rho(residuals_original, df.diff)

        final_x = X_p
        final_residuals = residuals_p

        if abs(new_rho - rho) < tolerance
            rho = new_rho
            break
        end

        rho = new_rho
    end

    vcov = cluster_vcov(final_x, final_residuals, df.stateyr)
    se = sqrt.(diag(vcov))

    return (; beta, se, rho)
end

function print_result(title, vars, result)
    println()
    println(title)
    println("rho: ", result.rho)

    for (name, coef, se) in zip(String.(vars), result.beta[1:length(vars)], result.se[1:length(vars)])
        println(rpad(name, 14), coef, "    (", se, ")")
    end
end

vars_emp = [:ioudrg, :mnc_post92, :mnc_post87, :anwage_util, :ln_mwhs, :fgddum]
vars_no_wage = [:ioudrg, :mnc_post92, :mnc_post87, :ln_mwhs, :fgddum]

enf = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
fuel = prepare_input_data(load_frw_data("frw1extract_f.dta"))

println("GLS basic columns for Tables 3, 4, and 5")

result_t3 = prais_gls(enf, :ln_emp, vars_emp)
result_t4 = prais_gls(enf, :ln_nfexp, vars_no_wage)
result_t5 = prais_gls(fuel, :ln_btu, vars_no_wage)

print_result("Table 3 column 1: ln_emp", vars_emp, result_t3)
print_result("Table 4 column 1: ln_nfexp", vars_no_wage, result_t4)
print_result("Table 5 column 1: ln_btu", vars_no_wage, result_t5)