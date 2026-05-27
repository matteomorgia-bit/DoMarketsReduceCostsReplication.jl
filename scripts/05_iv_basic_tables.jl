using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra

const YEAR_VARS = Symbol.("year", 82:99)

function dummy_matrix(values)
    levels = sort(unique(values))
    levels = levels[2:end]
    return hcat([Float64.(values .== level) for level in levels]...)
end

function fixed_controls(df)
    years = Matrix{Float64}(df[:, YEAR_VARS])
    plants = dummy_matrix(df.plant_num2)
    return hcat(years, plants)
end

function transform_matrix(X, diff, rho)
    out = similar(X, Float64)

    for j in axes(X, 2)
        out[:, j] = prais_transform(X[:, j], diff, rho)
    end

    return out
end

function fit_2sls(y, endog, instrument, exog, fixed, cons)
    z = hcat(exog, instrument, fixed, cons)
    first_beta = z \ endog
    endog_hat = z * first_beta

    xhat = hcat(exog, endog_hat, fixed, cons)
    beta = xhat \ y

    return beta, first_beta, xhat, endog_hat
end

function structural_eps(y, endog, instrument, exog, fixed, cons, beta, first_beta)
    z = hcat(exog, instrument, fixed, cons)
    mhat = z * first_beta

    x_actual = hcat(exog, endog, fixed, cons)
    que = y - x_actual * beta

    beta_endog = beta[size(exog, 2) + 1]

    return que - (endog - mhat) .* beta_endog
end

function cluster_vcov_iv(xhat, residuals, cluster)
    bread = inv(xhat' * xhat)
    meat = zeros(size(xhat, 2), size(xhat, 2))

    for g in unique(cluster)
        idx = findall(cluster .== g)
        xg = xhat[idx, :]
        ug = residuals[idx]
        score = xg' * ug
        meat .+= score * score'
    end

    n = size(xhat, 1)
    k = size(xhat, 2)
    g = length(unique(cluster))
    correction = (g / (g - 1)) * ((n - 1) / (n - k))

    return correction * bread * meat * bread
end

function prais_iv(df, depvar, exog_vars; maxiter = 7, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    endog = Float64.(df.ln_mwhs)
    instrument = reshape(Float64.(df.lstsales), :, 1)
    exog = Matrix{Float64}(df[:, exog_vars])
    fixed = fixed_controls(df)
    cons = ones(size(df, 1), 1)

    beta, first_beta, _, _ = fit_2sls(y, endog, instrument, exog, fixed, cons)
    eps = structural_eps(y, endog, instrument, exog, fixed, cons, beta, first_beta)
    rho = estimate_rho(eps, df.diff)

    final_xhat = nothing
    final_structural_residuals = nothing

    for iter in 1:maxiter
        y_p = prais_transform(y, df.diff, rho)
        endog_p = prais_transform(endog, df.diff, rho)
        instrument_p = reshape(prais_transform(vec(instrument), df.diff, rho), :, 1)
        exog_p = transform_matrix(exog, df.diff, rho)
        fixed_p = transform_matrix(fixed, df.diff, rho)
        cons_p = reshape(prais_transform(vec(cons), df.diff, rho), :, 1)

        beta, first_beta, xhat_p, _ = fit_2sls(
            y_p,
            endog_p,
            instrument_p,
            exog_p,
            fixed_p,
            cons_p,
        )

        x_actual_p = hcat(exog_p, endog_p, fixed_p, cons_p)
        structural_residuals_p = y_p - x_actual_p * beta

        final_xhat = xhat_p
        final_structural_residuals = structural_residuals_p

        eps = structural_eps(y, endog, instrument, exog, fixed, cons, beta, first_beta)
        new_rho = estimate_rho(eps, df.diff)

        if abs(new_rho - rho) < tolerance
            rho = new_rho
            break
        end

        rho = new_rho
    end

    vcov = cluster_vcov_iv(final_xhat, final_structural_residuals, df.stateyr)
    se = sqrt.(diag(vcov))

    return (; beta, se, rho)
end

function run_table(df, depvar, specs)
    results = Dict{String, Any}()

    for (name, exog_vars) in specs
        results[name] = prais_iv(df, depvar, exog_vars)
    end

    return results
end

function print_result(table, spec, vars, result)
    names = vcat(String.(vars), ["ln_mwhs"])

    println()
    println(table, " ", spec)
    println("rho: ", result.rho)

    for (name, coef, se) in zip(names, result.beta[1:length(names)], result.se[1:length(names)])
        println(rpad(name, 14), coef, "    (", se, ")")
    end
end

function print_table_results(table, results, specs)
    for (spec, vars) in specs
        print_result(table, spec, vars, results[spec])
    end
end

specs_emp = [
    ("col2_basic", [:ioudrg, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
    ("col3_law", [:ioulaw, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
    ("col4_retail", [:ioudrg, :iouretail, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
    ("col5_nug", [:iounug, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
]

specs_no_wage = [
    ("col2_basic", [:ioudrg, :mnc_post92, :mnc_post87, :fgddum]),
    ("col3_law", [:ioulaw, :mnc_post92, :mnc_post87, :fgddum]),
    ("col4_retail", [:ioudrg, :iouretail, :mnc_post92, :mnc_post87, :fgddum]),
    ("col5_nug", [:iounug, :mnc_post92, :mnc_post87, :fgddum]),
]

enf = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
fuel = prepare_input_data(load_frw_data("frw1extract_f.dta"))

println("GLS-IV columns 2-5 for Tables 3, 4, and 5")

table3 = run_table(enf, :ln_emp, specs_emp)
table4 = run_table(enf, :ln_nfexp, specs_no_wage)
table5 = run_table(fuel, :ln_btu, specs_no_wage)

print_table_results("Table 3", table3, specs_emp)
print_table_results("Table 4", table4, specs_no_wage)
print_table_results("Table 5", table5, specs_no_wage)