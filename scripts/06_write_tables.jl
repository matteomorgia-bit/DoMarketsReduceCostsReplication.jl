using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra
using Printf

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

function build_gls_x(df, vars)
    main = Matrix{Float64}(df[:, vars])
    fixed = fixed_controls(df)
    cons = ones(size(df, 1), 1)
    return hcat(main, fixed, cons)
end

function ols_fit(y, X)
    beta = X \ y
    residuals = y - X * beta
    return beta, residuals
end

function prais_gls(df, depvar, vars; maxiter = 20, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    X = build_gls_x(df, vars)

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

    return (; vars, beta, se, rho)
end

function fit_2sls(y, endog, instrument, exog, fixed, cons)
    z = hcat(exog, instrument, fixed, cons)
    first_beta = z \ endog
    endog_hat = z * first_beta

    xhat = hcat(exog, endog_hat, fixed, cons)
    beta = xhat \ y

    return beta, first_beta, xhat
end

function structural_eps(y, endog, instrument, exog, fixed, cons, beta, first_beta)
    z = hcat(exog, instrument, fixed, cons)
    mhat = z * first_beta

    x_actual = hcat(exog, endog, fixed, cons)
    que = y - x_actual * beta

    beta_endog = beta[size(exog, 2) + 1]
    return que - (endog - mhat) .* beta_endog
end

function prais_iv(df, depvar, exog_vars; maxiter = 7, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    endog = Float64.(df.ln_mwhs)
    instrument = reshape(Float64.(df.lstsales), :, 1)
    exog = Matrix{Float64}(df[:, exog_vars])
    fixed = fixed_controls(df)
    cons = ones(size(df, 1), 1)

    beta, first_beta, _ = fit_2sls(y, endog, instrument, exog, fixed, cons)
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

        beta, first_beta, xhat_p = fit_2sls(
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

    vcov = cluster_vcov(final_xhat, final_structural_residuals, df.stateyr)
    se = sqrt.(diag(vcov))

    return (; vars = vcat(exog_vars, [:ln_mwhs]), beta, se, rho)
end

function fmt(x)
    return @sprintf("%.3f", x)
end

function stars(beta, se)
    t = abs(beta / se)
    if t >= 2.576
        return "***"
    elseif t >= 1.960
        return "**"
    elseif t >= 1.645
        return "*"
    else
        return ""
    end
end

function cell(result, var)
    idx = findfirst(==(var), result.vars)
    isnothing(idx) && return ""
    return "$(fmt(result.beta[idx]))$(stars(result.beta[idx], result.se[idx]))<br>($(fmt(result.se[idx])))"
end

function rho_cell(result)
    return @sprintf("%.2f", result.rho)
end

function write_markdown_table(path, title, rows, results)
    open(path, "w") do io
        println(io, "### ", title)
        println(io)
        println(io, "| Variable | GLS Basic | GLS-IV Basic | GLS-IV Law Date | GLS-IV Retail Access | GLS-IV Nonutility Generation |")
        println(io, "|---|---:|---:|---:|---:|---:|")

        for (label, vars) in rows
            vals = [cell(results[col], var) for (col, var) in zip(keys(results), vars)]
            println(io, "| ", label, " | ", join(vals, " | "), " |")
        end

        rhos = [rho_cell(results[col]) for col in keys(results)]
        println(io, "| rho | ", join(rhos, " | "), " |")
        println(io)
        println(io, "Notes: cluster-robust standard errors are in parentheses. Stars are computed from absolute t-ratios: *** 1 percent, ** 5 percent, * 10 percent. Very small last-digit differences from the printed article can reflect rounding/truncation in the PDF/Stata output.")
    end
end

mkpath("tables")

enf = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
fuel = prepare_input_data(load_frw_data("frw1extract_f.dta"))

vars_emp_gls = [:ioudrg, :mnc_post92, :mnc_post87, :anwage_util, :ln_mwhs, :fgddum]
vars_nowage_gls = [:ioudrg, :mnc_post92, :mnc_post87, :ln_mwhs, :fgddum]

t3 = (
    gls = prais_gls(enf, :ln_emp, vars_emp_gls),
    iv_basic = prais_iv(enf, :ln_emp, [:ioudrg, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
    iv_law = prais_iv(enf, :ln_emp, [:ioulaw, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
    iv_retail = prais_iv(enf, :ln_emp, [:ioudrg, :iouretail, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
    iv_nug = prais_iv(enf, :ln_emp, [:iounug, :mnc_post92, :mnc_post87, :anwage_util, :fgddum]),
)

t4 = (
    gls = prais_gls(enf, :ln_nfexp, vars_nowage_gls),
    iv_basic = prais_iv(enf, :ln_nfexp, [:ioudrg, :mnc_post92, :mnc_post87, :fgddum]),
    iv_law = prais_iv(enf, :ln_nfexp, [:ioulaw, :mnc_post92, :mnc_post87, :fgddum]),
    iv_retail = prais_iv(enf, :ln_nfexp, [:ioudrg, :iouretail, :mnc_post92, :mnc_post87, :fgddum]),
    iv_nug = prais_iv(enf, :ln_nfexp, [:iounug, :mnc_post92, :mnc_post87, :fgddum]),
)

t5 = (
    gls = prais_gls(fuel, :ln_btu, vars_nowage_gls),
    iv_basic = prais_iv(fuel, :ln_btu, [:ioudrg, :mnc_post92, :mnc_post87, :fgddum]),
    iv_law = prais_iv(fuel, :ln_btu, [:ioulaw, :mnc_post92, :mnc_post87, :fgddum]),
    iv_retail = prais_iv(fuel, :ln_btu, [:ioudrg, :iouretail, :mnc_post92, :mnc_post87, :fgddum]),
    iv_nug = prais_iv(fuel, :ln_btu, [:iounug, :mnc_post92, :mnc_post87, :fgddum]),
)

table3_rows = [
    ("IOU * Restructured", [:ioudrg, :ioudrg, :ioudrg, :ioudrg, :ioudrg]),
    ("IOU * Law Passed", [:ioulaw, :ioulaw, :ioulaw, :ioulaw, :ioulaw]),
    ("IOU * Retail Access", [:iouretail, :iouretail, :iouretail, :iouretail, :iouretail]),
    ("IOU * High Nonutility Generation", [:iounug, :iounug, :iounug, :iounug, :iounug]),
    ("MUNI * Post 1992", [:mnc_post92, :mnc_post92, :mnc_post92, :mnc_post92, :mnc_post92]),
    ("MUNI * Post 1987", [:mnc_post87, :mnc_post87, :mnc_post87, :mnc_post87, :mnc_post87]),
    ("ln(WAGE)", [:anwage_util, :anwage_util, :anwage_util, :anwage_util, :anwage_util]),
    ("ln(NET MWH)", [:ln_mwhs, :ln_mwhs, :ln_mwhs, :ln_mwhs, :ln_mwhs]),
    ("SCRUBBER", [:fgddum, :fgddum, :fgddum, :fgddum, :fgddum]),
]

table_nowage_rows = [
    ("IOU * Restructured", [:ioudrg, :ioudrg, :ioudrg, :ioudrg, :ioudrg]),
    ("IOU * Law Passed", [:ioulaw, :ioulaw, :ioulaw, :ioulaw, :ioulaw]),
    ("IOU * Retail Access", [:iouretail, :iouretail, :iouretail, :iouretail, :iouretail]),
    ("IOU * High Nonutility Generation", [:iounug, :iounug, :iounug, :iounug, :iounug]),
    ("MUNI * Post 1992", [:mnc_post92, :mnc_post92, :mnc_post92, :mnc_post92, :mnc_post92]),
    ("MUNI * Post 1987", [:mnc_post87, :mnc_post87, :mnc_post87, :mnc_post87, :mnc_post87]),
    ("ln(NET MWH)", [:ln_mwhs, :ln_mwhs, :ln_mwhs, :ln_mwhs, :ln_mwhs]),
    ("SCRUBBER", [:fgddum, :fgddum, :fgddum, :fgddum, :fgddum]),
]

write_markdown_table("tables/table3.md", "Table 3: Labor Input Demand", table3_rows, pairs(t3))
write_markdown_table("tables/table4.md", "Table 4: Nonfuel Expense Input Demand", table_nowage_rows, pairs(t4))
write_markdown_table("tables/table5.md", "Table 5: Fuel Input Demand", table_nowage_rows, pairs(t5))

println("Wrote generated tables to tables/table3.md, tables/table4.md, and tables/table5.md")
