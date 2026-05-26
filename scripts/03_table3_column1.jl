using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra

const MAIN_VARS = [
    :ioudrg,
    :mnc_post92,
    :mnc_post87,
    :anwage_util,
    :ln_mwhs,
    :fgddum,
]

const YEAR_VARS = Symbol.("year", 82:99)

function dummy_matrix(values)
    levels = sort(unique(values))
    levels = levels[2:end]  # drop first plant dummy, as Stata drops one collinear dummy
    return hcat([Float64.(values .== level) for level in levels]...)
end

function build_x(df)
    main = Matrix{Float64}(df[:, MAIN_VARS])
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

df = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
y = Float64.(df.ln_emp)
X = build_x(df)

beta, residuals_raw = ols_fit(y, X)
rho = estimate_rho(residuals_raw, df.diff)

println("Table 3, column 1 matrix Prais-Winsten")
println("Observations: ", size(df, 1))
println("Plant-epochs: ", count_plant_epochs(df))
println("Initial rho: ", rho)

for iter in 1:20
    y_prais = prais_transform(y, df.diff, rho)
    X_prais = transform_matrix(X, df.diff, rho)

    global beta, residuals_prais = ols_fit(y_prais, X_prais)
    residuals_original = y - X * beta
    new_rho = estimate_rho(residuals_original, df.diff)

    println("Iteration ", iter, " rho: ", new_rho)

    if abs(new_rho - rho) < 0.005
        global rho = new_rho
        break
    end

    global rho = new_rho
end

names_main = String.(MAIN_VARS)

println()
println("Final rho: ", rho)
println("Main coefficients:")
for (name, coef) in zip(names_main, beta[1:length(MAIN_VARS)])
    println(rpad(name, 14), coef)
end