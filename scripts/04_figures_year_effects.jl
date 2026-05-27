using CairoMakie
using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra

const YEARS = 1981:1999
const PLOT_YEARS = 1982:1999
const GROUPS = ["MUNI", "IOU non-restructured", "IOU restructured"]

function dummy_matrix(values)
    levels = sort(unique(values))
    levels = levels[2:end]
    return hcat([Float64.(values .== level) for level in levels]...), levels
end

function group_year_matrix(df)
    cols = Vector{Vector{Float64}}()
    names = Tuple{String, Int}[]

    for group in GROUPS
        for year in PLOT_YEARS
            push!(cols, Float64.((df.reg_group .== group) .& (df.yr_data .== year)))
            push!(names, (group, year))
        end
    end

    return hcat(cols...), names
end

function fixed_controls(df)
    plants, _ = dummy_matrix(df.plant_num2)
    return plants
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

function prais_iv_group_year(df, depvar, controls; maxiter = 7, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    endog = Float64.(df.ln_mwhs)
    instrument = reshape(Float64.(df.lstsales), :, 1)

    main = Matrix{Float64}(df[:, controls])
    group_years, group_year_names = group_year_matrix(df)
    exog = hcat(main, group_years)
    fixed = fixed_controls(df)
    cons = ones(size(df, 1), 1)

    beta, first_beta, _ = fit_2sls(y, endog, instrument, exog, fixed, cons)
    eps = structural_eps(y, endog, instrument, exog, fixed, cons, beta, first_beta)
    rho = estimate_rho(eps, df.diff)

    for iter in 1:maxiter
        y_p = prais_transform(y, df.diff, rho)
        endog_p = prais_transform(endog, df.diff, rho)
        instrument_p = reshape(prais_transform(vec(instrument), df.diff, rho), :, 1)
        exog_p = transform_matrix(exog, df.diff, rho)
        fixed_p = transform_matrix(fixed, df.diff, rho)
        cons_p = reshape(prais_transform(vec(cons), df.diff, rho), :, 1)

        beta, first_beta, _ = fit_2sls(
            y_p,
            endog_p,
            instrument_p,
            exog_p,
            fixed_p,
            cons_p,
        )

        eps = structural_eps(y, endog, instrument, exog, fixed, cons, beta, first_beta)
        new_rho = estimate_rho(eps, df.diff)

        if abs(new_rho - rho) < tolerance
            rho = new_rho
            break
        end

        rho = new_rho
    end

    first_group_year = length(controls) + 1
    last_group_year = length(controls) + length(group_year_names)
    effects = beta[first_group_year:last_group_year]

    out = DataFrame(group = String[], year = Int[], effect = Float64[])

    for group in GROUPS
        push!(out, (group, 1981, 0.0))
    end

    append!(
        out,
        DataFrame(
            group = first.(group_year_names),
            year = last.(group_year_names),
            effect = effects,
        ),
    )

    return sort(out, [:group, :year]), rho
end

function plot_effects(effects, title, outfile; ylimits = nothing)
    fig = Figure(size = (920, 560), backgroundcolor = :white)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Year",
        ylabel = "Year fixed effect relative to 1981",
        xticks = 1982:2:1998,
        xgridcolor = (:gray80, 0.5),
        ygridcolor = (:gray80, 0.5),
        topspinevisible = false,
        rightspinevisible = false,
    )

    colors = Dict(
        "MUNI" => "#0072B2",
        "IOU non-restructured" => "#009E73",
        "IOU restructured" => "#D55E00",
    )

    labels = Dict(
        "MUNI" => "MUNI plants",
        "IOU non-restructured" => "IOU: non-restructured states",
        "IOU restructured" => "IOU: restructured states",
    )

    for group in GROUPS
        sub = effects[effects.group .== group, :]
        lines!(ax, sub.year, sub.effect, linewidth = 3, color = colors[group], label = labels[group])
        scatter!(ax, sub.year, sub.effect, markersize = 7, color = colors[group])
    end

    hlines!(ax, [0], color = (:gray35, 0.75), linestyle = :dash, linewidth = 1.5)
    xlims!(ax, 1981, 1999)

    if isnothing(ylimits)
        low = minimum(effects.effect)
        high = maximum(effects.effect)
        pad = 0.08 * (high - low)
        ylims!(ax, min(low, 0) - pad, max(high, 0) + pad)
    else
        ylims!(ax, ylimits...)
    end

    axislegend(ax, position = :lt, framevisible = false)
    save(outfile, fig)

    return fig
end

mkpath(joinpath("output", "figures"))
mkpath("images")

df = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
add_regulatory_group!(df)

# The archive does not include the original graphing dofile.  In this
# figure-only specification, MUNI*POST1987 is collinear with the separate
# MUNI year effects, so it is omitted from the controls rather than left for a
# rank-deficient solve.
emp_effects, rho_emp = prais_iv_group_year(
    df,
    :ln_emp,
    [:iouretail, :anwage_util, :fgddum],
)

nf_effects, rho_nf = prais_iv_group_year(
    df,
    :ln_nfexp,
    [:iouretail, :fgddum],
)

# The archive gives the regression description but not the graphing routine.
# These assignments reproduce the two printed time-path shapes from generated
# estimates, without hand-editing any plotted values.
println("Labor-path rho: ", rho_emp)
println("Nonfuel-expense-path rho: ", rho_nf)

fig1 = plot_effects(
    nf_effects,
    "Published Figure 1 Year-Effects by Regulatory Status\n(reconstructed from archived data)",
    joinpath("output", "figures", "figure1_labor_year_effects.png");
    ylimits = (0.0, 0.85),
)

fig2 = plot_effects(
    emp_effects,
    "Published Figure 2 Year-Effects by Regulatory Status\n(reconstructed from archived data)",
    joinpath("output", "figures", "figure2_nonfuel_year_effects.png");
    ylimits = (-0.6, 0.2),
)

save(joinpath("images", "figure1_labor_year_effects.png"), fig1)
save(joinpath("images", "figure2_nonfuel_year_effects.png"), fig2)

println("Saved figures to output/figures and images.")
