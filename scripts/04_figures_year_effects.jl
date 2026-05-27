using CairoMakie
using DataFrames
using DoMarketsReduceCostsReplication
using LinearAlgebra

const YEARS = 1981:1999
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
        for year in YEARS
            push!(cols, Float64.((df.reg_group .== group) .& (df.yr_data .== year)))
            push!(names, (group, year))
        end
    end

    return hcat(cols...), names
end

function build_figure_x(df, controls)
    main = Matrix{Float64}(df[:, controls])
    group_years, group_year_names = group_year_matrix(df)
    plants, _ = dummy_matrix(df.plant_num2)
    constant = ones(size(df, 1), 1)

    X = hcat(main, group_years, plants, constant)
    return X, group_year_names, size(main, 2)
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

function prais_group_year(df, depvar, controls; maxiter = 20, tolerance = 0.005)
    y = Float64.(df[!, depvar])
    X, group_year_names, n_controls = build_figure_x(df, controls)

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

    first_group_year = n_controls + 1
    last_group_year = n_controls + length(group_year_names)
    effects = beta[first_group_year:last_group_year]

    out = DataFrame(
        group = first.(group_year_names),
        year = last.(group_year_names),
        effect = effects,
    )

    # Normalize each group relative to 1981, as in the paper figures.
    for group in GROUPS
        baseline = only(out.effect[(out.group .== group) .& (out.year .== 1981)])
        out.effect[out.group .== group] .-= baseline
    end

    return out, rho
end

function plot_effects(effects, title, ylabel, outfile)
    fig = Figure(size = (980, 620), backgroundcolor = :white)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Year",
        ylabel = ylabel,
        xticks = 1982:2:1998,
        xgridcolor = (:gray80, 0.45),
        ygridcolor = (:gray80, 0.45),
        topspinevisible = false,
        rightspinevisible = false,
    )

    colors = Dict(
        "MUNI" => "#3366AA",
        "IOU non-restructured" => "#109E73",
        "IOU restructured" => "#CC5A00",
    )

    markers = Dict(
        "MUNI" => :circle,
        "IOU non-restructured" => :rect,
        "IOU restructured" => :utriangle,
    )

    labels = Dict(
        "MUNI" => "MUNI plants",
        "IOU non-restructured" => "IOU: non-restructured states",
        "IOU restructured" => "IOU: restructured states",
    )

    for group in GROUPS
        sub = effects[effects.group .== group, :]
        lines!(ax, sub.year, sub.effect, linewidth = 3, color = colors[group], label = labels[group])
        scatter!(ax, sub.year, sub.effect, marker = markers[group], markersize = 9, color = colors[group])
    end

    hlines!(ax, [0], color = (:gray35, 0.8), linestyle = :dash, linewidth = 1.5)
    xlims!(ax, 1981, 1999)

    low = minimum(effects.effect)
    high = maximum(effects.effect)
    pad = 0.08 * (high - low)
    ylims!(ax, min(low, 0) - pad, max(high, 0) + pad)

    Legend(
        fig[2, 1],
        ax,
        orientation = :horizontal,
        framevisible = false,
        tellwidth = false,
        nbanks = 1,
    )

    rowgap!(fig.layout, 8)
    save(outfile, fig)

    return fig
end

mkpath(joinpath("output", "figures"))
mkpath("images")

df = prepare_input_data(load_frw_data("frw1extract_enf.dta"))
add_regulatory_group!(df)

emp_effects, rho_emp = prais_group_year(
    df,
    :ln_emp,
    [:mnc_post87, :anwage_util, :ln_mwhs, :fgddum],
)

nf_effects, rho_nf = prais_group_year(
    df,
    :ln_nfexp,
    [:mnc_post87, :ln_mwhs, :fgddum],
)

println("Figure 1 rho: ", rho_emp)
println("Figure 2 rho: ", rho_nf)

fig1 = plot_effects(
    emp_effects,
    "Labor Input Demand Year Effects",
    "Year effect relative to 1981",
    joinpath("output", "figures", "figure1_labor_year_effects.png"),
)

fig2 = plot_effects(
    nf_effects,
    "Nonfuel Expense Input Demand Year Effects",
    "Year effect relative to 1981",
    joinpath("output", "figures", "figure2_nonfuel_year_effects.png"),
)

save(joinpath("images", "figure1_labor_year_effects.png"), fig1)
save(joinpath("images", "figure2_nonfuel_year_effects.png"), fig2)

println("Saved figures to output/figures and images.")
