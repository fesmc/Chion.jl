#!/usr/bin/env julia

"""
Run the complete prescribed-versus-parameterized MAR forcing factorial for
climatological Greenland spinups, then compare each member with MAR during one
additional diagnostic year.

The forcing file uses its native EPSG:3413 polar-stereographic x/y coordinates
(km). `AL2` is a two-sector MAR field; sector one is used as the surface
albedo, consistently with Chion's forcing loader convention for two-level
fields.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using Dates: day, month
import CairoMakie
using NCDatasets
using Plots
using Statistics: cor, mean, quantile

const CONFIG = (
    forcing_file="/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc",
    output_dir=get(
        ENV,
        "CHION_OUTPUT_DIR",
        joinpath(@__DIR__, "..", "plots", "gris_prescribed_parameterized", "forcing_factorial"),
    ),
    mask_threshold=50.0,
    spinup_years=parse(Int, get(ENV, "CHION_SPINUP_YEARS", "200")),
    backend=Symbol(get(ENV, "CHION_BACKEND", "gpu")),
    turbulent_flux_scheme=Symbol(get(ENV, "CHION_TURBULENT_FLUX_SCHEME", "semix")),
    semix_surface_height_m=parse(Float64, get(ENV, "CHION_SEMIX_SURFACE_HEIGHT_M", "10.0")),
    semix_z0m_snow_m=parse(Float64, get(ENV, "CHION_SEMIX_Z0M_SNOW_M", "0.001")),
    semix_zm_to_zh=parse(Float64, get(ENV, "CHION_SEMIX_ZM_TO_ZH", "10.0")),
    semix_sensible_exchange_factor=parse(Float64, get(ENV, "CHION_SEMIX_SENSIBLE_EXCHANGE_FACTOR", "2.50")),
    semix_stable_coefficient=parse(Float64, get(ENV, "CHION_SEMIX_STABLE_COEFFICIENT", "20.0")),
    diurnal_substeps=parse(Int, get(ENV, "CHION_DIURNAL_SUBSTEPS", "1")),
    diurnal_temperature_amplitude_c=parse(Float64, get(ENV, "CHION_DIURNAL_TEMPERATURE_AMPLITUDE_C", "1.0")),
    write_figures=lowercase(get(ENV, "CHION_WRITE_FIGURES", "true")) in ("1", "true", "yes"),
    # Hold wind fixed across the factorial instead of using MAR's variable
    # wind field; this isolates the five prescribed/parameterized products.
    wind_speed_m_s=5.0,
    surface_albedo_name="AL2",
    # Keep the scatter PDFs compact: vectorizing every grid-cell point makes
    # the six monthly budget panels unnecessarily large.  This deterministic
    # subsample is representative; panel statistics still use every valid
    # grid-cell value.
    scatter_max_points=20_000,
)

const COMPARISON_VARIABLES = (
    q_sw_net=("Net shortwave flux (W m⁻²)", "q_sw_net"),
    albedo=("Surface albedo (-)", "albedo"),
    q_lw_down=("Incoming longwave flux (W m⁻²)", "q_lw_down"),
    q_sh=("Sensible heat flux (W m⁻²)", "q_sh"),
    q_lh=("Latent heat flux (W m⁻²)", "q_lh"),
)

const BUDGET_VARIABLES = (
    melt="Melt",
    runoff="Runoff",
    smb_ice="Ice SMB",
    refreezing="Refreezing",
    vapor_mass="Vapor mass",
    surface_balance="SMB",
)

# Ice SMB remains in the mass-budget products, but is excluded from the
# MAR--Chion scatter figure at the requested 2 × 5 layout.
const BUDGET_SCATTER_VARIABLES = (:melt, :runoff, :refreezing, :vapor_mass, :surface_balance)

const CUMULATIVE_BUDGET_VARIABLES = (
    :melt,
    :runoff,
    :smb_ice,
    :refreezing,
    :vapor_mass,
    :sublimation,
)

const ENERGY_VARIABLES = (
    shortwave="Net shortwave",
    longwave="Net longwave",
    sensible_heat="Sensible heat",
    latent_heat="Latent heat",
    rain_heat="Rain heat",
)

# Fixed limits make a given field directly comparable among every factorial
# member. Values are annual mmWE.
const END_STATE_MAP_SPECS = (
    (name=:smb_ice, label="Annual ice SMB [mmWE]", colormap=:RdBu, clims=(-2000.0, 2000.0)),
    (name=:runoff, label="Annual runoff [mmWE]", colormap=:batlow, clims=(0.0, 2500.0)),
    (name=:melt, label="Annual melt [mmWE]", colormap=:batlow, clims=(0.0, 2500.0)),
    (name=:refreezing, label="Annual refreezing [mmWE]", colormap=:batlow, clims=(0.0, 1500.0)),
    (name=:sublimation, label="Annual sublimation [mmWE]", colormap=:batlow, clims=(0.0, 500.0)),
)

const END_STATE_DIFFERENCE_LIMITS = (
    smb_ice=(-500.0, 500.0),
    runoff=(-500.0, 500.0),
    melt=(-500.0, 500.0),
    refreezing=(-500.0, 500.0),
    sublimation=(-50.0, 50.0),
)

"""Materialize device arrays before CPU-side diagnostic accumulation and plotting."""
to_host(values) = values isa Array ? values : Array(values)

"""The five MAR products independently switched in the climatology factorial."""
const FORCING_SWITCHES = (:q_sw_net, :albedo, :q_lw_down, :q_sh, :q_lh)

"""P/M identifier in `FORCING_SWITCHES` order, suitable for filenames and Slurm selection."""
forcing_code(prescribed) = join(value ? "p" : "m" for value in prescribed)

# Keep physical and numerical settings fixed: this ensemble attributes
# differences solely to the prescribed/parameterized forcing products.
const FACTORIAL_DIURNAL_SUBSTEPS = CONFIG.diurnal_substeps
const FACTORIAL_NTOT = 15
1 <= FACTORIAL_DIURNAL_SUBSTEPS <= 24 || error("`CHION_DIURNAL_SUBSTEPS` must be between 1 and 24.")
CONFIG.diurnal_temperature_amplitude_c >= 0 ||
    error("`CHION_DIURNAL_TEMPERATURE_AMPLITUDE_C` must be non-negative.")
const CASES = Tuple(
    let prescribed = NamedTuple{FORCING_SWITCHES}(
            ntuple(index -> (mask & (1 << (index - 1))) != 0, length(FORCING_SWITCHES)),
        )
        (
        name=Symbol("climatology_$(forcing_code(values(prescribed)))"),
        title="Climatological forcing factorial ($(forcing_code(values(prescribed))))",
        prescribed=prescribed,
        diurnal_substeps=FACTORIAL_DIURNAL_SUBSTEPS,
        diurnal_temperature_cycle=true,
        diurnal_temperature_amplitude_c=CONFIG.diurnal_temperature_amplitude_c,
        ntot=FACTORIAL_NTOT,
        )
    end for mask in 0:(2^length(FORCING_SWITCHES) - 1)
)

"""Read a MAR field, select its surface sector when needed, and retain grid columns."""
function read_mar_columns(path, name, grid, ntime)
    rows = grid.is .+ (grid.js .- 1) .* length(grid.x)
    return NCDataset(dataset -> read_forcing_columns(dataset, name, rows, ntime, length(grid.y), length(grid.x)), path)
end

"""Read a static MAR grid field and retain the selected snowpack columns."""
function read_mar_static_columns(path, name, grid)
    field = NCDataset(path) do dataset
        raw, dim_names = Chion._read_variable_data(dataset, name)
        Chion._as_y_x(raw, dim_names, length(grid.y), length(grid.x), name)
    end
    return [field[grid.js[column], grid.is[column]] for column in eachindex(grid.js)]
end

"""Return a forcing object with the requested MAR products enabled or disabled."""
function forcing_for_case(forcing, reference, prescribed)
    return SnowpackForcing(
        time_values=forcing.time_values,
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=forcing.wind_speed,
        q_sw_net=reference.q_sw_net,
        has_q_sw_net=prescribed.q_sw_net,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=prescribed.q_lw_down,
        q_sh=forcing.q_sh,
        has_q_sh=prescribed.q_sh,
        q_lh=forcing.q_lh,
        has_q_lh=prescribed.q_lh,
        relative_humidity=forcing.relative_humidity,
        has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure,
        surface_height=forcing.surface_height,
        prescribed_albedo=reference.albedo,
        has_prescribed_albedo=prescribed.albedo,
        latitude_deg=forcing.latitude_deg,
    )
end

function comparison_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(COMPARISON_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer(), month_days=zeros(Float64, 12))
end

function budget_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(BUDGET_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer(), month_days=zeros(Float64, 12))
end

function energy_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(ENERGY_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer())
end

const SH_STABILITY_EDGES_C = [-Inf, -10.0, -5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0, 10.0, Inf]

function sh_stability_buffers()
    n = length(SH_STABILITY_EDGES_C) - 1
    return (weight=zeros(n), mar=zeros(n), chion=zeros(n), mar2=zeros(n),
        chion2=zeros(n), cross=zeros(n), error2=zeros(n), count=zeros(Int, n))
end

function accumulate_sh_stability!(buffers, delta_temperature, mar, chion, area_km2, dt_days)
    for column in eachindex(delta_temperature)
        d, x, y, area = delta_temperature[column], mar[column], chion[column], area_km2[column]
        isfinite(d) && isfinite(x) && isfinite(y) && isfinite(area) && area > 0 || continue
        bin = clamp(searchsortedlast(SH_STABILITY_EDGES_C, d), 1, length(buffers.weight))
        weight = area * dt_days
        buffers.weight[bin] += weight
        buffers.mar[bin] += weight * x
        buffers.chion[bin] += weight * y
        buffers.mar2[bin] += weight * x * x
        buffers.chion2[bin] += weight * y * y
        buffers.cross[bin] += weight * x * y
        buffers.error2[bin] += weight * (y - x)^2
        buffers.count[bin] += 1
    end
    return nothing
end

"""Accumulate a daily set of values into calendar-month means."""
function accumulate_comparison!(buffers, time, dt_days, truth, configured)
    month_index = month(time)
    for name in keys(COMPARISON_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name] .* dt_days
        @views buffers.configured[name][month_index, :] .+= configured[name] .* dt_days
    end
    buffers.month_days[month_index] += dt_days
    return nothing
end

function finalize_comparison(buffers)
    for name in keys(COMPARISON_VARIABLES), month_index in 1:12
        @views buffers.truth[name][month_index, :] ./= buffers.month_days[month_index]
        @views buffers.configured[name][month_index, :] ./= buffers.month_days[month_index]
    end
    yearly_truth = Dict(
        name => vec(sum(buffers.truth[name] .* reshape(buffers.month_days, :, 1); dims=1) ./ sum(buffers.month_days))
        for name in keys(COMPARISON_VARIABLES)
    )
    yearly_configured = Dict(
        name => vec(sum(buffers.configured[name] .* reshape(buffers.month_days, :, 1); dims=1) ./ sum(buffers.month_days))
        for name in keys(COMPARISON_VARIABLES)
    )
    return (; buffers..., yearly_truth, yearly_configured)
end

"""Accumulate interval mass changes into calendar-month totals."""
function accumulate_budget!(buffers, time, dt_days, truth, configured)
    month_index = month(time)
    for name in keys(BUDGET_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name]
        @views buffers.configured[name][month_index, :] .+= configured[name]
    end
    buffers.month_days[month_index] += dt_days
    return nothing
end

function finalize_budget(buffers)
    yearly_truth = Dict(name => vec(sum(buffers.truth[name]; dims=1)) for name in keys(BUDGET_VARIABLES))
    yearly_configured = Dict(name => vec(sum(buffers.configured[name]; dims=1)) for name in keys(BUDGET_VARIABLES))
    return (; buffers..., yearly_truth, yearly_configured)
end

"""Accumulate interval energy changes (J m⁻²) into calendar-month totals."""
function accumulate_energy!(buffers, time, truth, configured)
    month_index = month(time)
    for name in keys(ENERGY_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name]
        @views buffers.configured[name][month_index, :] .+= configured[name]
    end
    return nothing
end

function selected_scatter_points(x, y, max_points)
    valid = findall(isfinite.(x) .& isfinite.(y))
    isempty(valid) && return Float64[], Float64[]
    indices = length(valid) <= max_points ? valid : valid[round.(Int, range(1, length(valid); length=max_points))]
    return x[indices], y[indices]
end

function scatter_panel(x, y, label, period, max_points; unit=nothing, panel_label=nothing)
    valid = isfinite.(x) .& isfinite.(y)
    xmetrics, ymetrics = x[valid], y[valid]
    xplot, yplot = selected_scatter_points(x, y, max_points)
    isempty(xplot) && return plot(title="$label — $period (no valid data)")
    # Plot a representative subset, but calculate skill from every valid
    # grid-cell value so file-size optimization cannot change the diagnostics.
    rmse = sqrt(mean((ymetrics .- xmetrics) .^ 2))
    correlation = length(xmetrics) > 1 ? cor(xmetrics, ymetrics) : NaN
    r_squared = correlation^2
    lower = min(minimum(xplot), minimum(yplot))
    upper = max(maximum(xplot), maximum(yplot))
    pad = max((upper - lower) * 0.04, eps(Float64))
    limits = (lower - pad, upper + pad)
    panel = scatter(
        xplot,
        yplot;
        markersize=1.0,
        markerstrokewidth=0,
        markeralpha=0.16,
        label="Chion",
        xlabel="$period $label (MAR)$(isnothing(unit) ? "" : " [$unit]")",
        ylabel="$period $label (Chion)$(isnothing(unit) ? "" : " [$unit]")",
        title=isnothing(panel_label) ? "" : panel_label,
        titlelocation=:left,
        titlefont=Plots.font(11),
        legend=:bottomright,
        aspect_ratio=:equal,
        xlims=limits,
        ylims=limits,
    )
    plot!(panel, [limits[1], limits[2]], [limits[1], limits[2]]; color=:black, linewidth=1.2, linestyle=:dash, label="1:1")
    stats = length(xplot) > 1 ? "R²=$(round(r_squared; digits=3))\nr=$(round(correlation; digits=3))\nRMSE=$(round(rmse; digits=3))" : "One annual total"
    annotate!(panel, limits[1] + 0.04 * (limits[2] - limits[1]), limits[2] - 0.04 * (limits[2] - limits[1]),
        Plots.text(stats, 8, :left, :top))
    return panel
end

function comparison_figure(comparison, case)
    panels = Any[]
    for (name, (label, _)) in pairs(COMPARISON_VARIABLES)
        case.prescribed[name] && continue
        push!(panels, scatter_panel(vec(comparison.truth[name]), vec(comparison.configured[name]), label, "monthly", CONFIG.scatter_max_points))
        push!(panels, scatter_panel(comparison.yearly_truth[name], comparison.yearly_configured[name], label, "yearly", CONFIG.scatter_max_points))
    end
    isempty(panels) && return plot(title="$(case.title): all variables are prescribed")
    return plot(
        panels...;
        layout=(length(panels) ÷ 2, 2),
        size=(1500, 440 * (length(panels) ÷ 2)),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"""Add one rasterized budget scatter layer to a CairoMakie axis."""
function budget_scatter_panel!(position, x, y, label, period; unit, panel_label)
    valid = isfinite.(x) .& isfinite.(y)
    xmetrics, ymetrics = x[valid], y[valid]
    xplot, yplot = selected_scatter_points(x, y, CONFIG.scatter_max_points)
    axis = CairoMakie.Axis(
        position;
        xlabel="$period $label (MAR) [$unit]",
        ylabel="$period $label (Chion) [$unit]",
        title="",
        aspect=CairoMakie.DataAspect(),
    )
    isempty(xplot) && return axis

    rmse = sqrt(mean((ymetrics .- xmetrics) .^ 2))
    correlation = length(xmetrics) > 1 ? cor(xmetrics, ymetrics) : NaN
    lower = min(minimum(xplot), minimum(yplot))
    upper = max(maximum(xplot), maximum(yplot))
    pad = max((upper - lower) * 0.04, eps(Float64))
    limits = (lower - pad, upper + pad)
    CairoMakie.xlims!(axis, limits...)
    CairoMakie.ylims!(axis, limits...)
    # Only the dense marker series is rasterized in the PDF. Axes, labels,
    # annotations, and the 1:1 reference line remain vector graphics.
    CairoMakie.scatter!(
        axis,
        xplot,
        yplot;
        color=(:dodgerblue, 0.16),
        markersize=3.0,
        rasterize=2,
    )
    CairoMakie.lines!(axis, collect(limits), collect(limits); color=:black, linewidth=1.5, linestyle=:dash)
    stats = "R²=$(round(correlation^2; digits=3))\nr=$(round(correlation; digits=3))\nRMSE=$(round(rmse; digits=3))"
    CairoMakie.text!(
        axis,
        stats;
        position=(limits[1] + 0.04 * (limits[2] - limits[1]), limits[2] - 0.04 * (limits[2] - limits[1])),
        align=(:left, :top),
        fontsize=14,
    )
    return axis
end

"""Budget scatter figure with rasterized point layers and vector annotations."""
function budget_figure(comparison, area_km2, case)
    # 1 mmWE over 1 km² is 10⁻³ Mt. The monthly panel has 12 × ncell points.
    monthly_cell_total(values) = values .* reshape(area_km2, 1, :) .* 1e-3
    # Near-square cells preserve the 1:1 scatter aspect without unused
    # horizontal space; minimal layout gaps keep the 2 × 5 grid compact.
    figure = CairoMakie.Figure(size=(3000, 1200), fontsize=18, figure_padding=(2, 2, 2, 2))
    CairoMakie.colgap!(figure.layout, 2)
    CairoMakie.rowgap!(figure.layout, 4)
    for (index, name) in enumerate(BUDGET_SCATTER_VARIABLES)
        label = BUDGET_VARIABLES[name]
        monthly_panel_label = "($(Char('a' + index - 1)))"
        annual_panel_label = "($(Char('a' + length(BUDGET_SCATTER_VARIABLES) + index - 1)))"
        monthly_mar = monthly_cell_total(comparison.truth[name])
        monthly_chion = monthly_cell_total(comparison.configured[name])
        yearly_mar = vec(sum(monthly_mar; dims=1))
        yearly_chion = vec(sum(monthly_chion; dims=1))
        budget_scatter_panel!(
            figure[1, index], vec(monthly_mar), vec(monthly_chion), label, "Monthly";
            unit="Mt/month", panel_label=monthly_panel_label,
        )
        CairoMakie.Label(
            figure[1, index, CairoMakie.TopLeft()], monthly_panel_label;
            fontsize=20, halign=:left, valign=:bottom, padding=(0, 0, 4, 0),
        )
        budget_scatter_panel!(
            figure[2, index], yearly_mar, yearly_chion, label, "Annual";
            unit="Mt/yr", panel_label=annual_panel_label,
        )
        CairoMakie.Label(
            figure[2, index, CairoMakie.TopLeft()], annual_panel_label;
            fontsize=20, halign=:left, valign=:bottom, padding=(0, 0, 4, 0),
        )
    end
    return figure
end

"""Return annual per-grid-cell budget-scatter skill metrics in Mt/yr."""
function yearly_budget_scatter_metrics(comparison, area_km2, case)
    # 1 mmWE over 1 km² is 10⁻³ Mt.
    yearly_cell_total(values) = vec(sum(values .* reshape(area_km2, 1, :); dims=1)) .* 1e-3
    rows = NamedTuple[]
    for (name, _) in pairs(BUDGET_VARIABLES)
        mar = yearly_cell_total(comparison.truth[name])
        chion = yearly_cell_total(comparison.configured[name])
        selected = isfinite.(mar) .& isfinite.(chion)
        x, y = mar[selected], chion[selected]
        isempty(x) && continue
        push!(rows, (
            variable=name,
            diurnal_steps=case.diurnal_substeps,
            ntot=case.ntot,
            rmse=sqrt(mean((y .- x) .^ 2)),
            correlation=length(x) > 1 ? cor(x, y) : NaN,
        ))
    end
    return rows
end

function write_yearly_budget_scatter_metrics(path, rows)
    open(path, "w") do io
        println(io, "variable,diurnal_steps,ntot,rmse_mt_per_year,correlation")
        for row in rows
            println(io, "$(row.variable),$(row.diurnal_steps),$(row.ntot),$(row.rmse),$(row.correlation)")
        end
    end
    return nothing
end

"""Place an annual-total annotation above, rather than over, plotted series."""
function annotate_integrated_total!(panel, values, annotation)
    lower, upper = extrema(values)
    span = max(upper - lower, max(abs(lower), abs(upper), 1.0) * 0.1)
    padding = 0.18 * span
    ylims!(panel, (lower - 0.03 * span, upper + padding))
    annotate!(panel, 6.5, upper + 0.92 * padding, Plots.text(annotation, 9, :center, :top))
    return panel
end

"""Plot monthly area-integrated budgets in Gt/month, with annual totals in Gt/yr."""
function monthly_integrated_budget_figure(comparison, area_km2, case)
    # 1 mmWE over 1 km² is 10⁶ kg = 10⁻⁶ Gt.
    integrate_monthly(values) = vec(sum(values .* reshape(area_km2, 1, :); dims=2)) .* 1e-6
    panels = Any[]
    for (name, label) in pairs(BUDGET_VARIABLES)
        panel_index = length(panels) + 1
        panel_label = "($(Char('a' + panel_index - 1)))"
        mar = integrate_monthly(comparison.truth[name])
        chion = integrate_monthly(comparison.configured[name])
        mar_yearly = sum(mar)
        chion_yearly = sum(chion)
        labels = ["MAR" "Chion"]
        annotation = "MAR=$(round(mar_yearly; digits=1)) Gt/yr\nChion=$(round(chion_yearly; digits=1)) Gt/yr"
        if name == :smb_ice
            delta_firn = integrate_monthly(comparison.delta_firn)
            delta_firn_yearly = sum(delta_firn)
            chion_surface_smb = chion_yearly + delta_firn_yearly
            annotation = "MAR=$(round(mar_yearly; digits=1)) Gt/yr\nChion=$(round(chion_yearly; digits=1)) Gt/yr\nΔM firn=$(round(delta_firn_yearly; digits=1)) Gt/yr\nChion + ΔM=$(round(chion_surface_smb; digits=1)) Gt/yr"
            panel = plot(
                1:12,
                mar;
                label="MAR surface SMB",
                marker=:circle,
                markerstrokewidth=0,
                linewidth=2.8,
                xlabel="Month",
                ylabel="$label [Gt/month]",
                title=panel_label,
                titlelocation=:left,
                titlefont=Plots.font(11),
                xticks=1:12,
                framestyle=:box,
            )
            plot!(
                panel,
                1:12,
                chion;
                label="Chion ice SMB",
                marker=:square,
                markerstrokewidth=0,
                linewidth=2.8,
            )
            plot!(
                panel,
                1:12,
                chion .+ delta_firn;
                label="Chion ice SMB + ΔM firn",
                marker=:diamond,
                markerstrokewidth=0,
                linewidth=2.8,
            )
            annotate_integrated_total!(panel, vcat(mar, chion, chion .+ delta_firn), annotation)
            push!(panels, panel)
            continue
        end
        push!(
            panels,
            plot(
                1:12,
                hcat(mar, chion);
                label=labels,
                marker=[:circle :square],
                markerstrokewidth=0,
                linewidth=2.8,
                xlabel="Month",
                ylabel="$label [Gt/month]",
                title=panel_label,
                titlelocation=:left,
                titlefont=Plots.font(11),
                xticks=1:12,
                framestyle=:box,
            ),
        )
        annotate_integrated_total!(panels[end], vcat(mar, chion), annotation)
    end
    return plot(
        panels...;
        layout=(2, 3),
        size=(1800, 1150),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"""Plot final-year monthly integrated surface-energy components in PJ."""
function monthly_integrated_energy_figure(comparison, area_km2, case)
    # 1 J m⁻² over 1 km² is 10⁻⁹ PJ.
    # Model fields retain NaNs outside Greenland.  Treat these cells as
    # absent—not as a NaN contribution to a domain total—so the diagnostic
    # ledger and its plot remain defined for every factorial member.
    integrate_monthly(values) = vec(sum(ifelse.(isfinite.(values),
        values .* reshape(area_km2, 1, :), 0.0); dims=2)) .* 1e-9
    panels = Any[]
    for (name, label) in pairs(ENERGY_VARIABLES)
        panel_index = length(panels) + 1
        panel_label = "($(Char('a' + panel_index - 1)))"
        mar = integrate_monthly(comparison.truth[name])
        chion = integrate_monthly(comparison.configured[name])
        annotation = "MAR=$(round(sum(mar); digits=1)) PJ\nChion=$(round(sum(chion); digits=1)) PJ"
        panel = plot(
            1:12,
            hcat(mar, chion);
            label=["MAR" "Chion"],
            marker=[:circle :square],
            markerstrokewidth=0,
            linewidth=2.8,
            xlabel="Month",
            ylabel="$label [PJ]",
            title=panel_label,
            titlelocation=:left,
            titlefont=Plots.font(11),
            xticks=1:12,
            framestyle=:box,
        )
        annotate_integrated_total!(panel, vcat(mar, chion), annotation)
        push!(panels, panel)
    end
    return plot(
        panels...;
        layout=(2, 3),
        size=(1800, 1100),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"Write the monthly energy totals independently of plotting backends."
function write_monthly_integrated_energy(path, comparison, area_km2)
    integrate_monthly(values) = vec(sum(ifelse.(isfinite.(values),
        values .* reshape(area_km2, 1, :), 0.0); dims=2)) .* 1e-9
    open(path, "w") do io
        println(io, "component,month,mar_pj,chion_pj")
        for (name, _) in pairs(ENERGY_VARIABLES)
            mar = integrate_monthly(comparison.truth[name])
            chion = integrate_monthly(comparison.configured[name])
            for month_index in eachindex(mar)
                println(io, "$(name),$(month_index),$(mar[month_index]),$(chion[month_index])")
            end
        end
    end
    return nothing
end

function weighted_sh_statistics(mar, chion, weights)
    valid = isfinite.(mar) .& isfinite.(chion) .& isfinite.(weights) .& (weights .> 0)
    any(valid) || return (mar=NaN, chion=NaN, bias=NaN, rmse=NaN, correlation=NaN, count=0)
    x, y, w = mar[valid], chion[valid], weights[valid]
    total = sum(w)
    mx, my = sum(w .* x) / total, sum(w .* y) / total
    vx = sum(w .* (x .- mx).^2) / total
    vy = sum(w .* (y .- my).^2) / total
    covariance = sum(w .* (x .- mx) .* (y .- my)) / total
    correlation = vx > 0 && vy > 0 ? covariance / sqrt(vx * vy) : NaN
    return (mar=mx, chion=my, bias=my - mx,
        rmse=sqrt(sum(w .* (y .- x).^2) / total), correlation, count=count(valid))
end

function write_sh_diagnostics(prefix, comparison, stability, area_km2, surface_height)
    month_path = "$(prefix)_by_month.csv"
    open(month_path, "w") do io
        println(io, "month,mar_w_m2,chion_w_m2,bias_w_m2,rmse_w_m2,correlation,n")
        for month_index in 1:12
            stats = weighted_sh_statistics(
                view(comparison.truth[:q_sh], month_index, :),
                view(comparison.configured[:q_sh], month_index, :), area_km2,
            )
            println(io, "$(month_index),$(stats.mar),$(stats.chion),$(stats.bias),$(stats.rmse),$(stats.correlation),$(stats.count)")
        end
    end

    elevation_path = "$(prefix)_by_elevation.csv"
    elevation_edges = collect(0.0:500.0:3500.0)
    open(elevation_path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,mar_w_m2,chion_w_m2,bias_w_m2,rmse_w_m2,correlation,n")
        for index in 1:(length(elevation_edges) - 1)
            lower, upper = elevation_edges[index], elevation_edges[index + 1]
            selected = (surface_height .>= lower) .& (surface_height .< upper)
            stats = weighted_sh_statistics(
                comparison.yearly_truth[:q_sh][selected],
                comparison.yearly_configured[:q_sh][selected], area_km2[selected],
            )
            println(io, "$(lower),$(upper),$(stats.mar),$(stats.chion),$(stats.bias),$(stats.rmse),$(stats.correlation),$(stats.count)")
        end
    end

    stability_path = "$(prefix)_by_stability.csv"
    open(stability_path, "w") do io
        println(io, "delta_t_lower_c,delta_t_upper_c,mar_w_m2,chion_w_m2,bias_w_m2,rmse_w_m2,correlation,n")
        for index in eachindex(stability.weight)
            w = stability.weight[index]
            if w > 0
                mx, my = stability.mar[index] / w, stability.chion[index] / w
                vx = max(stability.mar2[index] / w - mx^2, 0.0)
                vy = max(stability.chion2[index] / w - my^2, 0.0)
                covariance = stability.cross[index] / w - mx * my
                correlation = vx > 0 && vy > 0 ? covariance / sqrt(vx * vy) : NaN
                rmse = sqrt(stability.error2[index] / w)
            else
                mx = my = correlation = rmse = NaN
            end
            println(io, "$(SH_STABILITY_EDGES_C[index]),$(SH_STABILITY_EDGES_C[index + 1]),$(mx),$(my),$(my - mx),$(rmse),$(correlation),$(stability.count[index])")
        end
    end
    return (month_path, elevation_path, stability_path)
end

"""Return NaN-separated coastline segments from the MAR land-mask cell edges."""
function coastline_segments(grid; threshold)
    mask = grid.mask .>= threshold
    ny, nx = size(mask)
    dx = length(grid.x) > 1 ? abs(grid.x[2] - grid.x[1]) : 1.0
    dy = length(grid.y) > 1 ? abs(grid.y[2] - grid.y[1]) : 1.0
    xs, ys = Float64[], Float64[]
    function append_edge!(x1, y1, x2, y2)
        append!(xs, (x1, x2, NaN))
        append!(ys, (y1, y2, NaN))
        return nothing
    end
    for j in 1:ny, i in 1:nx
        mask[j, i] || continue
        x, y = grid.x[i], grid.y[j]
        left = x - dx / 2
        right = x + dx / 2
        bottom = y - dy / 2
        top = y + dy / 2
        (i == 1 || !mask[j, i - 1]) && append_edge!(left, bottom, left, top)
        (i == nx || !mask[j, i + 1]) && append_edge!(right, bottom, right, top)
        (j == 1 || !mask[j - 1, i]) && append_edge!(left, bottom, right, bottom)
        (j == ny || !mask[j + 1, i]) && append_edge!(left, top, right, top)
    end
    return xs, ys
end

function map_panel(
    grid,
    values,
    title;
    coastline_level,
    colorbar_label=title,
    clims=nothing,
    clip_quantiles=nothing,
    colormap=nothing,
    show_title=true,
    panel_label=nothing,
)
    field = Chion.scatter_to_grid(Float64.(to_host(values)), grid.js, grid.is, size(grid.mask))
    finite_values = field[isfinite.(field)]
    isempty(finite_values) && return plot(
        title="$title (no finite values)",
        aspect_ratio=:equal,
        axis=false,
        grid=false,
        framestyle=:none,
    )
    nonnegative = all(>=(zero(eltype(finite_values))), finite_values)
    resolved_colormap = isnothing(colormap) ? (nonnegative ? :batlow : :RdBu) : colormap
    limits = if !isnothing(clims)
        clims
    elseif !isnothing(clip_quantiles)
        lower_quantile, upper_quantile = clip_quantiles
        if nonnegative
            (quantile(finite_values, lower_quantile), quantile(finite_values, upper_quantile))
        else
            extent = quantile(abs.(finite_values), upper_quantile)
            (-extent, extent)
        end
    elseif nonnegative
        extrema(finite_values)
    else
        extent = maximum(abs, finite_values)
        (-extent, extent)
    end
    if limits[1] == limits[2]
        padding = max(abs(limits[1]) * 0.05, one(limits[1]))
        limits = (limits[1] - padding, limits[2] + padding)
    end
    panel = heatmap(
        grid.x,
        grid.y,
        field;
        title=show_title ? title : "",
        aspect_ratio=:equal,
        c=resolved_colormap,
        clims=limits,
        colorbar=true,
        colorbar_title=colorbar_label,
        colorbar_tickfontsize=8,
        colorbar_titlefontsize=11,
        axis=false,
        grid=false,
        framestyle=:none,
        margin=0Plots.mm,
    )
    coastline_x, coastline_y = coastline_segments(grid; threshold=coastline_level)
    plot!(panel, coastline_x, coastline_y; color=:black, linewidth=0.35, label=false, primary=false)
    if !isnothing(panel_label)
        dx = length(grid.x) > 1 ? abs(grid.x[2] - grid.x[1]) : 1.0
        dy = length(grid.y) > 1 ? abs(grid.y[2] - grid.y[1]) : 1.0
        x_min, x_max = extrema(grid.x)
        y_min, y_max = extrema(grid.y)
        map_left = x_min - dx / 2
        map_top = y_max + dy / 2
        plot!(panel; xlims=(map_left - 2dx, x_max + dx / 2), ylims=(y_min - dy / 2, map_top + 1.5dy))
        annotate!(panel, map_left, map_top, Plots.text(panel_label, 13, :right, :bottom))
    end
    return panel
end

function end_state_figure(grid, state, annual, case)
    panels = [
        map_panel(
            grid,
            getfield(annual, spec.name),
            spec.label;
            coastline_level=CONFIG.mask_threshold,
            colorbar_label=spec.label,
            colormap=spec.colormap,
            clims=spec.clims,
            show_title=false,
            panel_label="($(Char('a' + index - 1)))",
        )
        for (index, spec) in enumerate(END_STATE_MAP_SPECS)
    ]
    return plot(
        panels...;
        # Two map rows by three columns, sized for the Greenland aspect ratio.
        layout=(2, 3),
        size=(1800, 1050),
        left_margin=4Plots.mm,
        right_margin=1Plots.mm,
        bottom_margin=0Plots.mm,
        top_margin=2Plots.mm,
    )
end

"""Annual MAR budget fields comparable with Chion's end-state diagnostics."""
function annual_mar_end_state(reference, dt_days)
    annual_sum(values) = vec(sum(values .* reshape(dt_days, 1, :); dims=2))
    return (
        smb_ice=annual_sum(reference.budgets.smb_ice),
        runoff=annual_sum(reference.budgets.runoff),
        melt=annual_sum(reference.budgets.melt),
        refreezing=annual_sum(reference.budgets.refreezing),
        # MAR SU is positive for sublimation; the reference vapor field uses
        # the opposite sign to match Chion's deposition-positive convention.
        sublimation=annual_sum(-reference.budgets.vapor_mass),
    )
end

"""Plot Chion minus MAR annual budget maps using Crameri's `bam` colormap."""
function end_state_difference_figure(grid, mar_annual, annual, case)
    panels = [
        map_panel(
            grid,
            getfield(annual, spec.name) .- getfield(mar_annual, spec.name),
            "Chion − MAR $(spec.label)";
            coastline_level=CONFIG.mask_threshold,
            colorbar_label="$(replace(spec.label, " [mmWE]" => "")) Difference [mmWE]",
            colormap=:bam,
            clims=getfield(END_STATE_DIFFERENCE_LIMITS, spec.name),
            show_title=false,
            panel_label="($(Char('a' + index - 1)))",
        )
        for (index, spec) in enumerate(END_STATE_MAP_SPECS)
    ]
    return plot(
        panels...;
        layout=(2, 3),
        size=(1800, 1050),
        left_margin=4Plots.mm,
        right_margin=1Plots.mm,
        bottom_margin=0Plots.mm,
        top_margin=2Plots.mm,
    )
end

"""Save map figures with narrow GR colourbars, without changing other figures."""
function save_map_figure(path, figure)
    original_width = Plots.gr_cbar_width[]
    Plots.gr_cbar_width[] = 0.009
    try
        savefig(figure, path)
    finally
        Plots.gr_cbar_width[] = original_width
    end
    return nothing
end

"""Summarize annual MAR--Chion refreezing by surface-elevation band in Gt."""
function refreezing_elevation_summary(mar_refreezing, chion_refreezing, mar_melt, chion_melt, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    rows = NamedTuple[]
    for (lower, upper) in zip(edges[1:end-1], edges[2:end])
        selected = (surface_height .>= lower) .& (surface_height .< upper)
        isempty(findall(selected)) && continue
        mar = sum(mar_refreezing[selected] .* area_km2[selected]) * 1e-6
        chion = sum(chion_refreezing[selected] .* area_km2[selected]) * 1e-6
        push!(rows, (
            lower=lower,
            upper=upper,
            mar=mar,
            chion=chion,
            deficit=mar - chion,
            mar_melt=sum(mar_melt[selected] .* area_km2[selected]) * 1e-6,
            chion_melt=sum(chion_melt[selected] .* area_km2[selected]) * 1e-6,
        ))
    end
    return rows
end

function write_refreezing_elevation_summary(path, rows)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,mar_refreezing_gt,chion_refreezing_gt,mar_minus_chion_gt,mar_melt_gt,chion_melt_gt")
        for row in rows
            println(io, "$(row.lower),$(row.upper),$(row.mar),$(row.chion),$(row.deficit),$(row.mar_melt),$(row.chion_melt)")
        end
    end
    return nothing
end

function refreezing_difference_figure(grid, mar_refreezing, chion_refreezing, case)
    shared_upper = max(maximum(mar_refreezing), maximum(chion_refreezing))
    shared_limits = (0.0, shared_upper)
    difference = chion_refreezing .- mar_refreezing
    difference_extent = maximum(abs, difference)
    panels = (
        map_panel(grid, mar_refreezing, "MAR annual refreezing";
            coastline_level=CONFIG.mask_threshold, colorbar_label="Annual refreezing [mmWE]", clims=shared_limits),
        map_panel(grid, chion_refreezing, "Chion annual refreezing";
            coastline_level=CONFIG.mask_threshold, colorbar_label="Annual refreezing [mmWE]", clims=shared_limits),
        map_panel(grid, difference, "Chion − MAR annual refreezing";
            coastline_level=CONFIG.mask_threshold, colorbar_label="Chion − MAR annual refreezing [mmWE]", clims=(-difference_extent, difference_extent)),
    )
    return plot(panels...; layout=(1, 3), size=(1800, 600), title="$(case.title): refreezing diagnosis")
end

function annual_budget_start(state)
    return NamedTuple(name => copy(to_host(getfield(state, name))) for name in CUMULATIVE_BUDGET_VARIABLES)
end

function annual_budget_change(state, start)
    return NamedTuple(name => to_host(getfield(state, name)) .- getfield(start, name) for name in keys(start))
end

"""Total solid plus liquid firn mass per column (mmWE)."""
function firn_mass(state)
    # `TransposedLayerMatrix` is a GPU storage wrapper whose elementwise
    # broadcast falls back to scalar indexing. Transfer each logical field
    # first, then form this CPU-side diagnostic.
    return vec(sum(to_host(state.mass) .+ to_host(state.mass_w); dims=1))
end

"""Diagnose the firn state immediately before the 1 May forcing interval."""
function pre_melt_firn_state(state)
    N = to_host(state.N)
    mass = to_host(state.mass)
    mass_w = to_host(state.mass_w)
    density = to_host(state.density)
    temperature = to_host(state.temperature)
    ncol = length(N)
    cold_content = zeros(Float64, ncol)
    pore_capacity = zeros(Float64, ncol)
    liquid_water = zeros(Float64, ncol)
    depth_to_dense_layer = fill(NaN, ncol)
    c = state.c
    dense_threshold = 810.0
    for column in 1:ncol
        depth = 0.0
        for layer in 1:N[column]
            layer_mass = mass[layer, column]
            layer_density = density[layer, column]
            if layer_mass <= 0 || layer_density <= 0
                continue
            end
            cold_content[column] += max(c.T0 - temperature[layer, column], 0.0) * c.ci * layer_mass / c.Lm
            pore_capacity[column] += max(layer_mass / layer_density - layer_mass / c.rho_i, 0.0) * c.rho_w
            liquid_water[column] += mass_w[layer, column]
            if isnan(depth_to_dense_layer[column]) && layer_density >= dense_threshold
                depth_to_dense_layer[column] = depth
            end
            depth += layer_mass / layer_density
        end
    end
    return (
        cold_content=cold_content,
        pore_capacity=pore_capacity,
        liquid_water=liquid_water,
        depth_to_dense_layer=depth_to_dense_layer,
    )
end

function firn_state_elevation_summary(firn, refreezing, surface_height, area_km2)
    edges = collect(0.0:500.0:3500.0)
    rows = NamedTuple[]
    for (lower, upper) in zip(edges[1:end-1], edges[2:end])
        selected = (surface_height .>= lower) .& (surface_height .< upper)
        !any(selected) && continue
        area = area_km2[selected]
        dense_depth = firn.depth_to_dense_layer[selected]
        finite_dense = isfinite.(dense_depth)
        push!(rows, (
            lower=lower,
            upper=upper,
            refreezing_gt=sum(refreezing[selected] .* area) * 1e-6,
            cold_content_gt=sum(firn.cold_content[selected] .* area) * 1e-6,
            pore_capacity_gt=sum(firn.pore_capacity[selected] .* area) * 1e-6,
            liquid_water_gt=sum(firn.liquid_water[selected] .* area) * 1e-6,
            dense_layer_area_fraction=sum(area[finite_dense]) / sum(area),
            dense_layer_depth_m=any(finite_dense) ? sum(dense_depth[finite_dense] .* area[finite_dense]) / sum(area[finite_dense]) : NaN,
        ))
    end
    return rows
end

function write_firn_state_summary(path, rows)
    open(path, "w") do io
        println(io, "elevation_lower_m,elevation_upper_m,chion_refreezing_gt,cold_content_capacity_gt,pore_capacity_gt,liquid_water_gt,dense_layer_area_fraction,mean_dense_layer_depth_m")
        for row in rows
            println(io, "$(row.lower),$(row.upper),$(row.refreezing_gt),$(row.cold_content_gt),$(row.pore_capacity_gt),$(row.liquid_water_gt),$(row.dense_layer_area_fraction),$(row.dense_layer_depth_m)")
        end
    end
    return nothing
end

function pre_melt_firn_figure(grid, firn, case)
    panels = (
        map_panel(grid, firn.cold_content, "1 May cold-content capacity (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, firn.pore_capacity, "1 May pore capacity (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, firn.liquid_water, "1 May liquid water (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, firn.depth_to_dense_layer, "1 May depth to density ≥810 kg m⁻³ (m)"; coastline_level=CONFIG.mask_threshold),
    )
    return plot(panels...; layout=(2, 2), size=(1200, 1100), title="$(case.title): pre-melt firn state")
end

"""Run a spinup, then collect one full year of diagnostic comparisons."""
function run_case(grid, base_forcing, reference, case)
    alpha_ice = parse(Float64, get(ENV, "CHION_ALPHA_ICE", "0.3"))
    0.0 <= alpha_ice <= 1.0 || error("`CHION_ALPHA_ICE` must be between 0 and 1.")
    forcing = forcing_for_case(base_forcing, reference, case.prescribed)
    model = BESSIModel(
        grid;
        Ntot=case.ntot,
        albedo=case.prescribed.albedo ? :prescribed : :dynamic,
        seb_scheme=:bessi,
        turbulent_flux_scheme=CONFIG.turbulent_flux_scheme,
        semix_surface_height=CONFIG.semix_surface_height_m,
        semix_z0m_snow=CONFIG.semix_z0m_snow_m,
        semix_zm_to_zh=CONFIG.semix_zm_to_zh,
        semix_sensible_exchange_factor=CONFIG.semix_sensible_exchange_factor,
        semix_stable_coefficient=CONFIG.semix_stable_coefficient,
        densification=:bessi,
        fresh_snow_density=:constant,
        alpha_ice=alpha_ice,
        diurnal_shortwave_substeps=hasproperty(case, :diurnal_substeps),
        diurnal_shortwave_max_substeps=hasproperty(case, :diurnal_substeps) ? case.diurnal_substeps : 1,
        diurnal_temperature_cycle=hasproperty(case, :diurnal_temperature_cycle) && case.diurnal_temperature_cycle,
        diurnal_temperature_amplitude_c=hasproperty(case, :diurnal_temperature_amplitude_c) ? case.diurnal_temperature_amplitude_c : 5.0,
    )
    spinup = Simulation(model; forcing, years=CONFIG.spinup_years, backend=CONFIG.backend, write_netcdf=false)
    run!(spinup)

    diagnostic = Simulation(model; forcing, state=spinup.now, years=1, backend=CONFIG.backend, write_netcdf=false)
    integrator = init_integrator(diagnostic; io=devnull)
    runtime = integrator.model_runtime.data
    buffers = comparison_buffers(length(grid.js))
    budgets = budget_buffers(length(grid.js))
    energy = energy_buffers(length(grid.js))
    initial_budget = annual_budget_start(runtime.state)
    previous_budget = annual_budget_start(runtime.state)
    previous_firn_mass = firn_mass(runtime.state)
    previous_net_longwave_energy = zeros(Float64, length(grid.js))
    previous_sensible_heat_energy = zeros(Float64, length(grid.js))
    sh_stability = sh_stability_buffers()
    monthly_delta_firn = zeros(Float64, 12, length(grid.js))
    firn_pre_melt = nothing
    c = model.c
    for time_index in eachindex(forcing.time_values)
        # Diagnose against the state entering this forcing interval. Using the
        # post-step temperature would compare MAR's flux at t with Chion's
        # surface state at t + Δt and artificially degrade agreement.
        state = runtime.state
        if month(forcing.time_values[time_index]) == 5 && day(forcing.time_values[time_index]) == 1
            firn_pre_melt = pre_melt_firn_state(state)
        end
        air_temperature = to_host(forcing.air_temperature[:, time_index])
        surface_temperature = to_host(state.Tsrf)
        surface_albedo = to_host(state.albedo)
        shortwave_down = to_host(forcing.shortwave_down[:, time_index])
        wind_speed = to_host(forcing.wind_speed[:, time_index])
        relative_humidity = to_host(forcing.relative_humidity[:, time_index])
        has_relative_humidity = to_host(forcing.has_relative_humidity[:, time_index])
        air_pressure = to_host(forcing.air_pressure[:, time_index])
        parameterized = (
            q_sw_net=(one(eltype(surface_albedo)) .- surface_albedo) .* shortwave_down,
            albedo=surface_albedo,
            q_lw_down=c.σ * c.ϵ_air .* air_temperature .^ 4,
            q_sh=zeros(eltype(surface_temperature), length(surface_temperature)),
            q_lh=[Chion._resolved_turbulent_latent_heat_flux(
                c,
                surface_temperature[column],
                air_temperature[column],
                false,
                zero(eltype(surface_temperature)),
                has_relative_humidity[column],
                relative_humidity[column],
                air_pressure[column],
                wind_speed[column],
            ) for column in eachindex(surface_temperature)],
        )
        truth = (
            q_sw_net=reference.q_sw_net[:, time_index],
            albedo=reference.albedo[:, time_index],
            q_lw_down=to_host(forcing.q_lw_down[:, time_index]),
            q_sh=to_host(forcing.q_sh[:, time_index]),
            q_lh=to_host(forcing.q_lh[:, time_index]),
        )
        dt_seconds = forcing.dt_days[time_index] * c.seconds_per_day
        mar_energy_flux = (
            shortwave=reference.energy.shortwave[:, time_index],
            longwave=reference.energy.longwave[:, time_index],
            sensible_heat=truth.q_sh,
            latent_heat=truth.q_lh,
            rain_heat=to_host(forcing.rainfall_rate[:, time_index]) .* c.cw .* (air_temperature .- c.T0),
        )
        Chion.step_model!(model, diagnostic.now, integrator.model_runtime, runtime.step_fields, time_index)
        net_longwave_energy = to_host(runtime.workspace.net_longwave_energy)
        interval_net_longwave_energy = net_longwave_energy .- previous_net_longwave_energy
        copyto!(previous_net_longwave_energy, net_longwave_energy)
        sensible_heat_energy = to_host(runtime.workspace.sensible_heat_energy)
        actual_sensible_heat_flux =
            (sensible_heat_energy .- previous_sensible_heat_energy) ./ dt_seconds
        copyto!(previous_sensible_heat_energy, sensible_heat_energy)
        resolved_surface_temperature = to_host(runtime.state.Tsrf)
        counterfactual_semix_sensible_heat = map(eachindex(resolved_surface_temperature)) do column
            constant, linear, _, _ = Chion._semix_turbulent_flux_linearized(
                resolved_surface_temperature[column], c, air_temperature[column],
                relative_humidity[column], air_pressure[column], wind_speed[column],
                c.semix_z0m_snow,
            )
            constant - linear * resolved_surface_temperature[column]
        end
        configured = NamedTuple(
            # For prescribed SH, retain the counterfactual SEMIX estimate in
            # the diagnostic tables. For parameterized SH, use the flux
            # actually integrated by the surface solver.
            name => name === :q_sh ?
                    (case.prescribed.q_sh ? counterfactual_semix_sensible_heat : actual_sensible_heat_flux) :
                    case.prescribed[name] ? truth[name] : parameterized[name]
            for name in keys(COMPARISON_VARIABLES)
        )
        accumulate_comparison!(buffers, forcing.time_values[time_index], forcing.dt_days[time_index], truth, configured)
        accumulate_sh_stability!(
            sh_stability, air_temperature .- resolved_surface_temperature,
            truth.q_sh, configured.q_sh, reference.area_km2,
            forcing.dt_days[time_index],
        )
        chion_energy_flux = (
            shortwave=configured.q_sw_net,
            longwave=interval_net_longwave_energy ./ dt_seconds,
            sensible_heat=actual_sensible_heat_flux,
            latent_heat=configured.q_lh,
            rain_heat=to_host(forcing.rainfall_rate[:, time_index]) .* c.cw .* (air_temperature .- c.T0),
        )
        cumulative_budget = NamedTuple(
            name => to_host(getfield(runtime.state, name)) .- getfield(previous_budget, name)
            for name in CUMULATIVE_BUDGET_VARIABLES
        )
        model_budget = (
            cumulative_budget...,
            surface_balance=(to_host(forcing.snowfall_rate[:, time_index]) .+
                             to_host(forcing.rainfall_rate[:, time_index])) .*
                            (forcing.dt_days[time_index] * model.c.seconds_per_day) .-
                            cumulative_budget.runoff .+
                            cumulative_budget.vapor_mass,
        )
        mar_budget = NamedTuple(
            name => reference.budgets[name][:, time_index] .* forcing.dt_days[time_index]
            for name in keys(BUDGET_VARIABLES)
        )
        accumulate_budget!(budgets, forcing.time_values[time_index], forcing.dt_days[time_index], mar_budget, model_budget)
        mar_energy = (
            mar_energy_flux...,
            net_surface_energy=mar_energy_flux.shortwave .+
                               mar_energy_flux.longwave .+
                               mar_energy_flux.sensible_heat .+
                               mar_energy_flux.latent_heat .+
                               mar_energy_flux.rain_heat,
            melt_energy=reference.budgets.melt[:, time_index] .* forcing.dt_days[time_index] .* c.Lm,
        )
        chion_energy = (
            chion_energy_flux...,
            net_surface_energy=chion_energy_flux.shortwave .+
                               chion_energy_flux.longwave .+
                               chion_energy_flux.sensible_heat .+
                               chion_energy_flux.latent_heat .+
                               chion_energy_flux.rain_heat,
            melt_energy=cumulative_budget.melt .* c.Lm,
        )
        interval_energy_truth = NamedTuple(name => mar_energy[name] .* dt_seconds for name in keys(ENERGY_VARIABLES))
        interval_energy_chion = NamedTuple(name => chion_energy[name] .* dt_seconds for name in keys(ENERGY_VARIABLES))
        interval_energy_truth = (; interval_energy_truth..., melt_energy=mar_energy.melt_energy)
        interval_energy_chion = (; interval_energy_chion..., melt_energy=chion_energy.melt_energy)
        accumulate_energy!(energy, forcing.time_values[time_index], interval_energy_truth, interval_energy_chion)
        firn_mass_change = firn_mass(runtime.state) .- previous_firn_mass
        @views monthly_delta_firn[month(forcing.time_values[time_index]), :] .+= firn_mass_change
        for name in CUMULATIVE_BUDGET_VARIABLES
            copyto!(getfield(previous_budget, name), getfield(runtime.state, name))
        end
        copyto!(previous_firn_mass, firn_mass(runtime.state))
    end
    update_diagnostics!(diagnostic.now)
    isnothing(firn_pre_melt) && error("No 1 May pre-melt diagnostic snapshot was found in the forcing calendar.")
    return (
        comparison=finalize_comparison(buffers),
        budgets=(; finalize_budget(budgets)..., delta_firn=monthly_delta_firn),
        energy=energy,
        sh_stability=sh_stability,
        state=diagnostic.now,
        annual=annual_budget_change(runtime.state, initial_budget),
        firn_pre_melt=firn_pre_melt,
    )
end

function main()
    mkpath(CONFIG.output_dir)
    requested_case_names = filter(!isempty, split(get(ENV, "CHION_CASES", ""), ','))
    selected_cases = isempty(requested_case_names) ? CASES : filter(CASES) do case
        String(case.name) in requested_case_names
    end
    isempty(selected_cases) && error("`CHION_CASES` did not select any known case.")
    loaded = load_forcing_file(
        CONFIG.forcing_file;
        time_name="TIME",
        air_temperature_name="TTZ",
        wind_speed_name=nothing,
        wind_default=CONFIG.wind_speed_m_s,
        mask_name="MSK",
        mask_threshold=CONFIG.mask_threshold,
    )
    grid, base_forcing = loaded.grid, loaded.forcing
    albedo = read_mar_columns(CONFIG.forcing_file, CONFIG.surface_albedo_name, grid, length(base_forcing.time_values))
    albedo .= clamp.(albedo, 0.0, 1.0)
    mar_shortwave_up = read_mar_columns(CONFIG.forcing_file, "SWU", grid, length(base_forcing.time_values))
    mar_longwave_up = read_mar_columns(CONFIG.forcing_file, "LWU", grid, length(base_forcing.time_values))
    reference = (
        albedo=albedo,
        # Use MAR's diagnosed net shortwave when it is prescribed. AL2 × SWD
        # is retained only as the independent albedo comparison target.
        q_sw_net=base_forcing.shortwave_down .- mar_shortwave_up,
        budgets=(
            melt=read_mar_columns(CONFIG.forcing_file, "ME", grid, length(base_forcing.time_values)),
            runoff=read_mar_columns(CONFIG.forcing_file, "RU", grid, length(base_forcing.time_values)),
            smb_ice=read_mar_columns(CONFIG.forcing_file, "SMB", grid, length(base_forcing.time_values)),
            refreezing=read_mar_columns(CONFIG.forcing_file, "RZ", grid, length(base_forcing.time_values)),
            vapor_mass=-read_mar_columns(CONFIG.forcing_file, "SU", grid, length(base_forcing.time_values)),
            surface_balance=read_mar_columns(CONFIG.forcing_file, "SMB", grid, length(base_forcing.time_values)),
        ),
        energy=(
            shortwave=base_forcing.shortwave_down .- mar_shortwave_up,
            longwave=base_forcing.q_lw_down .- mar_longwave_up,
        ),
        area_km2=read_mar_static_columns(CONFIG.forcing_file, "AREA", grid),
    )
    mar_annual = annual_mar_end_state(reference, base_forcing.dt_days)

    for case in selected_cases
        @info "Running Greenland comparison" case=case.name spinup_years=CONFIG.spinup_years
        result = run_case(grid, base_forcing, reference, case)
        scatter_path = joinpath(CONFIG.output_dir, "$(case.name)_flux_scatter.pdf")
        budget_path = joinpath(CONFIG.output_dir, "$(case.name)_budget_scatter.pdf")
        integrated_budget_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_budgets.pdf")
        energy_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_energy.pdf")
        state_path = joinpath(CONFIG.output_dir, "$(case.name)_end_state.pdf")
        state_difference_path = joinpath(CONFIG.output_dir, "$(case.name)_end_state_difference.pdf")
        refreezing_path = joinpath(CONFIG.output_dir, "$(case.name)_refreezing_difference.pdf")
        refreezing_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_refreezing_by_elevation.csv")
        firn_state_path = joinpath(CONFIG.output_dir, "$(case.name)_pre_melt_firn_state.pdf")
        firn_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_pre_melt_firn_by_elevation.csv")
        yearly_metrics_path = joinpath(CONFIG.output_dir, "$(case.name)_budget_scatter_yearly_metrics.csv")
        energy_table_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_energy.csv")
        sh_diagnostic_prefix = joinpath(CONFIG.output_dir, "$(case.name)_sensible_heat")
        mar_refreezing = vec(sum(reference.budgets.refreezing .* reshape(base_forcing.dt_days, 1, :); dims=2))
        mar_melt = vec(sum(reference.budgets.melt .* reshape(base_forcing.dt_days, 1, :); dims=2))
        refreezing_rows = refreezing_elevation_summary(
            mar_refreezing,
            result.annual.refreezing,
            mar_melt,
            result.annual.melt,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        firn_rows = firn_state_elevation_summary(
            result.firn_pre_melt,
            result.annual.refreezing,
            base_forcing.surface_height[:, 1],
            reference.area_km2,
        )
        yearly_metrics = yearly_budget_scatter_metrics(result.budgets, reference.area_km2, case)
        write_refreezing_elevation_summary(refreezing_summary_path, refreezing_rows)
        write_firn_state_summary(firn_summary_path, firn_rows)
        write_yearly_budget_scatter_metrics(yearly_metrics_path, yearly_metrics)
        write_monthly_integrated_energy(energy_table_path, result.energy, reference.area_km2)
        sh_diagnostic_paths = write_sh_diagnostics(
            sh_diagnostic_prefix, result.comparison, result.sh_stability,
            reference.area_km2, base_forcing.surface_height[:, 1],
        )
        if CONFIG.write_figures
            savefig(comparison_figure(result.comparison, case), scatter_path)
            CairoMakie.save(budget_path, budget_figure(result.budgets, reference.area_km2, case))
            savefig(monthly_integrated_budget_figure(result.budgets, reference.area_km2, case), integrated_budget_path)
            savefig(monthly_integrated_energy_figure(result.energy, reference.area_km2, case), energy_path)
            save_map_figure(state_path, end_state_figure(grid, result.state, result.annual, case))
            save_map_figure(state_difference_path, end_state_difference_figure(grid, mar_annual, result.annual, case))
            savefig(refreezing_difference_figure(grid, mar_refreezing, result.annual.refreezing, case), refreezing_path)
            savefig(pre_melt_firn_figure(grid, result.firn_pre_melt, case), firn_state_path)
        end
        @info "Refreezing by elevation" rows=refreezing_rows
        @info "Pre-melt firn state by elevation" rows=firn_rows
        @info "Wrote comparison diagnostics" figures=CONFIG.write_figures refreezing_summary_path firn_summary_path yearly_metrics_path energy_table_path sh_diagnostic_paths
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
