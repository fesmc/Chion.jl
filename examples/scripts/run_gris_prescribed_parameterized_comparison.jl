#!/usr/bin/env julia

"""
Run four non-SEMIX BESSI Greenland spinups that selectively prescribe MAR
surface fluxes, then compare the effective Chion fluxes with MAR during one
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
using NCDatasets
using Plots
using Statistics: cor, mean

const CONFIG = (
    forcing_file="/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc",
    output_dir=joinpath(@__DIR__, "..", "plots", "gris_prescribed_parameterized"),
    mask_threshold=50.0,
    spinup_years=200,
    backend=:threads,
    Ntot=20,
    surface_albedo_name="AL2",
    scatter_max_points=100_000,
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
    vapor_mass="Vapor mass: deposition (+), sublimation (−)",
    surface_balance="Surface balance: snowfall + rainfall − runoff + vapor mass",
)

const CUMULATIVE_BUDGET_VARIABLES = (
    :melt,
    :runoff,
    :smb_ice,
    :refreezing,
    :vapor_mass,
    :sublimation,
)

const ENERGY_VARIABLES = (
    net_surface_energy="Net surface energy",
    shortwave="Net shortwave",
    longwave="Net longwave",
    sensible_heat="Sensible heat",
    latent_heat="Latent heat",
    rain_heat="Rain heat",
    melt_energy="Melt energy",
)

const CASES = (
    (
        name=:all_prescribed_1_step_50layers,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
        ntot=50,
    ),
    (
        name=:all_prescribed_4_step_50layers,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
        diurnal_substeps=4,
        diurnal_temperature_cycle=true,
        diurnal_temperature_amplitude_c=5.0,
        ntot=50
    ),
    (
        name=:all_prescribed,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
    ),
    (
        name=:all_prescribed_4_step_false,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
        diurnal_substeps=4,
        diurnal_temperature_cycle=false,
        diurnal_temperature_amplitude_c=5.0
    ),
        (
        name=:all_prescribed_4_step_15,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
        diurnal_substeps=4,
        diurnal_temperature_cycle=true,
        diurnal_temperature_amplitude_c=15.0
    ),
            (
        name=:all_prescribed_12_step,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
        diurnal_substeps=12,
        diurnal_temperature_cycle=true,
        diurnal_temperature_amplitude_c=5.0
    ),
    (
        name=:all_prescribed_24_step,
        title="All prescribed",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=true, q_sh=true, q_lh=true),
        diurnal_substeps=24,
        diurnal_temperature_cycle=true,
        diurnal_temperature_amplitude_c=5.0
    ),
    (
        name=:age_albedo,
        title="Age-only albedo parameterized",
        prescribed=(q_sw_net=true, albedo=false, q_lw_down=true, q_sh=true, q_lh=true),
    ),
    (
        name=:fluxes_parameterized,
        title="Longwave and turbulent fluxes parameterized",
        prescribed=(q_sw_net=true, albedo=true, q_lw_down=false, q_sh=false, q_lh=false),
    ),
    (
        name=:all_parameterized_diurnal_1,
        title="All parameterized (1 diurnal substep)",
        prescribed=(q_sw_net=false, albedo=false, q_lw_down=false, q_sh=false, q_lh=false),
        diurnal_substeps=1,
    ),
    (
        name=:all_parameterized_diurnal_2,
        title="All parameterized (2 diurnal substeps)",
        prescribed=(q_sw_net=false, albedo=false, q_lw_down=false, q_sh=false, q_lh=false),
        diurnal_substeps=2,
    ),
    (
        name=:all_parameterized_diurnal_4,
        title="All parameterized (4 diurnal substeps)",
        prescribed=(q_sw_net=false, albedo=false, q_lw_down=false, q_sh=false, q_lh=false),
        diurnal_substeps=4,
    ),
    (
        name=:all_parameterized_diurnal_12,
        title="All parameterized (12 diurnal substeps)",
        prescribed=(q_sw_net=false, albedo=false, q_lw_down=false, q_sh=false, q_lh=false),
        diurnal_substeps=12,
    ),
)

"""Read a MAR field, select its surface sector when needed, and retain grid columns."""
function read_mar_columns(path, name, grid, ntime)
    data = NCDataset(path) do dataset
        raw, dim_names = Chion._read_variable_data(dataset, name)
        ny, nx = length(grid.y), length(grid.x)
        Chion._column_matrix(Chion._as_time_y_x(raw, dim_names, ntime, ny, nx, name))
    end
    rows = grid.is .+ (grid.js .- 1) .* length(grid.x)
    return data[rows, :]
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
    return (truth=make_buffer(), configured=make_buffer())
end

function energy_buffers(ncol)
    make_buffer() = Dict(name => zeros(Float64, 12, ncol) for name in keys(ENERGY_VARIABLES))
    return (truth=make_buffer(), configured=make_buffer())
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
function accumulate_budget!(buffers, time, truth, configured)
    month_index = month(time)
    for name in keys(BUDGET_VARIABLES)
        @views buffers.truth[name][month_index, :] .+= truth[name]
        @views buffers.configured[name][month_index, :] .+= configured[name]
    end
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

function scatter_panel(x, y, label, period, max_points)
    xplot, yplot = selected_scatter_points(x, y, max_points)
    isempty(xplot) && return plot(title="$label — $period (no valid data)")
    rmse = sqrt(mean((yplot .- xplot) .^ 2))
    correlation = length(xplot) > 1 ? cor(xplot, yplot) : NaN
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
        label=false,
        xlabel="MAR ground truth",
        ylabel="Chion configured value",
        title="$label — $period\nR²=$(round(r_squared; digits=3)), r=$(round(correlation; digits=3)), RMSE=$(round(rmse; digits=3))",
        aspect_ratio=:equal,
        xlims=limits,
        ylims=limits,
    )
    plot!(panel, [limits[1], limits[2]], [limits[1], limits[2]]; color=:black, linewidth=1.2, linestyle=:dash, label=false)
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

function budget_figure(comparison, case)
    panels = Any[]
    for (name, label) in pairs(BUDGET_VARIABLES)
        push!(panels, scatter_panel(vec(comparison.truth[name]), vec(comparison.configured[name]), label, "monthly", CONFIG.scatter_max_points))
        push!(panels, scatter_panel(comparison.yearly_truth[name], comparison.yearly_configured[name], label, "yearly", CONFIG.scatter_max_points))
    end
    return plot(
        panels...;
        layout=(length(BUDGET_VARIABLES), 2),
        size=(1500, 440 * length(BUDGET_VARIABLES)),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

"""Plot area-integrated monthly budgets for the final diagnostic year in Gt."""
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
        values = hcat(mar, chion)
        labels = ["MAR" "Chion"]
        annotation = "MAR=$(round(mar_yearly; digits=1)) Gt\nChion=$(round(chion_yearly; digits=1)) Gt"
        if name == :smb_ice
            delta_firn = integrate_monthly(comparison.delta_firn)
            delta_firn_yearly = sum(delta_firn)
            chion_surface_smb = chion_yearly + delta_firn_yearly
            annotation = "MAR=$(round(mar_yearly; digits=1)) Gt\nChion=$(round(chion_yearly; digits=1)) Gt\nΔM firn=$(round(delta_firn_yearly; digits=1)) Gt\nChion + ΔM=$(round(chion_surface_smb; digits=1)) Gt"
            panel = plot(
                1:12,
                mar;
                label="MAR surface SMB",
                marker=:circle,
                markerstrokewidth=0,
                linewidth=2.8,
                xlabel="Month",
                ylabel="$label [Gt]",
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
            annotation_y = maximum(vcat(mar, chion, chion .+ delta_firn))
            annotate!(panel, 1.1, annotation_y, Plots.text(annotation, 7, :left, :top))
            push!(panels, panel)
            continue
        end
        annotation_y = maximum(values)
        push!(
            panels,
            plot(
                1:12,
                values;
                label=labels,
                marker=[:circle :square],
                markerstrokewidth=0,
                linewidth=2.8,
                xlabel="Month",
                ylabel="$label [Gt]",
                title=panel_label,
                titlelocation=:left,
                titlefont=Plots.font(11),
                xticks=1:12,
                framestyle=:box,
            ),
        )
        annotate!(panels[end], 1.1, annotation_y, Plots.text(annotation, 7, :left, :top))
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
    integrate_monthly(values) = vec(sum(values .* reshape(area_km2, 1, :); dims=2)) .* 1e-9
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
        annotate!(panel, 1.1, maximum(vcat(mar, chion)), Plots.text(annotation, 7, :left, :top))
        push!(panels, panel)
    end
    return plot(
        panels...;
        layout=(3, 3),
        size=(1800, 1500),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=8Plots.mm,
    )
end

function map_panel(grid, values, title; coastline_level)
    field = Chion.scatter_to_grid(Float64.(values), grid.js, grid.is, size(grid.mask))
    finite_values = field[isfinite.(field)]
    nonnegative = all(>=(zero(eltype(finite_values))), finite_values)
    colormap = nonnegative ? :viridis : :RdBu
    limits = if nonnegative
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
        title=title,
        aspect_ratio=:equal,
        c=colormap,
        clims=limits,
        colorbar=true,
        axis=false,
        grid=false,
        framestyle=:none,
    )
    contour!(panel, grid.x, grid.y, grid.mask; levels=[coastline_level], color=:black, linewidth=1.2, label=false)
    return panel
end

function end_state_figure(grid, state, annual, case)
    fields = (
        (state.thickness, "End-of-year snow thickness (m)"),
        (annual.smb_ice, "Annual ice SMB (mmWE)"),
        (annual.runoff, "Annual runoff (mmWE)"),
        (annual.melt, "Annual melt (mmWE)"),
        (annual.refreezing, "Annual refreezing (mmWE)"),
        (annual.sublimation, "Annual sublimation (mmWE)"),
        (state.albedo, "End-of-year surface albedo (-)"),
    )
    panels = [map_panel(grid, values, title; coastline_level=CONFIG.mask_threshold) for (values, title) in fields]
    return plot(
        panels...;
        layout=(3, 3),
        size=(1800, 1500),
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        top_margin=8Plots.mm,
    )
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
    panels = (
        map_panel(grid, mar_refreezing, "MAR annual refreezing (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, chion_refreezing, "Chion annual refreezing (mmWE)"; coastline_level=CONFIG.mask_threshold),
        map_panel(grid, chion_refreezing .- mar_refreezing, "Chion − MAR refreezing (mmWE)"; coastline_level=CONFIG.mask_threshold),
    )
    return plot(panels...; layout=(1, 3), size=(1800, 600), title="$(case.title): refreezing diagnosis")
end

function annual_budget_start(state)
    return NamedTuple(name => copy(getfield(state, name)) for name in CUMULATIVE_BUDGET_VARIABLES)
end

function annual_budget_change(state, start)
    return NamedTuple(name => getfield(state, name) .- getfield(start, name) for name in keys(start))
end

"""Total solid plus liquid firn mass per column (mmWE)."""
function firn_mass(state)
    return vec(sum(state.mass .+ state.mass_w; dims=1))
end

"""Diagnose the firn state immediately before the 1 May forcing interval."""
function pre_melt_firn_state(state)
    ncol = length(state.N)
    cold_content = zeros(Float64, ncol)
    pore_capacity = zeros(Float64, ncol)
    liquid_water = zeros(Float64, ncol)
    depth_to_dense_layer = fill(NaN, ncol)
    c = state.c
    dense_threshold = 810.0
    for column in 1:ncol
        depth = 0.0
        for layer in 1:state.N[column]
            mass = state.mass[layer, column]
            density = state.density[layer, column]
            if mass <= 0 || density <= 0
                continue
            end
            cold_content[column] += max(c.T0 - state.temperature[layer, column], 0.0) * c.ci * mass / c.Lm
            pore_capacity[column] += max(mass / density - mass / c.rho_i, 0.0) * c.rho_w
            liquid_water[column] += state.mass_w[layer, column]
            if isnan(depth_to_dense_layer[column]) && density >= dense_threshold
                depth_to_dense_layer[column] = depth
            end
            depth += mass / density
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
    forcing = forcing_for_case(base_forcing, reference, case.prescribed)
    model = BESSIModel(
        grid;
        Ntot=CONFIG.Ntot,
        albedo=case.prescribed.albedo ? :prescribed : :dynamic,
        seb_scheme=:bessi,
        densification=:bessi,
        fresh_snow_density=:constant,
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
        air_temperature = forcing.air_temperature[:, time_index]
        surface_temperature = state.Tsrf
        parameterized = (
            q_sw_net=(one(eltype(state.albedo)) .- state.albedo) .* forcing.shortwave_down[:, time_index],
            albedo=state.albedo,
            q_lw_down=c.σ * c.ϵ_air .* air_temperature .^ 4,
            q_sh=c.D_sh .* (air_temperature .- surface_temperature),
            q_lh=[Chion._resolved_turbulent_latent_heat_flux(
                c,
                surface_temperature[column],
                air_temperature[column],
                false,
                zero(eltype(surface_temperature)),
                forcing.has_relative_humidity[column, time_index],
                forcing.relative_humidity[column, time_index],
                forcing.air_pressure[column, time_index],
                forcing.wind_speed[column, time_index],
            ) for column in eachindex(surface_temperature)],
        )
        truth = (
            q_sw_net=reference.q_sw_net[:, time_index],
            albedo=reference.albedo[:, time_index],
            q_lw_down=forcing.q_lw_down[:, time_index],
            q_sh=forcing.q_sh[:, time_index],
            q_lh=forcing.q_lh[:, time_index],
        )
        configured = NamedTuple(name => case.prescribed[name] ? truth[name] : parameterized[name] for name in keys(COMPARISON_VARIABLES))
        accumulate_comparison!(buffers, forcing.time_values[time_index], forcing.dt_days[time_index], truth, configured)
        dt_seconds = forcing.dt_days[time_index] * c.seconds_per_day
        mar_energy_flux = (
            shortwave=reference.energy.shortwave[:, time_index],
            longwave=reference.energy.longwave[:, time_index],
            sensible_heat=forcing.q_sh[:, time_index],
            latent_heat=forcing.q_lh[:, time_index],
            rain_heat=forcing.rainfall_rate[:, time_index] .* c.cw .* (air_temperature .- c.T0),
        )
        chion_energy_flux = (
            shortwave=configured.q_sw_net,
            longwave=forcing.q_lw_down[:, time_index] .- c.σ * c.ϵ_snow .* surface_temperature .^ 4,
            sensible_heat=configured.q_sh,
            latent_heat=configured.q_lh,
            rain_heat=forcing.rainfall_rate[:, time_index] .* c.cw .* (air_temperature .- c.T0),
        )
        Chion.step_model!(model, diagnostic.now, integrator.model_runtime, runtime.step_fields, time_index)
        cumulative_budget = NamedTuple(
            name => getfield(runtime.state, name) .- getfield(previous_budget, name)
            for name in CUMULATIVE_BUDGET_VARIABLES
        )
        model_budget = (
            cumulative_budget...,
            surface_balance=(forcing.snowfall_rate[:, time_index] .+
                             forcing.rainfall_rate[:, time_index]) .* 
                            (forcing.dt_days[time_index] * model.c.seconds_per_day) .-
                            cumulative_budget.runoff .+
                            cumulative_budget.vapor_mass,
        )
        mar_budget = NamedTuple(
            name => reference.budgets[name][:, time_index] .* forcing.dt_days[time_index]
            for name in keys(BUDGET_VARIABLES)
        )
        accumulate_budget!(budgets, forcing.time_values[time_index], mar_budget, model_budget)
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

    for case in selected_cases
        @info "Running Greenland comparison" case=case.name spinup_years=CONFIG.spinup_years
        result = run_case(grid, base_forcing, reference, case)
        scatter_path = joinpath(CONFIG.output_dir, "$(case.name)_flux_scatter.pdf")
        budget_path = joinpath(CONFIG.output_dir, "$(case.name)_budget_scatter.pdf")
        integrated_budget_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_budgets.pdf")
        energy_path = joinpath(CONFIG.output_dir, "$(case.name)_monthly_integrated_energy.pdf")
        state_path = joinpath(CONFIG.output_dir, "$(case.name)_end_state.pdf")
        refreezing_path = joinpath(CONFIG.output_dir, "$(case.name)_refreezing_difference.pdf")
        refreezing_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_refreezing_by_elevation.csv")
        firn_state_path = joinpath(CONFIG.output_dir, "$(case.name)_pre_melt_firn_state.pdf")
        firn_summary_path = joinpath(CONFIG.output_dir, "$(case.name)_pre_melt_firn_by_elevation.csv")
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
        savefig(comparison_figure(result.comparison, case), scatter_path)
        savefig(budget_figure(result.budgets, case), budget_path)
        savefig(monthly_integrated_budget_figure(result.budgets, reference.area_km2, case), integrated_budget_path)
        savefig(monthly_integrated_energy_figure(result.energy, reference.area_km2, case), energy_path)
        savefig(end_state_figure(grid, result.state, result.annual, case), state_path)
        savefig(refreezing_difference_figure(grid, mar_refreezing, result.annual.refreezing, case), refreezing_path)
        savefig(pre_melt_firn_figure(grid, result.firn_pre_melt, case), firn_state_path)
        write_refreezing_elevation_summary(refreezing_summary_path, refreezing_rows)
        write_firn_state_summary(firn_summary_path, firn_rows)
        @info "Refreezing by elevation" rows=refreezing_rows
        @info "Pre-melt firn state by elevation" rows=firn_rows
        @info "Wrote comparison figures" scatter_path budget_path integrated_budget_path energy_path state_path refreezing_path refreezing_summary_path firn_state_path firn_summary_path
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
