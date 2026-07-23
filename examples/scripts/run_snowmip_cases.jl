#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using Dates
using NCDatasets
using Plots
using Printf
using Statistics

const DATA_DIR = joinpath(@__DIR__, "..", "..", "data", "ESM-SnowMIP_all")
const PLOT_DIR = joinpath(@__DIR__, "..", "plots")

arg(name, default) = begin
    prefix = "--$(name)="
    for value in ARGS
        startswith(value, prefix) && return split(value, "="; limit=2)[2]
    end
    default
end

has_flag(name) = any(==("--$(name)"), ARGS)

function matching_file(data_dir, prefix, site)
    matches = filter(name -> startswith(name, "$(prefix)_$(site)_") && endswith(name, ".nc"), readdir(data_dir))
    length(matches) == 1 || error("Expected one $(prefix)_$(site)_*.nc file in $(data_dir), found $(length(matches)).")
    joinpath(data_dir, only(matches))
end

function read_values(ds, name)
    haskey(ds, name) || error("Variable $(name) is missing from $(path(ds)).")
    Float64.(coalesce.(ds[name][:], NaN))
end

function read_met(path)
    NCDataset(path) do ds
        dates = DateTime.(ds["time"][:])
        return (
            dates=dates,
            tair=read_values(ds, "Tair"),
            snow=max.(read_values(ds, "Snowf"), 0.0),
            rain=max.(read_values(ds, "Rainf"), 0.0),
            sw=max.(read_values(ds, "SWdown"), 0.0),
            lw=read_values(ds, "LWdown"),
            wind=max.(read_values(ds, "Wind"), 0.0),
        )
    end
end

function first_available(ds, names)
    for name in names
        haskey(ds, name) && return name
    end
    nothing
end

function read_obs(path)
    NCDataset(path) do ds
        depth_name = first_available(ds, ["snd_auto", "snd_can_auto", "snd_gap_auto", "snd_gap1_auto", "snd_gap2_auto", "snd_man"])
        swe_name = first_available(ds, ["snw_auto", "snw_man"])
        albedo_name = haskey(ds, "albs") ? "albs" : nothing
        albedo_long_name = isnothing(albedo_name) ? "" : lowercase(String(ds[albedo_name].attrib["long_name"]))
        isnothing(depth_name) && error("No supported snow-depth variable in $(path).")
        dates = DateTime.(ds["time"][:])
        return (
            dates=dates,
            depth=read_values(ds, depth_name),
            swe=isnothing(swe_name) ? fill(NaN, length(dates)) : read_values(ds, swe_name),
            albedo=isnothing(albedo_name) ? fill(NaN, length(dates)) : read_values(ds, albedo_name),
            has_surface_albedo=!isnothing(albedo_name) && !occursin("above-canopy", albedo_long_name),
        )
    end
end

function align_observations(target_dates, obs)
    index = Dict(date => i for (i, date) in pairs(obs.dates))
    depth = fill(NaN, length(target_dates))
    swe = fill(NaN, length(target_dates))
    albedo = fill(NaN, length(target_dates))
    for (i, date) in pairs(target_dates)
        j = get(index, date, 0)
        j == 0 && continue
        depth[i] = obs.depth[j]
        swe[i] = obs.swe[j]
        albedo[i] = obs.albedo[j]
    end
    depth, swe, albedo
end

function sanitize(values, fallback; nonnegative=false)
    output = copy(values)
    previous = fallback
    for i in eachindex(output)
        valid = isfinite(output[i]) && (!nonnegative || output[i] >= 0)
        valid ? (previous = output[i]) : (output[i] = previous)
    end
    output
end

function forcing_from_met(met)
    n = length(met.dates)
    n > 1 || error("SnowMIP forcing needs at least two time records.")
    dt_days = [Dates.value(met.dates[min(i + 1, n)] - met.dates[i]) / 86_400_000 for i in 1:(n - 1)]
    push!(dt_days, dt_days[end])

    SnowpackForcing(
        dt_days=dt_days,
        ncol=1,
        air_temperature=sanitize(met.tair, 268.0),
        snowfall_rate=replace(met.snow, NaN => 0.0),
        rainfall_rate=replace(met.rain, NaN => 0.0),
        shortwave_down=replace(met.sw, NaN => 0.0),
        q_lw_down=replace(met.lw, NaN => 0.0),
        has_q_lw_down=isfinite.(met.lw),
        wind_speed=sanitize(met.wind, 5.0; nonnegative=true),
        time_values=met.dates,
    )
end

function snapshot(state)
    values = get_state(state)
    depth = max(values["total_thickness"], 0.0)
    swe = values["total_mass"]
    (depth=depth, swe=swe, density=depth > 0 ? swe / depth : NaN, albedo=values["albedo"])
end

function run_case(met; ntot=60)
    forcing = forcing_from_met(met)
    model = BESSIModel(SnowpackGrid(1); Ntot=ntot, densification=:htessel, fresh_snow_density=:htessel)
    simulation = Simulation(model; forcing, backend=:cpu, write_netcdf=false)
    workspace = Chion.ColumnarStepWorkspace(simulation.now)
    n = length(met.dates)
    depth, swe, density, albedo = fill(NaN, n), fill(NaN, n), fill(NaN, n), fill(NaN, n)

    initial = snapshot(simulation.now)
    depth[1], swe[1], density[1], albedo[1] = initial.depth, initial.swe, initial.density, initial.albedo
    for i in 2:n
        Chion.step!(simulation.now, forcing, i - 1, workspace)
        state = snapshot(simulation.now)
        depth[i], swe[i], density[i], albedo[i] = state.depth, state.swe, state.density, state.albedo
    end
    (depth=depth, swe=swe, density=density, albedo=albedo)
end

function rmse(observed, simulated)
    valid = isfinite.(observed) .& isfinite.(simulated)
    any(valid) ? sqrt(mean((simulated[valid] .- observed[valid]) .^ 2)) : NaN
end

function daily_means(dates, values)
    sums = Dict{Date,Float64}()
    counts = Dict{Date,Int}()
    for (datetime, value) in zip(dates, values)
        isfinite(value) || continue
        date = Date(datetime)
        sums[date] = get(sums, date, 0.0) + value
        counts[date] = get(counts, date, 0) + 1
    end
    Dict(date => total / counts[date] for (date, total) in sums)
end

function normalized_rmse(dates, observed, simulated)
    observed_daily = daily_means(dates, observed)
    simulated_daily = daily_means(dates, simulated)
    common_dates = sort!(collect(intersect(keys(observed_daily), keys(simulated_daily))))
    length(common_dates) > 1 || return NaN
    observed_values = [observed_daily[date] for date in common_dates]
    simulated_values = [simulated_daily[date] for date in common_dates]
    observed_std = std(observed_values; corrected=false)
    observed_std > 0 || return NaN
    sqrt(mean((simulated_values .- observed_values) .^ 2)) / observed_std
end

function trim_upper_percentile(values, probability=0.99)
    finite = filter(isfinite, values)
    isempty(finite) && return copy(values)
    cutoff = quantile(finite, probability)
    [isfinite(value) && value <= cutoff ? value : NaN for value in values]
end

function climatology_day(date)
    month(date) == 2 && day(date) == 29 && return 0
    reference = month(date) >= 10 ? Date(2000, month(date), day(date)) : Date(2001, month(date), day(date))
    Dates.value(reference - Date(2000, 10, 1)) + 1
end

water_year(date) = month(date) >= 10 ? year(date) + 1 : year(date)

function water_year_curves(dates, values)
    sums = Dict{Int,Vector{Float64}}()
    counts = Dict{Int,Vector{Int}}()
    for (datetime, value) in zip(dates, values)
        isfinite(value) || continue
        date = Date(datetime)
        day_index = climatology_day(date)
        day_index == 0 && continue
        wy = water_year(date)
        sum_curve = get!(sums, wy, zeros(365))
        count_curve = get!(counts, wy, zeros(Int, 365))
        sum_curve[day_index] += value
        count_curve[day_index] += 1
    end
    curves = Dict{Int,Vector{Float64}}()
    for wy in sort(collect(keys(sums)))
        curves[wy] = [counts[wy][day] > 0 ? sums[wy][day] / counts[wy][day] : NaN for day in 1:365]
    end
    curves
end

function curve_mean(curves)
    output = fill(NaN, 365)
    for day in 1:365
        day_values = [curve[day] for curve in Base.values(curves) if isfinite(curve[day])]
        isempty(day_values) || (output[day] = mean(day_values))
    end
    output
end

function add_panel_label!(panel, label)
    xlo, xhi = xlims(panel)
    ylo, yhi = ylims(panel)
    annotate!(panel, xlo - 0.08 * (xhi - xlo), yhi + 0.05 * (yhi - ylo), text(label, 16, :black, :left, :bottom))
end

function save_plot(path, site, met, obs_depth, obs_swe, sim)
    obs_density = [isfinite(swe) && isfinite(depth) && depth > 0 ? swe / depth : NaN for (swe, depth) in zip(obs_swe, obs_depth)]
    depth_rmse = rmse(obs_depth, sim.depth)
    swe_rmse = rmse(obs_swe, sim.swe)
    density_rmse = rmse(obs_density, sim.density)
    obs_density_plot = trim_upper_percentile(obs_density)
    sim_density_plot = trim_upper_percentile(sim.density)

    common = (framestyle=:box, legend=:topright, left_margin=12Plots.mm,
        top_margin=7Plots.mm, bottom_margin=5Plots.mm, guidefontsize=12,
        tickfontsize=10, titlefontsize=13, legendfontsize=10)
    p1 = scatter(met.dates, obs_depth; label="Observed", ms=1.5, markerstrokewidth=0, color=:steelblue,
        ylabel="Snow depth (m)", title=@sprintf("RMSE = %.3f m", depth_rmse), common...)
    plot!(p1, met.dates, sim.depth; label="Simulated", lw=1.5, color=:firebrick)
    p2 = scatter(met.dates, obs_swe; label="Observed", ms=2, markerstrokewidth=0, color=:steelblue,
        ylabel="SWE (kg m⁻²)", title=@sprintf("RMSE = %.1f kg m⁻²", swe_rmse), common...)
    plot!(p2, met.dates, sim.swe; label="Simulated", lw=1.5, color=:firebrick)
    p3 = scatter(met.dates, obs_density_plot; label="Observed", ms=2, markerstrokewidth=0, color=:steelblue,
        ylabel="Bulk density (kg m⁻³)", title=@sprintf("RMSE = %.1f kg m⁻³", density_rmse), common...)
    plot!(p3, met.dates, sim_density_plot; label="Simulated", lw=1.5, color=:firebrick)
    p4 = scatter(obs_depth, sim.depth; label=false, ms=2, markerstrokewidth=0, alpha=0.4, color=:firebrick,
        xlabel="Observed snow depth (m)", ylabel="Simulated snow depth (m)", common...)
    finite_depth = filter(isfinite, [obs_depth; sim.depth])
    if !isempty(finite_depth)
        limits = extrema(finite_depth)
        plot!(p4, collect(limits), collect(limits); label="1:1", color=:black, linestyle=:dash)
    end
    p5 = plot(met.dates, met.tair .- 273.15; label=false, color=:teal, lw=1,
        ylabel="Air temperature (°C)", common...)
    p6 = plot(met.dates, met.sw; label="SWdown", color=:purple, lw=1,
        ylabel="Radiation (W m⁻²)", common...)
    plot!(p6, met.dates, met.lw; label="LWdown", color=:dodgerblue, lw=1)

    obs_years = water_year_curves(met.dates, obs_depth)
    sim_years = water_year_curves(met.dates, sim.depth)
    p7 = plot(; xlabel="Month (October–September)", ylabel="Snow-depth climatology (m)",
        xticks=([1, 93, 183, 274, 336], ["Oct", "Jan", "Apr", "Jul", "Sep"]), common...)
    for curve in values(obs_years)
        plot!(p7, 1:365, curve; label=false, color=:steelblue, alpha=0.18, lw=0.7)
    end
    for curve in values(sim_years)
        plot!(p7, 1:365, curve; label=false, color=:firebrick, alpha=0.18, lw=0.7)
    end
    plot!(p7, 1:365, curve_mean(obs_years); label="Observed mean", color=:steelblue, lw=2.5)
    plot!(p7, 1:365, curve_mean(sim_years); label="Simulated mean", color=:firebrick, lw=2.5)

    p8 = plot(met.dates, 86_400 .* met.snow; label="Snowfall", color=:steelblue, lw=0.8,
        ylabel="Precipitation (mm d⁻¹)", common...)
    plot!(p8, met.dates, 86_400 .* met.rain; label="Rainfall", color=:darkorange, lw=0.8)

    panels = [p1, p2, p3, p4, p5, p6, p7, p8]
    for (index, panel) in enumerate(panels)
        add_panel_label!(panel, "($(Char('a' + index - 1)))")
    end

    figure = plot(panels...; layout=(4, 2), size=(1500, 1450),
        plot_title="ESM-SnowMIP $(uppercase(site)) — in situ forcing")
    savefig(figure, path)
end

function save_site_climatology_figure(path, cases, field, ylabel;
    ylims_value=nothing, show_observed=case -> true)
    panels = Any[]
    for (index, case) in enumerate(cases)
        observed_mean = getproperty(case.observed, field)
        simulated_mean = getproperty(case.simulated, field)
        panel = plot(;
            title=uppercase(case.site),
            ylabel=ylabel,
            xticks=([1, 93, 183, 274, 336], ["Oct", "Jan", "Apr", "Jul", "Sep"]),
            xlims=(1, 365),
            framestyle=:box,
            legend=index == 1 ? :topright : false,
            left_margin=12Plots.mm,
            top_margin=7Plots.mm,
            bottom_margin=5Plots.mm,
            guidefontsize=13,
            tickfontsize=11,
            titlefontsize=14,
            legendfontsize=11,
        )
        if show_observed(case)
            plot!(panel, 1:365, observed_mean;
                label=index == 1 ? "Observed" : false,
                color=:steelblue,
                lw=2.2,
            )
        end
        plot!(panel, 1:365, simulated_mean;
            label=index == 1 ? "Simulated" : false,
            color=:firebrick,
            lw=2.2,
        )
        isnothing(ylims_value) || ylims!(panel, ylims_value)
        add_panel_label!(panel, "($(Char('a' + index - 1)))")
        push!(panels, panel)
    end
    ncols = min(2, length(panels))
    nrows = cld(length(panels), ncols)
    figure = plot(panels...; layout=(nrows, ncols), size=(1500, 360 * nrows))
    savefig(figure, path)
    println("Wrote ", abspath(path))
end

function save_swe_nrmse_figure(cases, out_dir)
    sites = uppercase.([case.site for case in cases])
    values = [case.swe_nrmse for case in cases]
    # Approximate values digitized from the visible white markers in Fig. 4b
    # of Krinner et al. (2018), doi:10.5194/gmd-11-5027-2018. Markers that
    # overlap exactly in the source raster cannot be resolved independently.
    snowmip_models = Dict(
        "CDP" => [1.65, 1.60, 1.56, 1.41, 1.06, 0.94, 0.91, 0.88, 0.85,
                  0.80, 0.73, 0.66, 0.60, 0.54, 0.49, 0.44, 0.40],
        "OAS" => [1.93, 1.78, 1.74, 1.59, 1.56, 1.52, 1.47, 1.36, 1.25,
                  1.09, 1.04, 0.96, 0.92, 0.88, 0.83, 0.78],
        "OBS" => [2.03, 1.94, 1.90, 1.69, 1.65, 1.60, 1.55, 1.49, 1.39,
                  1.34, 1.28, 1.21, 1.16, 1.10, 1.04, 0.98, 0.92, 0.87,
                  0.82, 0.68],
        "OJP" => [1.72, 1.68, 1.47, 1.42, 1.38, 1.29, 1.23, 1.17, 1.11,
                  1.05, 0.99, 0.94, 0.90, 0.84, 0.68, 0.64, 0.60, 0.57],
        "RME" => [1.23, 1.18, 0.89, 0.73, 0.68, 0.63, 0.57, 0.52, 0.47,
                  0.42, 0.38, 0.34, 0.30],
        "SAP" => [2.04, 2.01, 1.85, 1.75, 1.40, 1.11, 1.01, 0.81, 0.78,
                  0.75, 0.71, 0.67, 0.64, 0.60, 0.57, 0.54],
        "SNB" => [2.06, 2.02, 1.61, 0.92, 0.88, 0.84, 0.80, 0.76, 0.72,
                  0.68, 0.64, 0.61],
        "SOD" => [1.31, 1.10, 0.69, 0.66, 0.63, 0.54, 0.49, 0.44, 0.40,
                  0.36],
        "SWA" => [2.36, 1.99, 1.21, 1.18, 0.85, 0.82, 0.78, 0.75, 0.71,
                  0.67, 0.63, 0.59, 0.56, 0.52, 0.49, 0.46, 0.42],
        "WFJ" => [1.58, 1.54, 1.49, 0.95, 0.91, 0.80, 0.69, 0.65, 0.60,
                  0.56, 0.53, 0.50, 0.47, 0.44, 0.40, 0.36],
    )
    x = collect(eachindex(sites))
    forest_sites = Set(["OAS", "OBS", "OJP"])
    open_x = Float64[]
    open_values = Float64[]
    forest_x = Float64[]
    forest_values = Float64[]
    for (site_index, site) in enumerate(sites)
        site_values = get(snowmip_models, site, Float64[])
        if site in forest_sites
            append!(forest_x, fill(site_index, length(site_values)))
            append!(forest_values, site_values)
        else
            append!(open_x, fill(site_index, length(site_values)))
            append!(open_values, site_values)
        end
    end
    panel = scatter(
        open_x,
        open_values;
        label="SnowMIP models: open sites (digitized)",
        marker=:circle,
        markersize=5,
        markercolor=:white,
        markerstrokecolor=:gray35,
        markerstrokewidth=1.3,
        xlabel="Site",
        ylabel="SWE normalised RMSE",
        framestyle=:box,
        ylims=(0, 2.5),
        size=(1200, 650),
        left_margin=14Plots.mm,
        bottom_margin=9Plots.mm,
        guidefontsize=14,
        tickfontsize=12,
        legendfontsize=11,
        legend=:topright,
        xticks=(x, sites),
    )
    scatter!(
        panel,
        forest_x,
        forest_values;
        label="SnowMIP models: forested sites (digitized)",
        marker=:utriangle,
        markersize=6,
        markercolor=:white,
        markerstrokecolor=:gray35,
        markerstrokewidth=1.3,
    )
    scatter!(
        panel,
        x .+ 0.10,
        values;
        label="Simulated",
        marker=:diamond,
        markersize=7,
        markercolor=:firebrick,
        markerstrokecolor=:firebrick,
    )
    hline!(panel, [1.0]; label=false, color=:gray45, linestyle=:dash, linewidth=1.5)
    output = joinpath(out_dir, "validation_esm_snowmip_swe_nrmse.pdf")
    savefig(panel, output)
    println("Wrote ", abspath(output))
end

function save_climatology_figures(cases, out_dir)
    save_site_climatology_figure(
        joinpath(out_dir, "validation_esm_snowmip_snow_depth_climatology.pdf"),
        cases,
        :depth,
        "Snow-depth climatology (m)",
    )
    save_site_climatology_figure(
        joinpath(out_dir, "validation_esm_snowmip_swe_climatology.pdf"),
        cases,
        :swe,
        "SWE climatology (kg m⁻²)",
    )
    save_site_climatology_figure(
        joinpath(out_dir, "validation_esm_snowmip_albedo_climatology.pdf"),
        cases,
        :albedo,
        "Albedo climatology (-)";
        ylims_value=(0, 1),
        show_observed=case -> case.has_surface_albedo,
    )
    save_swe_nrmse_figure(cases, out_dir)
end

function run_site(site, data_dir, out_dir, ntot; site_plot=true)
    met_path = matching_file(data_dir, "met_insitu", site)
    obs_path = matching_file(data_dir, "obs_insitu", site)
    met, obs = read_met(met_path), read_obs(obs_path)
    obs_depth, obs_swe, obs_albedo = align_observations(met.dates, obs)
    sim = run_case(met; ntot)
    if site_plot
        output = joinpath(out_dir, "validation_esm_$(site)_insitu_ldhtessel_evolution.pdf")
        save_plot(output, site, met, obs_depth, obs_swe, sim)
        println("Wrote ", abspath(output))
    end
    return (
        site=site,
        has_surface_albedo=obs.has_surface_albedo,
        swe_nrmse=normalized_rmse(met.dates, obs_swe, sim.swe),
        observed=(
            depth=curve_mean(water_year_curves(met.dates, obs_depth)),
            swe=curve_mean(water_year_curves(met.dates, obs_swe)),
            albedo=curve_mean(water_year_curves(met.dates, obs_albedo)),
        ),
        simulated=(
            depth=curve_mean(water_year_curves(met.dates, sim.depth)),
            swe=curve_mean(water_year_curves(met.dates, sim.swe)),
            albedo=curve_mean(water_year_curves(met.dates, sim.albedo)),
        ),
    )
end

function main()
    data_dir = abspath(arg("data-dir", DATA_DIR))
    out_dir = abspath(arg("out-dir", PLOT_DIR))
    ntot = parse(Int, arg("ntot", "60"))
    site_plot = !has_flag("comparison-only")
    requested = lowercase.(strip.(split(arg("sites", "all"), ",")))
    sites = requested == ["all"] ? sort([match(r"^met_insitu_([a-z0-9]+)_", name).captures[1]
        for name in readdir(data_dir) if occursin(r"^met_insitu_[a-z0-9]+_", name)]) : requested
    mkpath(out_dir)
    cases = NamedTuple[]
    for site in sites
        println("Running SnowMIP site ", uppercase(site), " …")
        push!(cases, run_site(site, data_dir, out_dir, ntot; site_plot))
    end
    save_climatology_figures(cases, out_dir)
end

main()
