#!/usr/bin/env julia

"""
Spin up a 15-layer Greenland BESSI state with fully prescribed MAR surface
forcing, continue that exact state through the available annual MAR files, and
plot ice-sheet-integrated annual mass-budget rates.  The transient includes
both a fully prescribed Chion configuration and a fully parameterized Chion
configuration (including prognostic, aging surface albedo).

The default forcing window is 1940--1980. Results are written as both CSV and
PDF so the plotted values remain easy to reuse.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using Dates: year
using NCDatasets
using Plots
using Statistics: mean

# CHION_DOMAIN selects the MAR domain. Antarctica has no daily climatology, so its
# spin-up repeats the first forcing year (1979); it also lacks AREA and RZ.
const DOMAIN = Symbol(get(ENV, "CHION_DOMAIN", "greenland"))
DOMAIN in (:greenland, :antarctica) || error("CHION_DOMAIN must be greenland or antarctica.")
const DOMAIN_SETUP = DOMAIN == :antarctica ? (
    x_name="X18_215", y_name="Y15_176", file_tag="27.5km", first_year=1979,
    forcing_dir="/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily_antarctica",
    spinup_file="/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily_antarctica/MARv3.14.3-27.5km-daily-ERA5-1979.nc",
    # Older Antarctic files lack AREA; all years share the grid, so take it from 2025.
    area_file="/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily_antarctica/MARv3.14.3-27.5km-daily-ERA5-2025.nc",
) : (
    x_name="x", y_name="y", file_tag="10km", first_year=1940,
    forcing_dir="/p/projects/ou/labs/ai/Nils/MAR_daily/mar_daily",
    spinup_file="/p/projects/ou/labs/ai/Nils/MAR3.14/MARv3.14.3-10km-daily-ERA5-1940-1980_daily_climatology.nc",
    area_file=nothing,
)

const CONFIG = (
    forcing_dir="/p/projects/ou/labs/ai/Nils/MAR3.14",
    spinup_file=get(ENV, "CHION_SPINUP_FILE", DOMAIN_SETUP.spinup_file),
    output_dir=get(ENV, "CHION_OUTPUT_DIR", joinpath(@__DIR__, "..", "plots", "gris_fully_prescribed_transient")),
    mask_threshold=50.0,
    wind_speed_m_s=5.0,
    spinup_years=parse(Int, get(ENV, "CHION_SPINUP_YEARS", "200")),
    backend=Symbol(get(ENV, "CHION_BACKEND", "gpu")),
    snow_layers=15,
    diurnal_substeps=parse(Int, get(ENV, "CHION_DIURNAL_SUBSTEPS", "8")),
)

const PLOT_VARIABLES = (
    surface_smb="Surface SMB",
    melt="Melt",
    runoff="Runoff",
    refreezing="Refreezing",
    vapor_mass="Vapor mass",
)

function mar_air_temperature_name(path)
    return NCDataset(path) do dataset
        haskey(dataset, "TTZ") && return "TTZ"
        haskey(dataset, "TT") && return "TT"
        error("Neither TTZ nor TT is available in $(basename(path)).")
    end
end

function load_mar(path)
    return load_forcing_file(
        path;
        x_name=DOMAIN_SETUP.x_name,
        y_name=DOMAIN_SETUP.y_name,
        time_name="TIME",
        air_temperature_name=mar_air_temperature_name(path),
        wind_speed_name=nothing,
        wind_default=CONFIG.wind_speed_m_s,
        mask_name="MSK",
        mask_threshold=CONFIG.mask_threshold,
    )
end

function read_mar_columns(path, name, grid, ntime)
    rows = grid.is .+ (grid.js .- 1) .* length(grid.x)
    return NCDataset(path) do dataset
        read_forcing_columns(dataset, name, rows, ntime, length(grid.y), length(grid.x);
            x_name=DOMAIN_SETUP.x_name, y_name=DOMAIN_SETUP.y_name)
    end
end

mar_has(path, name) = NCDataset(dataset -> haskey(dataset, name), path)

function read_mar_static_columns(path, name, grid)
    if name == "AREA" && !mar_has(path, name)
        isnothing(DOMAIN_SETUP.area_file) || return read_mar_static_columns(DOMAIN_SETUP.area_file, name, grid)
        error("No AREA in $(basename(path)).")
    end
    field = NCDataset(path) do dataset
        raw, dim_names = Chion._read_variable_data(dataset, name)
        Chion._as_y_x(raw, dim_names, length(grid.y), length(grid.x), name)
    end
    return [field[grid.js[column], grid.is[column]] for column in eachindex(grid.js)]
end

"""Attach every MAR surface-energy field and MAR albedo to the base forcing."""
function fully_prescribed_forcing(path, loaded)
    grid, forcing = loaded.grid, loaded.forcing
    ntime = length(forcing.time_values)
    albedo = clamp.(read_mar_columns(path, "AL2", grid, ntime), 0.0, 1.0)
    # Recent MAR files (e.g. Antarctica from 2025) omit SWU; derive it from AL2.
    shortwave_up = mar_has(path, "SWU") ? read_mar_columns(path, "SWU", grid, ntime) : forcing.shortwave_down .* albedo
    return SnowpackForcing(
        time_values=forcing.time_values,
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=forcing.wind_speed,
        q_sw_net=forcing.shortwave_down .- shortwave_up,
        has_q_sw_net=true,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=true,
        q_sh=forcing.q_sh,
        has_q_sh=true,
        q_lh=forcing.q_lh,
        has_q_lh=true,
        relative_humidity=forcing.relative_humidity,
        has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure,
        surface_height=forcing.surface_height,
        prescribed_albedo=albedo,
        has_prescribed_albedo=true,
        latitude_deg=forcing.latitude_deg,
    )
end

"""Forcing for the fully parameterized GrIS case used by the comparison suite.

Retain MAR meteorology and precipitation, but explicitly disable every MAR
surface-energy product and prescribed albedo.  The raw forcing file contains
LWD, SHF, and LHF; passing it through unchanged would silently make the run a
mixed-forcing experiment rather than `parameterized_constant_2k`.
"""
function parameterized_forcing(forcing; use_prescribed_albedo=false, wind_speed=forcing.wind_speed)
    return SnowpackForcing(
        time_values=forcing.time_values,
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=wind_speed,
        q_sw_net=forcing.q_sw_net,
        has_q_sw_net=false,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=false,
        q_sh=forcing.q_sh,
        has_q_sh=false,
        q_lh=forcing.q_lh,
        has_q_lh=false,
        relative_humidity=forcing.relative_humidity,
        has_relative_humidity=forcing.has_relative_humidity,
        air_pressure=forcing.air_pressure,
        surface_height=forcing.surface_height,
        prescribed_albedo=forcing.prescribed_albedo,
        has_prescribed_albedo=use_prescribed_albedo,
        latitude_deg=forcing.latitude_deg,
    )
end

"""Daily MAR wind speed at its lowest level (10 m; UVZ also has 50 and 100 m)."""
read_lowest_wind_level(path, grid, ntime) = read_mar_columns(path, "UVZ", grid, ntime)

"""Wind for the parameterized turbulent fluxes: CHION_WIND_SOURCE=mar uses MAR 10 m
wind where the file has UVZ, otherwise the constant default wind."""
function forcing_wind(path, loaded)
    get(ENV, "CHION_WIND_SOURCE", "constant") == "mar" || return loaded.forcing.wind_speed
    has_uvz = NCDataset(dataset -> haskey(dataset, "UVZ"), path)
    has_uvz || (@warn "No UVZ in $(basename(path)); using constant wind."; return loaded.forcing.wind_speed)
    return read_lowest_wind_level(path, loaded.grid, length(loaded.forcing.time_values))
end

function annual_files(directory)
    files = filter(
        path -> occursin(Regex("MARv3\\.14\\.3-$(DOMAIN_SETUP.file_tag)-daily-ERA5-\\d{4}\\.nc\$"), basename(path)),
        readdir(directory; join=true),
    )
    isempty(files) && error("No annual MAR forcing files found in $directory")
    return sort(files; by=path -> parse(Int, only(match(r"-(\d{4})\.nc$", basename(path)).captures)))
end

same_columns(a, b) = a.x == b.x && a.y == b.y && a.js == b.js && a.is == b.is

parse_float(x) = parse(Float64, x)
parse_profile(x) = (v = parse.(Float64, split(x, ',')); length(v) == 4 ||
    error("CHION_NEAR_SURFACE_LAYER_MAX_THICKNESSES_M must provide four comma-separated values."); Tuple(v))

# Optional model overrides from the environment. Unset variables keep the
# BESSIModel defaults, which are the calibrated GrIS setup.
const MODEL_ENV_OVERRIDES = (
    ("CHION_SEB_SCHEME", :seb_scheme, Symbol),
    ("CHION_LONGWAVE_SCHEME", :longwave_scheme, Symbol),
    ("CHION_EPSILON_AIR", :ϵ_air, parse_float),
    ("CHION_SEMIX_SENSIBLE_EXCHANGE_FACTOR", :semix_sensible_exchange_factor, parse_float),
    ("CHION_SEMIX_STABLE_COEFFICIENT", :semix_stable_coefficient, parse_float),
    ("CHION_ALPHA_DRY", :alpha_dry, parse_float),
    ("CHION_ALPHA_WET", :alpha_wet, parse_float),
    ("CHION_ALPHA_ICE", :alpha_ice, parse_float),
    ("CHION_MAX_LWC_ALBEDO", :max_lwc_albedo, parse_float),
    ("CHION_AGING_COLD_DAYS", :aging_cold_timescale_days, parse_float),
    ("CHION_AGING_MELTING_DAYS", :aging_melting_timescale_days, parse_float),
    ("CHION_ICE_SUBSTRATE_LAYERS", :ice_substrate_layers, x -> parse(Int, x)),
    ("CHION_ICE_SUBSTRATE_TOP_THICKNESS_M", :ice_substrate_top_thickness_m, parse_float),
    ("CHION_NEAR_SURFACE_LAYER_MAX_THICKNESSES_M", :near_surface_layer_max_thicknesses_m, parse_profile),
)

function model_env_overrides()
    pairs = [key => parser(ENV[name]) for (name, key, parser) in MODEL_ENV_OVERRIDES if !isempty(get(ENV, name, ""))]
    amplitude = get(ENV, "CHION_DIURNAL_TEMPERATURE_CONSTANT_AMPLITUDE_C", "")
    if !isempty(amplitude)
        push!(pairs, :diurnal_temperature_amplitude_c => parse_float(amplitude),
            :diurnal_temperature_amplitude_max_c => parse_float(amplitude))
    end
    return NamedTuple(pairs)
end

"""Chion with every MAR surface flux and albedo prescribed."""
build_fully_prescribed_model(grid) = BESSIModel(grid; Ntot=CONFIG.snow_layers, albedo=:prescribed,
    diurnal_shortwave_max_substeps=CONFIG.diurnal_substeps, model_env_overrides()...)

"""Chion with the parameterized SEB and dynamic albedo used in the GrIS comparison."""
build_parameterized_model(grid) = BESSIModel(grid; Ntot=CONFIG.snow_layers,
    albedo=Symbol(get(ENV, "CHION_ALBEDO_SCHEME", "dynamic")),
    diurnal_shortwave_max_substeps=CONFIG.diurnal_substeps, model_env_overrides()...)

"""Write the resolved surface and snowpack configuration of `model` as `name,value` lines."""
function write_model_configuration(io, model)
    p, c = model.parameters, model.parameters.c
    println(io, "seb_scheme,$(Chion._uses_semix_seb(c) ? "semix" : "bessi")")
    println(io, "longwave_scheme,$(c.longwave_scheme == Chion.LONGWAVE_CLOUD_PROXY ? "cloud_proxy" : "graybody")")
    println(io, "semix_sensible_exchange_factor,$(c.semix_sensible_exchange_factor)\nsemix_stable_coefficient,$(c.semix_stable_coefficient)")
    println(io, "alpha_dry,$(c.alpha_dry)\nalpha_wet,$(c.alpha_wet)\nalpha_ice,$(c.alpha_ice)")
    println(io, "wind_source,$(get(ENV, "CHION_WIND_SOURCE", "constant"))")
    println(io, "snow_layers,$(p.Ntot)\nnear_surface_layer_max_thicknesses_m,$(join(p.near_surface_layer_max_thicknesses_m, '|'))")
    println(io, "ice_substrate_layers,$(p.ice_substrate_layers)\nice_substrate_top_thickness_m,$(p.ice_substrate_top_thickness_m)")
    println(io, "diurnal_substeps,$(p.diurnal_shortwave_max_substeps)\ndiurnal_temperature_amplitude_c,$(p.diurnal_temperature_amplitude)")
end

function state_budget(state)
    return (
        smb_ice=Array(state.smb_ice),
        melt=Array(state.melt),
        runoff=Array(state.runoff),
        refreezing=Array(state.refreezing),
        vapor_mass=Array(state.vapor_mass),
    )
end

function integrated_gt(values, area_km2)
    return sum(values .* area_km2) * 1e-6
end

# MAR's annual totals never change, so they are computed once per file and grid
# and kept in a small CSV next to the MAR forcing.
const MAR_BUDGET_CACHE = get(ENV, "CHION_MAR_BUDGET_CACHE",
    joinpath(@__DIR__, "..", "..", "output", "mar_forcing", "mar_annual_budget_cache.csv"))
const MAR_BUDGET_TERMS = (:surface_smb, :melt, :runoff, :refreezing, :vapor_mass)

mar_budget_key(path, area_km2) = "$(basename(path))|$(round(Int, mtime(path)))|$(length(area_km2))|$(sum(area_km2))"

function cached_mar_budget(key)
    isfile(MAR_BUDGET_CACHE) || return nothing
    for line in eachline(MAR_BUDGET_CACHE)
        fields = split(line, ',')
        length(fields) == 1 + length(MAR_BUDGET_TERMS) && fields[1] == key || continue
        values = tryparse.(Float64, fields[2:end])
        any(isnothing, values) && continue
        return NamedTuple{MAR_BUDGET_TERMS}(Tuple(values))
    end
    return nothing
end

function mar_annual_budget(path, grid, forcing, area_km2)
    key = mar_budget_key(path, area_km2)
    cached = cached_mar_budget(key)
    isnothing(cached) || return cached
    ntime = length(forcing.time_values)
    annual(name) = vec(sum(
        read_mar_columns(path, name, grid, ntime) .* reshape(forcing.dt_days, 1, :);
        dims=2,
    ))
    budget = (
        surface_smb=integrated_gt(annual("SMB"), area_km2),
        melt=integrated_gt(annual("ME"), area_km2),
        runoff=integrated_gt(annual("RU"), area_km2),
        # MAR Antarctica has no refreezing (RZ) output.
        refreezing=mar_has(path, "RZ") ? integrated_gt(annual("RZ"), area_km2) : NaN,
        vapor_mass=integrated_gt(-annual("SU"), area_km2),
    )
    try
        mkpath(dirname(MAR_BUDGET_CACHE))
        open(io -> println(io, join((key, (repr(getproperty(budget, t)) for t in MAR_BUDGET_TERMS)...), ',')),
            MAR_BUDGET_CACHE, "a")
    catch err
        @warn "Could not write the MAR budget cache" MAR_BUDGET_CACHE exception=err
    end
    return budget
end

function chion_annual_budget(before, after, forcing, area_km2)
    change(name) = getfield(after, name) .- getfield(before, name)
    dt_seconds = forcing.dt_days .* 86_400.0
    accumulation = vec(sum(
        (forcing.snowfall_rate .+ forcing.rainfall_rate) .* reshape(dt_seconds, 1, :);
        dims=2,
    ))
    runoff = change(:runoff)
    vapor_mass = change(:vapor_mass)
    return (
        surface_smb=integrated_gt(accumulation .- runoff .+ vapor_mass, area_km2),
        melt=integrated_gt(change(:melt), area_km2),
        runoff=integrated_gt(runoff, area_km2),
        refreezing=integrated_gt(change(:refreezing), area_km2),
        vapor_mass=integrated_gt(vapor_mass, area_km2),
    )
end

function write_results(path, rows)
    open(path, "w") do io
        columns = String["year"]
        for name in keys(PLOT_VARIABLES)
            append!(columns, (
                "mar_$(name)_gt_per_year",
                "chion_prescribed_$(name)_gt_per_year",
                "chion_parameterized_$(name)_gt_per_year",
            ))
        end
        println(io, join(columns, ','))
        for row in rows
            values = Any[row.year]
            for name in keys(PLOT_VARIABLES)
                append!(values, (
                    getfield(row.mar, name),
                    getfield(row.chion_prescribed, name),
                    getfield(row.chion_parameterized, name),
                ))
            end
            println(io, join(values, ','))
        end
    end
    return path
end

"""Read a previously completed annual-rate CSV without rerunning the model."""
function read_results(path)
    lines = readlines(path)
    isempty(lines) && error("Annual-rate CSV is empty: $path")
    header = split(first(lines), ',')
    column = Dict(name => index for (index, name) in enumerate(header))
    value(fields, name) = parse(Float64, fields[column[name]])
    rows = NamedTuple[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        fields = split(line, ',')
        budget(source) = NamedTuple{keys(PLOT_VARIABLES)}(
            Tuple(value(fields, "$(source)_$(name)_gt_per_year") for name in keys(PLOT_VARIABLES)),
        )
        push!(rows, (
            year=parse(Int, fields[column["year"]]),
            mar=budget("mar"),
            chion_prescribed=budget(haskey(column, "chion_prescribed_surface_smb_gt_per_year") ? "chion_prescribed" : "chion"),
            chion_parameterized=haskey(column, "chion_parameterized_surface_smb_gt_per_year") ? budget("chion_parameterized") : nothing,
        ))
    end
    return rows
end

function plot_results(path, rows)
    years = getfield.(rows, :year)
    panels = map(enumerate(pairs(PLOT_VARIABLES))) do (index, (name, label))
        mar = [getfield(row.mar, name) for row in rows]
        prescribed = [getfield(row.chion_prescribed, name) for row in rows]
        parameterized_available = !isnothing(first(rows).chion_parameterized)
        parameterized = parameterized_available ? [getfield(row.chion_parameterized, name) for row in rows] : Float64[]
        rmse = sqrt(mean((prescribed .- mar) .^ 2))
        statistics = (
            "MAR mean=$(round(mean(mar); digits=1)) Gt yr⁻¹\n" *
            "Prescribed mean=$(round(mean(prescribed); digits=1)), RMSE=$(round(rmse; digits=1))"
        )
        if parameterized_available
            parameterized_rmse = sqrt(mean((parameterized .- mar) .^ 2))
            statistics *= "\nParameterized mean=$(round(mean(parameterized); digits=1)), RMSE=$(round(parameterized_rmse; digits=1))"
        end
        values = parameterized_available ? hcat(mar, prescribed, parameterized) : hcat(mar, prescribed)
        labels = parameterized_available ? ["MAR" "Chion: prescribed" "Chion: parameterized"] : ["MAR" "Chion: prescribed"]
        panel = plot(
            years,
            values;
            label=labels,
            linewidth=2,
            xlabel="Year",
            ylabel="$label [Gt yr⁻¹]",
            title="($(Char('a' + index - 1)))",
            titlelocation=:left,
            framestyle=:box,
        )
        xspan = maximum(years) - minimum(years)
        # MAR refreezing is missing (NaN) for most Antarctic years.
        finite_values = filter(isfinite, values)
        ymin, ymax = isempty(finite_values) ? (-1.0, 1.0) : extrema(finite_values)
        yspan = max(ymax - ymin, eps(Float64))
        # Allocate a dedicated annotation band above all values. This avoids
        # obscuring curves and leaves the legend in its own corner.
        ylims!(panel, ymin - 0.06yspan, ymax + 0.38yspan)
        annotate!(
            panel,
            minimum(years) + 0.03xspan,
            ymax + 0.34yspan,
            Plots.text(statistics, 8, :left, :top),
        )
        return panel
    end
    figure = plot(
        panels...;
        layout=(3, 2),
        size=(1500, 1500),
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
    )
    savefig(figure, path)
    return path
end

function main()
    CONFIG.spinup_years > 0 || error("CHION_SPINUP_YEARS must be positive.")
    mkpath(CONFIG.output_dir)

    spinup_loaded = load_mar(CONFIG.spinup_file)
    spinup_forcing = fully_prescribed_forcing(CONFIG.spinup_file, spinup_loaded)
    prescribed_model = build_fully_prescribed_model(spinup_loaded.grid)
    parameterized_model = build_parameterized_model(spinup_loaded.grid)
    prescribed_spinup = Simulation(
        prescribed_model;
        forcing=spinup_forcing,
        years=CONFIG.spinup_years,
        backend=CONFIG.backend,
        write_netcdf=false,
        compute_year_metrics=false,
        name="gris_fully_prescribed_spinup",
    )
    result = run!(prescribed_spinup)
    result.status == :complete || error("Spin-up ended with status $(result.status)")
    prescribed_state = prescribed_spinup.now
    parameterized_spinup = Simulation(
        parameterized_model;
        forcing=parameterized_forcing(spinup_loaded.forcing),
        years=CONFIG.spinup_years,
        backend=CONFIG.backend,
        write_netcdf=false,
        compute_year_metrics=false,
        name="gris_parameterized_spinup",
    )
    result = run!(parameterized_spinup)
    result.status == :complete || error("Parameterized spin-up ended with status $(result.status)")
    parameterized_state = parameterized_spinup.now
    @info "Spin-up complete" years=CONFIG.spinup_years

    files = annual_files(CONFIG.forcing_dir)
    area_km2 = read_mar_static_columns(first(files), "AREA", prescribed_model.grid)
    rows = NamedTuple[]
    for path in files
        loaded = load_mar(path)
        same_columns(prescribed_model.grid, loaded.grid) || error(
            "Selected MAR columns changed in $(basename(path)); cannot continue the spun-up state.",
        )
        forcing = fully_prescribed_forcing(path, loaded)
        prescribed_before = state_budget(prescribed_state)
        prescribed_transient = Simulation(
            prescribed_model;
            forcing,
            state=prescribed_state,
            years=1,
            backend=CONFIG.backend,
            write_netcdf=false,
            compute_year_metrics=false,
            name="gris_fully_prescribed_$(basename(path))",
        )
        result = run!(prescribed_transient)
        result.status == :complete || error("$(basename(path)) ended with status $(result.status)")
        prescribed_state = prescribed_transient.now
        prescribed_after = state_budget(prescribed_state)
        parameterized_before = state_budget(parameterized_state)
        parameterized_transient = Simulation(
            parameterized_model;
            forcing=parameterized_forcing(loaded.forcing),
            state=parameterized_state,
            years=1,
            backend=CONFIG.backend,
            write_netcdf=false,
            compute_year_metrics=false,
            name="gris_parameterized_$(basename(path))",
        )
        result = run!(parameterized_transient)
        result.status == :complete || error("Parameterized $(basename(path)) ended with status $(result.status)")
        parameterized_state = parameterized_transient.now
        parameterized_after = state_budget(parameterized_state)
        push!(rows, (
            year=year(first(forcing.time_values)),
            mar=mar_annual_budget(path, prescribed_model.grid, forcing, area_km2),
            chion_prescribed=chion_annual_budget(prescribed_before, prescribed_after, forcing, area_km2),
            chion_parameterized=chion_annual_budget(parameterized_before, parameterized_after, loaded.forcing, area_km2),
        ))
        @info "Completed transient year" year=rows[end].year
    end

    csv_path = write_results(joinpath(CONFIG.output_dir, "yearly_integrated_rates.csv"), rows)
    figure_path = plot_results(joinpath(CONFIG.output_dir, "yearly_integrated_rates.pdf"), rows)
    @info "Wrote fully prescribed transient diagnostics" csv_path figure_path
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
