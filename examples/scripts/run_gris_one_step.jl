#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using HDF5
import Plots
using Printf
using Statistics
using Base.Threads
using Chion

include(joinpath(@__DIR__, "..", "shared", "forcing_file_problem_loader.jl"))
using .ChionForcingFileProblemLoader

const SM = Chion.SnowpackModel
const DEFAULT_NC_PATH = begin
    candidate = "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc"
    isfile(candidate) ? candidate : ""
end
const DEFAULT_OUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_one_step")

plots_module() = Plots

function state_arrays(domain::SM.SnowpackDomain)
    return (
        N=domain.N,
        mass=domain.mass,
        mass_w=domain.mass_w,
        density=domain.density,
        temperature=domain.temperature,
        mass_base=domain.mass_base,
        smb_ice=domain.smb_ice,
        runoff=domain.runoff,
        Tsrf=domain.Tsrf,
        snow_cover=domain.snow_cover,
        albedo_dynamic=domain.albedo_dynamic,
    )
end

function print_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_one_step.jl [options]")
    println()
    println("Options:")
    println("  --nc=PATH                    Prepared HDF5/NetCDF forcing file")
    println("  --time-index=N               1-based time index to run (default: 1)")
    println("  --date=YYYYMMDDHH            Pick a time step by DATE value")
    println("  --out-dir=PATH               Output directory (default: examples/plots/gris_one_step)")
    println("  --ntot=N                     Chion maximum active layers (default: 80)")
    println("  --help                       Show this message")
end

function arg_value(args::Vector{String}, name::String, default::String)
    prefix = "--" * name * "="
    for arg in args
        startswith(arg, prefix) && return arg[(length(prefix) + 1):end]
    end
    return default
end

function has_flag(args::Vector{String}, name::String)
    return any(arg -> arg == "--" * name, args)
end

function parse_config(args::Vector{String})
    nc_path = arg_value(args, "nc", DEFAULT_NC_PATH)
    isempty(nc_path) && error("Pass --nc=PATH or place the prepared forcing file at $(DEFAULT_NC_PATH).")
    return (
        nc_path = nc_path,
        time_index = parse(Int, arg_value(args, "time-index", "1")),
        date_code = arg_value(args, "date", ""),
        out_dir = arg_value(args, "out-dir", DEFAULT_OUT_DIR),
        ntot = parse(Int, arg_value(args, "ntot", "80")),
    )
end

ensure_tools!() = nothing

function _summarize_column_state(
    n::Int,
    mass::AbstractVector{<:Real},
    mass_w::AbstractVector{<:Real},
    density::AbstractVector{<:Real},
    temperature::AbstractVector{<:Real},
    mass_base::Float64,
    smb_ice::Float64,
    runoff::Float64,
    T0::Float64,
)
    if n <= 0
        return (
            snow_mass = 0.0,
            liquid_mass = 0.0,
            wet_mass = 0.0,
            thickness = 0.0,
            bulk_density = NaN,
            base_mass = mass_base,
            smb_ice = smb_ice,
            runoff = runoff,
            snow_cover = 0.0,
            surface_temperature_c = NaN,
        )
    end

    snow_mass = 0.0
    liquid_mass = 0.0
    thickness = 0.0
    @inbounds for k in 1:n
        m = max(Float64(mass[k]), 0.0)
        mw = max(Float64(mass_w[k]), 0.0)
        rho = Float64(density[k])
        snow_mass += m
        liquid_mass += mw
        if m > 0.0 && rho > 0.0
            thickness += m / rho
        end
    end

    wet_mass = snow_mass + liquid_mass
    bulk_density = thickness > 0.0 ? snow_mass / thickness : NaN
    snow_cover = if wet_mass <= 0.0 || !isfinite(bulk_density) || bulk_density <= 0.0
        0.0
    else
        min(1.0, (wet_mass / bulk_density) / 0.1)
    end

    return (
        snow_mass = snow_mass,
        liquid_mass = liquid_mass,
        wet_mass = wet_mass,
        thickness = thickness,
        bulk_density = bulk_density,
        base_mass = mass_base,
        smb_ice = smb_ice,
        runoff = runoff,
        snow_cover = snow_cover,
        surface_temperature_c = Float64(temperature[1]) - T0,
    )
end

function summarize_column(domain::SM.SnowpackDomain, idx::Int)
    n = domain.N[idx]
    return _summarize_column_state(
        n,
        @view(domain.mass[:, idx]),
        @view(domain.mass_w[:, idx]),
        @view(domain.density[:, idx]),
        @view(domain.temperature[:, idx]),
        Float64(domain.mass_base[idx]),
        Float64(domain.smb_ice[idx]),
        Float64(domain.runoff[idx]),
        Float64(domain.c.T0),
    )
end

function masked_field(field::AbstractMatrix{<:Real}, valid_mask::BitMatrix)
    out = fill(NaN, size(field))
    @inbounds for I in eachindex(field)
        if valid_mask[I] && isfinite(field[I])
            out[I] = Float64(field[I])
        end
    end
    return out
end

function finite_extrema(field::AbstractMatrix{<:Real})
    vals = Float64[]
    sizehint!(vals, length(field))
    for x in field
        if isfinite(x)
            push!(vals, x)
        end
    end
    isempty(vals) && return (-1.0, 1.0)
    vmin = minimum(vals)
    vmax = maximum(vals)
    if vmin == vmax
        return (vmin - 1.0, vmax + 1.0)
    end
    return (vmin, vmax)
end

function symmetric_clims(field::AbstractMatrix{<:Real})
    vals = Float64[]
    sizehint!(vals, length(field))
    for x in field
        if isfinite(x)
            push!(vals, abs(Float64(x)))
        end
    end
    isempty(vals) && return (-1.0, 1.0)
    vmax = maximum(vals)
    vmax = vmax > 0.0 ? vmax : 1.0
    return (-vmax, vmax)
end

function heatmap_panel(x, y, field; title::String, unit::String="", clim=nothing, color=:viridis)
    P = plots_module()
    plot_clim = isnothing(clim) ? finite_extrema(field) : clim
    return P.heatmap(
        x,
        y,
        field;
        title=title,
        xlabel="x (km, EPSG:3413)",
        ylabel="y (km, EPSG:3413)",
        aspect_ratio=:equal,
        color=color,
        clims=plot_clim,
        colorbar_title=unit,
        framestyle=:box,
    )
end

function domain_mean(field::AbstractMatrix{<:Real})
    total = 0.0
    n = 0
    for x in field
        if isfinite(x)
            total += Float64(x)
            n += 1
        end
    end
    return n == 0 ? NaN : total / n
end

function write_summary(
    out_path::AbstractString,
    config,
    date_label::AbstractString,
    dt_days::Float64,
    valid_mask::BitMatrix,
    forcing::NamedTuple,
    response::NamedTuple,
)
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, "Chion GrIS one-step run")
        println(io, "Forcing file  : ", abspath(config.nc_path))
        println(io, "Date          : ", date_label)
        println(io, @sprintf("dt (days)     : %.3f", dt_days))
        println(io, "GrIS cells    : ", count(valid_mask), " / ", length(valid_mask))
        println(io)
        println(io, "Domain means over valid GrIS cells")
        println(io, @sprintf("T2m (C)               : %.3f", domain_mean(forcing.t2m_c)))
        println(io, @sprintf("Snowfall (mmWE/day)   : %.3f", domain_mean(forcing.snowfall_mm_day)))
        println(io, @sprintf("Rainfall (mmWE/day)   : %.3f", domain_mean(forcing.rainfall_mm_day)))
        println(io, @sprintf("SWD (W m-2)           : %.3f", domain_mean(forcing.swd)))
        println(io, @sprintf("LWD (W m-2)           : %.3f", domain_mean(forcing.lwd)))
        println(io, @sprintf("SHF+LHF used (W m-2)  : %.3f", domain_mean(forcing.turbulent_heat)))
        println(io, @sprintf("Initial snow h (m)    : %.3f", domain_mean(response.initial_thickness)))
        println(io, @sprintf("Snow h change (m)     : %.6f", domain_mean(response.delta_thickness)))
        println(io, @sprintf("SWE change (mmWE)     : %.6f", domain_mean(response.delta_wet_mass)))
        println(io, @sprintf("Bulk density (kg m-3) : %.3f", domain_mean(response.bulk_density)))
        println(io, @sprintf("Liquid water (mmWE)   : %.3f", domain_mean(response.liquid_water)))
        println(io, @sprintf("Runoff (mmWE)         : %.3f", domain_mean(response.runoff)))
    end
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_help()
        return
    end

    ensure_tools!()
    config = parse_config(args)
    shapes = read_dataset_shapes(config.nc_path)

    date_codes, time_values = read_forcing_times(config.nc_path, shapes)
    time_index = choose_time_index(config, date_codes)
    selected_time = time_values[time_index]
    dt_days = infer_dt_days(time_values, time_index)

    x = vec(read_hdf5_full(config.nc_path, "x", shapes))
    y = vec(read_hdf5_full(config.nc_path, "y", shapes))
    outlay_bounds = read_hdf5_full(config.nc_path, "OUTLAY_bnds", shapes)
    mask = read_hdf5_full(config.nc_path, "MSK", shapes)

    tt_c = read_timeslice_2d(config.nc_path, "TT", time_index, shapes)
    sf_mm_day = read_timeslice_2d(config.nc_path, "SF", time_index, shapes)
    rf_mm_day = read_timeslice_2d(config.nc_path, "RF", time_index, shapes)
    swd = read_timeslice_2d(config.nc_path, "SWD", time_index, shapes)
    lwd = read_timeslice_2d(config.nc_path, "LWD", time_index, shapes)
    shf = read_timeslice_2d(config.nc_path, "SHF", time_index, shapes)
    lhf = read_timeslice_2d(config.nc_path, "LHF", time_index, shapes)
    zn3 = read_timeslice_2d(config.nc_path, "ZN3", time_index, shapes)
    ro1 = read_timeslice_3d(config.nc_path, "RO1", time_index, shapes)
    ti1 = read_timeslice_3d(config.nc_path, "TI1", time_index, shapes)
    wa1 = read_timeslice_3d(config.nc_path, "WA1", time_index, shapes)

    ny, nx = size(mask)
    valid_mask = falses(ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        valid_mask[j, i] =
            isfinite(tt_c[j, i]) &&
            isfinite(sf_mm_day[j, i]) &&
            isfinite(rf_mm_day[j, i]) &&
            isfinite(swd[j, i])
    end

    initial_thickness = fill(NaN, ny, nx)
    delta_thickness = fill(NaN, ny, nx)
    delta_wet_mass = fill(NaN, ny, nx)
    bulk_density = fill(NaN, ny, nx)
    liquid_water = fill(NaN, ny, nx)
    runoff = fill(NaN, ny, nx)

    valid_indices = findall(valid_mask)
    nvalid = length(valid_indices)
    domain = SM.SnowpackDomain(ncol=nvalid, Ntot=config.ntot)
    workspaces = SM.threaded_workspaces(domain)
    initial_wet_mass = fill(NaN, ny, nx)
    T2m_vec = Vector{Float64}(undef, nvalid)
    P_snow_vec = Vector{Float64}(undef, nvalid)
    P_rain_vec = Vector{Float64}(undef, nvalid)
    S_boa_vec = Vector{Float64}(undef, nvalid)
    q_lw_vec = fill(NaN, nvalid)
    q_sh_vec = fill(NaN, nvalid)
    q_lh_vec = fill(NaN, nvalid)
    @threads :static for idx in eachindex(valid_indices)
        I = valid_indices[idx]
        j, i = Tuple(I)

        populate_domain_column_from_forcing_file!(
            domain,
            idx,
            Float64(zn3[j, i]),
            @view(ro1[:, j, i]),
            @view(ti1[:, j, i]),
            @view(wa1[:, j, i]),
            outlay_bounds,
        )

        initial = summarize_column(domain, idx)
        initial_thickness[j, i] = initial.thickness
        initial_wet_mass[j, i] = initial.wet_mass
        T2m_vec[idx] = valid_or(-15.0, Float64(tt_c[j, i])) + domain.c.T0
        P_snow_vec[idx] = mmwe_day_to_kgm2s(Float64(sf_mm_day[j, i]))
        P_rain_vec[idx] = mmwe_day_to_kgm2s(Float64(rf_mm_day[j, i]))
        S_boa_vec[idx] = valid_or(0.0, Float64(swd[j, i]))
        q_lw_vec[idx] = isfinite(lwd[j, i]) ? Float64(lwd[j, i]) : NaN
        q_sh_vec[idx] = isfinite(shf[j, i]) ? Float64(shf[j, i]) : NaN
        q_lh_vec[idx] = isfinite(lhf[j, i]) ? Float64(lhf[j, i]) : NaN
    end

    @threads :static for idx in eachindex(valid_indices)
        q_lw_down = isfinite(q_lw_vec[idx]) ? q_lw_vec[idx] : nothing
        q_sh_now = isfinite(q_sh_vec[idx]) ? q_sh_vec[idx] : nothing
        q_lh_now = isfinite(q_lh_vec[idx]) ? q_lh_vec[idx] : nothing
        SM.step!(
            domain,
            idx,
            T2m_vec[idx],
            P_snow_vec[idx] + P_rain_vec[idx],
            dt_days;
            workspace=workspaces[threadid()],
            snowfall_rate=P_snow_vec[idx],
            rainfall_rate=P_rain_vec[idx],
            shortwave_down=S_boa_vec[idx],
            wind_speed=10.0,
            q_lw_down=q_lw_down,
            q_sh=q_sh_now,
            q_lh=q_lh_now,
        )
    end

    @threads :static for idx in eachindex(valid_indices)
        I = valid_indices[idx]
        j, i = Tuple(I)
        final = summarize_column(domain, idx)
        delta_thickness[j, i] = final.thickness - initial_thickness[j, i]
        delta_wet_mass[j, i] = final.wet_mass - initial_wet_mass[j, i]
        bulk_density[j, i] = final.bulk_density
        liquid_water[j, i] = final.liquid_mass
        runoff[j, i] = final.runoff
    end

    forcing = (
        t2m_c = masked_field(tt_c, valid_mask),
        snowfall_mm_day = masked_field(sf_mm_day, valid_mask),
        rainfall_mm_day = masked_field(rf_mm_day, valid_mask),
        swd = masked_field(swd, valid_mask),
        lwd = masked_field(lwd, valid_mask),
        turbulent_heat = masked_field(shf .+ lhf, valid_mask),
    )
    response = (
        initial_thickness = masked_field(initial_thickness, valid_mask),
        delta_thickness = masked_field(delta_thickness, valid_mask),
        delta_wet_mass = masked_field(delta_wet_mass, valid_mask),
        bulk_density = masked_field(bulk_density, valid_mask),
        liquid_water = masked_field(liquid_water, valid_mask),
        runoff = masked_field(runoff, valid_mask),
    )

    mkpath(config.out_dir)
    date_label = Dates.format(selected_time, dateformat"yyyy-mm-dd HH:MM")
    P = plots_module()
    forcing_plot = P.plot(
        heatmap_panel(x, y, forcing.t2m_c; title="Air temperature", unit="C", color=:thermal),
        heatmap_panel(x, y, forcing.snowfall_mm_day; title="Snowfall", unit="mmWE/day", color=:ice),
        heatmap_panel(x, y, forcing.rainfall_mm_day; title="Rainfall", unit="mmWE/day", color=:blues),
        heatmap_panel(x, y, forcing.swd; title="SW down", unit="W/m^2", color=:solar),
        heatmap_panel(x, y, forcing.lwd; title="LW down", unit="W/m^2", color=:matter),
        heatmap_panel(x, y, forcing.turbulent_heat; title="SHF + LHF used", unit="W/m^2", clim=symmetric_clims(forcing.turbulent_heat), color=P.cgrad([:navy, :white, :firebrick])),
        layout=(2, 3),
        size=(1800, 1100),
        plot_title="Chion one-step forcing-file fields, $date_label",
    )

    response_plot = P.plot(
        heatmap_panel(x, y, response.initial_thickness; title="Initial snow height", unit="m", color=:ice),
        heatmap_panel(x, y, response.delta_thickness; title="Snow height change", unit="m", clim=symmetric_clims(response.delta_thickness), color=P.cgrad([:navy, :white, :firebrick])),
        heatmap_panel(x, y, response.delta_wet_mass; title="SWE change", unit="mmWE", clim=symmetric_clims(response.delta_wet_mass), color=P.cgrad([:navy, :white, :firebrick])),
        heatmap_panel(x, y, response.bulk_density; title="Bulk density", unit="kg/m^3", color=:dense),
        heatmap_panel(x, y, response.liquid_water; title="Liquid water", unit="mmWE", color=:amp),
        heatmap_panel(x, y, response.runoff; title="Runoff", unit="mmWE", color=:rainbow),
        layout=(2, 3),
        size=(1800, 1100),
        plot_title="Chion one-step response over the GrIS, $date_label",
    )

    forcing_plot_path = joinpath(config.out_dir, "gris_one_step_forcing.png")
    response_plot_path = joinpath(config.out_dir, "gris_one_step_response.png")
    summary_path = joinpath(config.out_dir, "gris_one_step_summary.txt")

    P.savefig(forcing_plot, forcing_plot_path)
    P.savefig(response_plot, response_plot_path)
    write_summary(summary_path, config, date_label, dt_days, valid_mask, forcing, response)

    println("GrIS one-step run complete.")
    println("Forcing file   : $(abspath(config.nc_path))")
    println("Time index     : $(time_index)")
    println("Date           : $(date_label)")
    println(@sprintf("dt (days)      : %.3f", dt_days))
    println("GrIS cells     : $(count(valid_mask)) / $(length(valid_mask))")
    println("Threads        : $(nthreads())")
    println("Forcing plot   : $(abspath(forcing_plot_path))")
    println("Response plot  : $(abspath(response_plot_path))")
    println("Summary        : $(abspath(summary_path))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
