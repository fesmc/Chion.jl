#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using Dates
using NCDatasets
include(joinpath(@__DIR__, "..", "shared", "script_helpers.jl"))
using .ChionExampleScriptHelpers

const DEFAULT_GRIS_MONTHLY_PDD_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_monthly_pdd_simulation")

function default_gris_monthly_pdd_forcing_file_path()
    for candidate in (
        get(ENV, "FORCING_PATH", ""),
        "/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc",
    )
        isfile(candidate) && return candidate
    end
    return ""
end

function apply_gris_mask(loaded, forcing_file::AbstractString; threshold::Float64=50.0)
    mask = NCDataset(forcing_file) do ds
        Float64.(ds["MSK"][:])
    end

    ny = length(loaded.grid.y)
    nx = length(loaded.grid.x)
    if length(mask) == nx * ny && ndims(mask) == 1
        mask = reshape(mask, nx, ny)
    end
    if size(mask) == (nx, ny)
        mask = permutedims(mask, (2, 1))
    end
    size(mask) == (ny, nx) || error("Unexpected MSK shape $(size(mask)); expected ($ny, $nx).")

    rows = Int[]
    js = Int[]
    is = Int[]
    for j in 1:ny, i in 1:nx
        row = i + (j - 1) * nx
        if isfinite(mask[j, i]) &&
           mask[j, i] >= threshold &&
           isfinite(loaded.forcing.air_temperature[row, 1])
            push!(rows, row)
            push!(js, j)
            push!(is, i)
        end
    end

    f = loaded.forcing
    forcing = SnowpackForcing(
        time_values=f.time_values,
        dt_days=f.dt_days,
        air_temperature=f.air_temperature[rows, :],
        snowfall_rate=f.snowfall_rate[rows, :],
        rainfall_rate=f.rainfall_rate[rows, :],
        shortwave_down=f.shortwave_down[rows, :],
        wind_speed=f.wind_speed[rows, :],
        q_lw_down=f.q_lw_down[rows, :],
        has_q_lw_down=f.has_q_lw_down[rows, :],
        q_sh=f.q_sh[rows, :],
        has_q_sh=f.has_q_sh[rows, :],
        q_lh=f.q_lh[rows, :],
        has_q_lh=f.has_q_lh[rows, :],
    )
    grid = SnowpackGrid(length(rows);
        x=loaded.grid.x,
        y=loaded.grid.y,
        js=js,
        is=is,
        mask=mask,
    )

    println("Applied GrIS mask: kept $(length(rows)) / $(nx * ny) columns")
    return (grid=grid, forcing=forcing)
end

function month_groups(time_values::Vector{DateTime})
    groups = Vector{UnitRange{Int}}()
    isempty(time_values) && return groups
    start = 1
    current = (year(time_values[1]), month(time_values[1]))
    for idx in 2:length(time_values)
        key = (year(time_values[idx]), month(time_values[idx]))
        if key != current
            push!(groups, start:idx-1)
            start = idx
            current = key
        end
    end
    push!(groups, start:length(time_values))
    return groups
end

function weighted_monthly_mean(field, dt_days, groups)
    ncol = size(field, 1)
    out = Matrix{Float64}(undef, ncol, length(groups))
    @views for (month_idx, indices) in pairs(groups)
        total_days = sum(dt_days[indices])
        out[:, month_idx] .= 0.0
        for time_idx in indices
            out[:, month_idx] .+= field[:, time_idx] .* dt_days[time_idx]
        end
        out[:, month_idx] ./= total_days
    end
    return out
end

function monthly_mean_optional_flux(field, has_field, dt_days, groups)
    ncol = size(field, 1)
    values = Matrix{Float64}(undef, ncol, length(groups))
    present = Matrix{Bool}(undef, ncol, length(groups))
    weights = Vector{Float64}(undef, ncol)
    @views for (month_idx, indices) in pairs(groups)
        values[:, month_idx] .= 0.0
        fill!(weights, 0.0)
        for time_idx in indices
            mask = has_field[:, time_idx]
            values[:, month_idx] .+= field[:, time_idx] .* mask .* dt_days[time_idx]
            weights .+= mask .* dt_days[time_idx]
        end
        present[:, month_idx] .= weights .> 0.0
        for col in 1:ncol
            values[col, month_idx] = present[col, month_idx] ? values[col, month_idx] / weights[col] : 0.0
        end
    end
    return values, present
end

function monthly_mean_precipitation_rate(rate, dt_days, groups)
    ncol = size(rate, 1)
    out = Matrix{Float64}(undef, ncol, length(groups))
    @views for (month_idx, indices) in pairs(groups)
        total_days = sum(dt_days[indices])
        out[:, month_idx] .= 0.0
        for time_idx in indices
            out[:, month_idx] .+= max.(rate[:, time_idx], 0.0) .* dt_days[time_idx]
        end
        out[:, month_idx] ./= total_days
    end
    return out
end

function aggregate_forcing_to_monthly(forcing::SnowpackForcing)
    groups = month_groups(forcing.time_values)
    isempty(groups) && error("Cannot aggregate empty forcing to monthly forcing.")
    monthly_time_values = DateTime[
        DateTime(year(forcing.time_values[first(indices)]), month(forcing.time_values[first(indices)]), 15, 12)
        for indices in groups
    ]
    monthly_dt_days = Float64[sum(forcing.dt_days[indices]) for indices in groups]
    q_lw_down, has_q_lw_down = monthly_mean_optional_flux(forcing.q_lw_down, forcing.has_q_lw_down, forcing.dt_days, groups)
    q_sh, has_q_sh = monthly_mean_optional_flux(forcing.q_sh, forcing.has_q_sh, forcing.dt_days, groups)
    q_lh, has_q_lh = monthly_mean_optional_flux(forcing.q_lh, forcing.has_q_lh, forcing.dt_days, groups)

    return SnowpackForcing(
        time_values=monthly_time_values,
        dt_days=monthly_dt_days,
        air_temperature=weighted_monthly_mean(forcing.air_temperature, forcing.dt_days, groups),
        snowfall_rate=monthly_mean_precipitation_rate(forcing.snowfall_rate, forcing.dt_days, groups),
        rainfall_rate=monthly_mean_precipitation_rate(forcing.rainfall_rate, forcing.dt_days, groups),
        shortwave_down=weighted_monthly_mean(forcing.shortwave_down, forcing.dt_days, groups),
        wind_speed=weighted_monthly_mean(forcing.wind_speed, forcing.dt_days, groups),
        q_lw_down=q_lw_down,
        has_q_lw_down=has_q_lw_down,
        q_sh=q_sh,
        has_q_sh=has_q_sh,
        q_lh=q_lh,
        has_q_lh=has_q_lh,
    )
end

function print_gris_monthly_pdd_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_monthly_pdd_case.jl [options]")
    println()
    println("Core options:")
    println("  --forcing-file=PATH              Prepared daily/subdaily NetCDF forcing file")
    println("  --mask-threshold=VALUE           Apply MAR-style MSK threshold before running (default: 50)")
    println("  --output-dir=PATH                Output directory")
    println("  --netcdf-path=PATH               Output NetCDF path")
    println("  --no-output                      Skip summary/CSV file output")
    println("  --no-nc                          Skip NetCDF output")
    println("  --netcdf-vars=SPEC               NetCDF variables to write (default: final,history,step)")
    println("  --years=N                        Number of repeated forcing years (default: 10)")
    println("  --backend=threads|cpu|gpu        Execution backend (default: threads)")
    println()
    println("PDD options:")
    println("  --pdd-ddf-snow=VALUE             Snow degree-day factor in mmWE d-1 C-1 (default: 3)")
    println("  --pdd-ddf-ice=VALUE              Ice degree-day factor in mmWE d-1 C-1 (default: 8)")
    println("  --pdd-refreezing-fraction=X      Refreezing fraction (default: 0.6)")
    println("  --pdd-temperature-sigma=VALUE    PISM-style monthly temperature variability in degC (default: 5)")
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_gris_monthly_pdd_help()
        return
    end

    default_forcing_file = default_gris_monthly_pdd_forcing_file_path()
    forcing_file = arg_value(args, "forcing-file", default_forcing_file)
    isempty(forcing_file) && error("Pass --forcing-file=PATH or set FORCING_PATH to the prepared forcing file.")

    loaded = load_forcing_file(forcing_file)
    loaded = apply_gris_mask(
        loaded,
        forcing_file;
        threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
    )
    monthly_forcing = aggregate_forcing_to_monthly(loaded.forcing)
    println("Aggregated forcing: $(length(loaded.forcing.time_values)) source steps -> $(length(monthly_forcing.time_values)) monthly steps")
    println("Monthly dt range: $(minimum(monthly_forcing.dt_days)) to $(maximum(monthly_forcing.dt_days)) days")

    model = PDDModel(
        loaded.grid;
        ddf_snow=parse(Float64, arg_value(args, "pdd-ddf-snow", env_value("PDD_DDF_SNOW", "8.0"))),
        ddf_ice=parse(Float64, arg_value(args, "pdd-ddf-ice", env_value("PDD_DDF_ICE", "16.0"))),
        refreezing_fraction=parse(Float64, arg_value(args, "pdd-refreezing-fraction", env_value("PDD_REFREEZING_FRACTION", "0.6"))),
        monthly_method=StochasticMonthlyPDD(
            temperature_sigma=parse(Float64, arg_value(args, "pdd-temperature-sigma", env_value("PDD_TEMPERATURE_SIGMA", "5.0"))),
        ),
    )

    netcdf_vars = has_flag(args, "no-nc") ? Symbol[] : arg_value(args, "netcdf-vars", "final,history,step")
    println("Model: pdd monthly")
    println("NetCDF variables: ", isempty(netcdf_vars) ? "(none)" : netcdf_vars)

    sim = Simulation(
        model;
        forcing=monthly_forcing,
        output=OutputOptions(
            save=netcdf_vars,
            output_dir=arg_value(args, "output-dir", DEFAULT_GRIS_MONTHLY_PDD_OUTPUT_DIR),
            netcdf_path=arg_value(args, "netcdf-path", ""),
            write_outputs=!has_flag(args, "no-output"),
        ),
        options=SimulationOptions(
            years=parse(Int, arg_value(args, "years", "10")),
            backend=arg_value(args, "backend", "threads"),
        ),
    )
    result = run!(sim; io=stdout)
    println("Status: $(result.status)")
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
