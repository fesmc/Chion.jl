#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion
using NCDatasets
include(joinpath(@__DIR__, "..", "shared", "script_helpers.jl"))
using .ChionExampleScriptHelpers

const DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_forcing_file_simulation")

function default_gris_forcing_file_path()
    for candidate in (
        get(ENV, "FORCING_PATH", ""),
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

function print_gris_api_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl [options]")
    println()
    println("Core options:")
    println("  --forcing-file=PATH          Prepared NetCDF forcing file")
    println("  --mask-threshold=VALUE       Apply MAR-style MSK threshold before running (default: 50)")
    println("  --output-dir=PATH            Output directory (default: examples/plots/gris_forcing_file_simulation)")
    println("  --netcdf-path=PATH           Output NetCDF path")
    println("  --no-output                  Skip summary/CSV file output")
    println("  --no-nc                      Skip NetCDF output")
    println("  --netcdf-vars=SPEC           NetCDF variables to write")
    println("  --model=NAME                 bessi|pdd (default: bessi)")
    println("  --ntot=N                     Maximum active layers (default: 20)")
    println("  --years=N                    Number of forcing years (default: 10)")
    println("  --backend=threads|cpu|gpu    Execution backend (default: threads)")
    println()
    println("Physics options:")
    println("  --albedo=NAME                constant|dynamic (default: dynamic)")
    println("  --densification=NAME         bessi|htessel (default: bessi)")
    println("  --fresh-snow-density=NAME    constant|parameterized (default: constant)")
    println("  --pdd-ddf-snow=VALUE         PDD snow degree-day factor in mmWE d-1 C-1 (default: 3)")
    println("  --pdd-ddf-ice=VALUE          PDD ice degree-day factor in mmWE d-1 C-1 (default: 8)")
    println("  --pdd-refreezing-fraction=X  PDD refreezing fraction (default: 0.6)")
    println()
    println("NetCDF variable groups:")
    println("  final, layers, history, monthly, step")
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_gris_api_help()
        return
    end

    default_forcing_file = default_gris_forcing_file_path()
    forcing_file = arg_value(args, "forcing-file", default_forcing_file)
    isempty(forcing_file) && error("Pass --forcing-file=PATH or set FORCING_PATH to the prepared forcing file.")

    loaded = load_forcing_file(forcing_file)
    loaded = apply_gris_mask(
        loaded,
        forcing_file;
        threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
    )
    model_name = lowercase(arg_value(args, "model", env_value("MODEL", "bessi")))
    model = if model_name == "bessi"
        build_model(
            :bessi,
            loaded.grid;
            Ntot=parse(Int, arg_value(args, "ntot", "20")),
            albedo=_albedo_scheme(lowercase(arg_value(args, "albedo", "dynamic"))),
            densification=_densification_scheme(lowercase(arg_value(args, "densification", "bessi"))),
            fresh_snow_density=_fresh_snow_scheme(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
        )
    elseif model_name == "pdd"
        build_model(
            :pdd,
            loaded.grid;
            ddf_snow=parse(Float64, arg_value(args, "pdd-ddf-snow", env_value("PDD_DDF_SNOW", "3.0"))),
            ddf_ice=parse(Float64, arg_value(args, "pdd-ddf-ice", env_value("PDD_DDF_ICE", "8.0"))),
            refreezing_fraction=parse(Float64, arg_value(args, "pdd-refreezing-fraction", env_value("PDD_REFREEZING_FRACTION", "0.6"))),
        )
    else
        error("Unsupported --model=$(model_name). Use bessi or pdd.")
    end
    println("Model: $(model_name)")
    default_netcdf_vars = model_name == "pdd" ? "final,history" : "all"
    netcdf_vars = has_flag(args, "no-nc") ? Symbol[] : arg_value(args, "netcdf-vars", default_netcdf_vars)
    println("NetCDF variables: ", isempty(netcdf_vars) ? "(none)" : netcdf_vars)
    sim = Simulation(
        model;
        forcing=loaded.forcing,
        output=OutputOptions(
            save=netcdf_vars,
            output_dir=arg_value(args, "output-dir", DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR),
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
