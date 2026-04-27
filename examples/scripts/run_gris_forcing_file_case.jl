#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

const DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR = joinpath(@__DIR__, "..", "plots", "gris_forcing_file_case")

function default_gris_forcing_file_path()
    for candidate in (
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2025.nc",
        "/Users/niboch001/Downloads/MARv3.14.3-10km-daily-ERA5-2026.nc",
    )
        isfile(candidate) && return candidate
    end
    return ""
end

function arg_value(args::Vector{String}, name::String, default::String)
    prefix = "--" * name * "="
    for arg in args
        startswith(arg, prefix) && return arg[length(prefix)+1:end]
    end
    return default
end

function has_flag(args::Vector{String}, name::String)
    needle = "--" * name
    return any(arg -> arg == needle, args)
end

function print_gris_api_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl [options]")
    println()
    println("Core options:")
    println("  --forcing-file=PATH          Prepared HDF5/NetCDF forcing file")
    println("  --mask-threshold=VALUE       Apply MAR-style masking in this script before loading (default: 50)")
    println("  --output-dir=PATH            Output directory (default: examples/plots/gris_forcing_file_case)")
    println("  --netcdf-path=PATH           Output NetCDF path (default: OUTPUT_DIR/gris_forcing_file_case_final_state.nc)")
    println("  --no-output                  Skip summary/CSV file output")
    println("  --no-nc                      Skip NetCDF output")
    println("  --netcdf-vars=SPEC           NetCDF variables to write: all, none, group names, or comma-separated variables")
    println("  --ntot=N                     Chion maximum active layers (default: 20)")
    println("  --cycles=N                   Number of forcing cycles to run (default: 10)")
    println("  --backend=threads|cpu|gpu    Execution backend (default: threads)")
    println()
    println("Physics options:")
    println("  --albedo=NAME                constant|dynamic|legacy|bessi (default: dynamic)")
    println("  --densification=NAME         bessi|htessel (default: bessi)")
    println("  --fresh-snow-density=NAME    constant|parameterized|bessi|htessel (default: constant)")
    println()
    println("NetCDF variable groups:")
    println("  final, layers, history, monthly, step")
    println()
    println("Examples:")
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl --forcing-file=$(default_gris_forcing_file_path()) --backend=threads --cycles=2 --no-output --no-nc")
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl --backend=gpu --no-output --netcdf-vars=final,history")
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_gris_api_help()
        return
    end

    default_forcing_file = default_gris_forcing_file_path()
    forcing_file = arg_value(args, "forcing-file", default_forcing_file)
    isempty(forcing_file) && error("Pass --forcing-file=PATH or place the prepared forcing file at $(default_forcing_file).")

    physics = Chion.physics(
        albedo=Symbol(lowercase(arg_value(args, "albedo", "dynamic"))),
        densification=Symbol(lowercase(arg_value(args, "densification", "bessi"))),
        fresh_snow_density=Symbol(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
    )

    problem = Chion.load_gris_forcing_file_problem(
        forcing_file;
        mask_threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
        ntot=parse(Int, arg_value(args, "ntot", "20")),
        physics=physics,
    )
    for note in problem.notes
        println(note)
    end
    fill!(problem.domain.N, 0)
    fill!(problem.domain.mass, 0.0)
    fill!(problem.domain.mass_w, 0.0)
    fill!(problem.domain.density, 0.0)
    fill!(problem.domain.temperature, problem.domain.c.T0)
    fill!(problem.domain.mass_base, 0.0)
    fill!(problem.domain.smb_ice, 0.0)
    fill!(problem.domain.runoff, 0.0)
    fill!(problem.domain.snow_cover, 0.0)
    fill!(problem.domain.albedo_dynamic, problem.domain.c.alpha_ice)
    fill!(problem.domain.Tsrf, problem.domain.c.T0)
    Chion.run!(
        problem;
        save=has_flag(args, "no-nc") ? Symbol[] : arg_value(args, "netcdf-vars", "all"),
        output_dir=arg_value(args, "output-dir", DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR),
        netcdf_path=arg_value(args, "netcdf-path", ""),
        write_outputs=!has_flag(args, "no-output"),
        cycles=parse(Int, arg_value(args, "cycles", "10")),
        backend=arg_value(args, "backend", "threads"),
        history_stride=1,
        io=stdout,
    )
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
