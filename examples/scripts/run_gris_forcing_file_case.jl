#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include("gris_forcing_file_case_backend.jl")

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
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl --forcing-file=$(DEFAULT_GRIS_FORCING_FILE_PATH) --backend=threads --cycles=2 --no-output --no-nc")
    println("  julia --project=. examples/scripts/run_gris_forcing_file_case.jl --backend=gpu --no-output --netcdf-vars=final,history")
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_gris_api_help()
        return
    end

    forcing_file = arg_value(args, "forcing-file", DEFAULT_GRIS_FORCING_FILE_PATH)
    isempty(forcing_file) && error("Pass --forcing-file=PATH or place the prepared forcing file at $(DEFAULT_GRIS_FORCING_FILE_PATH).")

    physics = Chion.physics(
        albedo=Symbol(lowercase(arg_value(args, "albedo", "dynamic"))),
        densification=Symbol(lowercase(arg_value(args, "densification", "bessi"))),
        fresh_snow_density=Symbol(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
    )

    run = Chion.RunConfig(
        name="GrIS forcing file case",
        output_dir=arg_value(args, "output-dir", DEFAULT_GRIS_FORCING_FILE_OUTPUT_DIR),
        netcdf_path=arg_value(args, "netcdf-path", ""),
        write_outputs=!has_flag(args, "no-output"),
        write_netcdf=!has_flag(args, "no-nc"),
        netcdf_variables=arg_value(args, "netcdf-vars", "all"),
        cycles=parse(Int, arg_value(args, "cycles", "10")),
        backend=arg_value(args, "backend", "threads"),
    )

    run_gris_forcing_file_case(
        forcing_file;
        io=stdout,
        mask_threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
        ntot=parse(Int, arg_value(args, "ntot", "20")),
        run=run,
        physics=physics,
    )
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
