#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

include("gris_mar_case_backend.jl")

function print_gris_api_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_mar_case.jl [options]")
    println()
    println("Core options:")
    println("  --nc=PATH                    MAR NetCDF/HDF5 file")
    println("  --output-dir=PATH            Output directory (default: examples/plots/gris_mar_case)")
    println("  --netcdf-path=PATH           Output NetCDF path (default: OUTPUT_DIR/gris_mar_case_final_state.nc)")
    println("  --no-output                  Skip summary/CSV file output")
    println("  --no-nc                      Skip NetCDF output")
    println("  --netcdf-vars=SPEC           NetCDF variables to write: all, none, group names, or comma-separated variables")
    println("  --mask-threshold=VALUE       Minimum MSK value for GrIS cells (default: 50)")
    println("  --ntot=N                     Chion maximum active layers (default: 20)")
    println("  --cycles=N                   Number of forcing cycles to run (default: 10)")
    println("  --backend=threads|cpu|gpu    Execution backend (default: threads)")
    println("  --flip-turbulent-fluxes      Multiply SHF and LHF by -1 before forcing Chion")
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
    println("  julia --project=. examples/scripts/run_gris_mar_case.jl --nc=$(DEFAULT_GRIS_MAR_NC_PATH) --backend=threads --cycles=2 --no-output --no-nc")
    println("  julia --project=. examples/scripts/run_gris_mar_case.jl --backend=gpu --no-output --netcdf-vars=final,history")
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_gris_api_help()
        return
    end

    nc_path = arg_value(args, "nc", DEFAULT_GRIS_MAR_NC_PATH)
    isempty(nc_path) && error("Pass --nc=PATH or place the MAR file at $(DEFAULT_GRIS_MAR_NC_PATH).")

    physics = Chion.physics(
        albedo=Symbol(lowercase(arg_value(args, "albedo", "dynamic"))),
        densification=Symbol(lowercase(arg_value(args, "densification", "bessi"))),
        fresh_snow_density=Symbol(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
    )

    run = Chion.RunConfig(
        name="GrIS MAR case",
        input_label=abspath(nc_path),
        output_dir=arg_value(args, "output-dir", DEFAULT_GRIS_MAR_OUTPUT_DIR),
        netcdf_path=arg_value(args, "netcdf-path", ""),
        write_outputs=!has_flag(args, "no-output"),
        write_netcdf=!has_flag(args, "no-nc"),
        netcdf_variables=arg_value(args, "netcdf-vars", "all"),
        cycles=parse(Int, arg_value(args, "cycles", "10")),
        backend=arg_value(args, "backend", "threads"),
    )

    run_gris_mar_case(
        nc_path;
        io=stdout,
        mask_threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
        turbulent_flux_sign=has_flag(args, "flip-turbulent-fluxes") ? -1.0 : 1.0,
        ntot=parse(Int, arg_value(args, "ntot", "20")),
        run=run,
        physics=physics,
    )
    return
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
