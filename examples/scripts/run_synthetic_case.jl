#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

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

function print_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_synthetic_case.jl [options]")
    println()
    println("Options:")
    println("  --backend=threads|cpu|gpu    Execution backend (default: threads)")
    println("  --cycles=N                   Number of forcing cycles to run (default: 3)")
    println("  --nx=N                       Synthetic multi-column grid size in x (default: 2)")
    println("  --ny=N                       Synthetic multi-column grid size in y (default: 2)")
    println("  --output-dir=PATH            Output directory")
    println("  --netcdf-path=PATH           Output NetCDF path")
    println("  --no-output                  Skip summary/CSV output")
    println("  --no-nc                      Skip NetCDF output")
    println("  --netcdf-vars=SPEC           NetCDF variables to write")
    println("  --albedo=NAME                constant|dynamic|legacy|bessi")
    println("  --densification=NAME         bessi|htessel")
    println("  --fresh-snow-density=NAME    constant|parameterized|bessi|htessel")
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_help()
        return
    end

    case = Chion.synthetic_case(
        variant=:multi_column,
        physics=Chion.physics(
            albedo=Symbol(lowercase(arg_value(args, "albedo", "dynamic"))),
            densification=Symbol(lowercase(arg_value(args, "densification", "bessi"))),
            fresh_snow_density=Symbol(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
        ),
        ntot=5,
        ntime=12,
        nx=parse(Int, arg_value(args, "nx", "2")),
        ny=parse(Int, arg_value(args, "ny", "2")),
        run=Chion.RunConfig(
            name="synthetic_case",
            output_dir=arg_value(args, "output-dir", joinpath(@__DIR__, "..", "plots", "synthetic_case")),
            netcdf_path=arg_value(args, "netcdf-path", ""),
            write_outputs=!has_flag(args, "no-output"),
            write_netcdf=!has_flag(args, "no-nc"),
            netcdf_variables=arg_value(args, "netcdf-vars", "final,history"),
            cycles=parse(Int, arg_value(args, "cycles", "3")),
            backend=arg_value(args, "backend", "threads"),
        ),
    )
    Chion.run_case(case; io=stdout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
