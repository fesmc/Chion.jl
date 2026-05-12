#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using NCDatasets
import Chion

include(joinpath(@__DIR__, "..", "shared", "script_helpers.jl"))
using .ChionExampleScriptHelpers

const SM = Chion
const DEFAULT_NC_PATH = get(
    ENV,
    "NC_PATH",
    get(ENV, "FORCING_PATH", "/p/projects/ou/labs/ai/Nils/MARv3.14.3-10km-daily-ERA5-2025.nc"),
)
const DEFAULT_OUT_DIR_EQUIL = joinpath(@__DIR__, "..", "plots", "gris_equilibrium")
const DEFAULT_NTOT = 20

function print_spinup_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_equilibrium.jl [options]")
    println()
    println("Options:")
    println("  --nc=PATH                    Prepared NetCDF forcing file")
    println("  --out-dir=PATH               Output directory (default: examples/plots/gris_equilibrium)")
    println("  --out-nc=PATH                Output NetCDF path (default: OUT_DIR/gris_equilibrium_final_state.nc)")
    println("  --no-nc                      Skip NetCDF output")
    println("  --no-output                  Skip summary/CSV file output and NetCDF")
    println("  --mask-threshold=VALUE       Apply MAR-style MSK threshold before running (default: 50)")
    println("  --ntot=N                     Chion maximum active layers (default: $(DEFAULT_NTOT))")
    println("  --max-years=N                Number of repeated forcing years (default: 10)")
    println("  --backend=threads|gpu        Execution backend (default: threads)")
    println("  --help                       Show this message")
end

function parse_spinup_config(args::Vector{String})
    nc_path = arg_value(args, "nc", DEFAULT_NC_PATH)
    isempty(nc_path) && error("Pass --nc=PATH or set NC_PATH/FORCING_PATH.")
    write_outputs = !has_flag(args, "no-output")
    return (
        nc_path=nc_path,
        out_dir=arg_value(args, "out-dir", DEFAULT_OUT_DIR_EQUIL),
        out_nc=arg_value(args, "out-nc", ""),
        write_outputs=write_outputs,
        write_netcdf=write_outputs && !has_flag(args, "no-nc"),
        mask_threshold=parse(Float64, arg_value(args, "mask-threshold", "50.0")),
        ntot=parse(Int, arg_value(args, "ntot", string(DEFAULT_NTOT))),
        max_years=parse(Int, arg_value(args, "max-years", "10")),
        backend=arg_value(args, "backend", "threads"),
    )
end

function read_gris_mask(path::AbstractString, ny::Int, nx::Int)
    mask = NCDataset(path) do ds
        haskey(ds, "MSK") ? Float64.(ds["MSK"][:]) : ones(Float64, ny, nx)
    end
    if length(mask) == nx * ny && ndims(mask) == 1
        mask = reshape(mask, nx, ny)
    end
    if size(mask) == (nx, ny)
        mask = permutedims(mask, (2, 1))
    end
    size(mask) == (ny, nx) || error("Unexpected MSK shape $(size(mask)); expected ($ny, $nx).")
    return Matrix{Float64}(mask)
end

function compact_gris_case(loaded, forcing_file::AbstractString, mask_threshold::Float64)
    ny = length(loaded.grid.y)
    nx = length(loaded.grid.x)
    mask = read_gris_mask(forcing_file, ny, nx)
    forcing = loaded.forcing

    rows = Int[]
    js = Int[]
    is = Int[]
    @inbounds for j in 1:ny, i in 1:nx
        row = i + (j - 1) * nx
        valid =
            isfinite(mask[j, i]) &&
            mask[j, i] >= mask_threshold &&
            all(isfinite, @view forcing.air_temperature[row, :]) &&
            all(isfinite, @view forcing.snowfall_rate[row, :]) &&
            all(isfinite, @view forcing.rainfall_rate[row, :]) &&
            all(isfinite, @view forcing.shortwave_down[row, :])
        if valid
            push!(rows, row)
            push!(js, j)
            push!(is, i)
        end
    end
    isempty(rows) && error("No valid GrIS columns found after applying MSK threshold $(mask_threshold).")

    compact_forcing = SM.SnowpackForcing(
        time_values=forcing.time_values,
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature[rows, :],
        snowfall_rate=forcing.snowfall_rate[rows, :],
        rainfall_rate=forcing.rainfall_rate[rows, :],
        shortwave_down=forcing.shortwave_down[rows, :],
        wind_speed=forcing.wind_speed[rows, :],
        q_lw_down=forcing.q_lw_down[rows, :],
        has_q_lw_down=forcing.has_q_lw_down[rows, :],
        q_sh=forcing.q_sh[rows, :],
        has_q_sh=forcing.has_q_sh[rows, :],
        q_lh=forcing.q_lh[rows, :],
        has_q_lh=forcing.has_q_lh[rows, :],
    )
    compact_grid = SM.SnowpackGrid(length(rows);
        x=loaded.grid.x,
        y=loaded.grid.y,
        js=js,
        is=is,
        mask=mask,
    )
    println("Applied GrIS mask: kept $(length(rows)) / $(nx * ny) columns")
    return (grid=compact_grid, forcing=compact_forcing)
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_spinup_help()
        return
    end

    config = parse_spinup_config(args)
    loaded = compact_gris_case(
        SM.load_forcing_file(config.nc_path),
        config.nc_path,
        config.mask_threshold,
    )
    model = SM.BESSIModel(loaded.grid; Ntot=config.ntot)
    nc_path = isempty(config.out_nc) ? joinpath(config.out_dir, "gris_equilibrium_final_state.nc") : config.out_nc
    sim = SM.Simulation(
        model;
        forcing=loaded.forcing,
        options=SM.SimulationOptions(
            name="gris_equilibrium",
            input_label=abspath(config.nc_path),
            years=config.max_years,
            backend=config.backend,
        ),
        output=SM.OutputOptions(
            save=config.write_netcdf ? "all" : "none",
            output_dir=config.out_dir,
            netcdf_path=nc_path,
            write_outputs=config.write_outputs,
        ),
    )
    result = SM.run!(sim; io=stdout)
    println("GrIS equilibrium public API run status: $(result.status)")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
