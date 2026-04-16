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

function seed_surface_domain!(
    domain::Chion.SnowpackDomain;
    surface_mass,
    density::Real=320.0,
    temperature_c::Real=-12.0,
)
    fill!(domain.N, 1)
    fill!(domain.mass, 0.0)
    fill!(domain.mass_w, 0.0)
    fill!(domain.density, 0.0)
    fill!(domain.temperature, domain.c.T0)
    fill!(domain.mass_base, 0.0)
    fill!(domain.smb_ice, 0.0)
    fill!(domain.runoff, 0.0)
    fill!(domain.snow_cover, 0.0)
    fill!(domain.albedo_dynamic, domain.c.alpha_dry)

    surface_values = surface_mass isa AbstractVector ? surface_mass : fill(Float64(surface_mass), domain.ncol)
    length(surface_values) == domain.ncol || error("`surface_mass` must match the domain column count.")
    @inbounds for idx in eachindex(surface_values)
        if surface_values[idx] > 0
            domain.N[idx] = 1
            domain.mass[1, idx] = Float64(surface_values[idx])
            domain.density[1, idx] = Float64(density)
            domain.temperature[1, idx] = domain.c.T0 + Float64(temperature_c)
            domain.Tsrf[idx] = domain.temperature[1, idx]
        else
            domain.Tsrf[idx] = domain.c.T0
        end
    end

    Chion.compute_auxiliary!(domain)
    return domain
end

function regular_layout(nx::Int, ny::Int)
    x = collect(1:nx)
    y = collect(1:ny)
    js = [j for j in 1:ny for _ in 1:nx]
    is = [i for _ in 1:ny for i in 1:nx]
    mask = ones(Float64, ny, nx)
    return Chion.GridLayout(Float64.(x), Float64.(y), js, is, mask)
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_help()
        return
    end

    nx = parse(Int, arg_value(args, "nx", "2"))
    ny = parse(Int, arg_value(args, "ny", "2"))
    ncol = nx * ny
    physics = Chion.physics(
        albedo=Symbol(lowercase(arg_value(args, "albedo", "dynamic"))),
        densification=Symbol(lowercase(arg_value(args, "densification", "bessi"))),
        fresh_snow_density=Symbol(lowercase(arg_value(args, "fresh-snow-density", "constant"))),
    )

    domain = Chion.SnowpackDomain(c=physics, Ntot=5, ncol=ncol)
    surface_mass = [
        220.0 + 35.0 * ((i - 1) / max(nx - 1, 1)) + 20.0 * ((j - 1) / max(ny - 1, 1))
        for j in 1:ny for i in 1:nx
    ]
    seed_surface_domain!(domain; surface_mass=surface_mass, density=320.0, temperature_c=-12.0)

    ntime = 12
    step = collect(1:ntime)
    forcing = Chion.ForcingData(
        dt_days=fill(1.0, ntime),
        air_temperature_c=(-18.0 .+ 0.02 .* step),
        snowfall_mm_day=(fill(0.15, ntime) .+ 0.1 .* sin.(0.3 .* step)),
        rainfall_mm_day=zeros(Float64, ntime),
        shortwave_down=fill(160.0, ntime),
        wind_speed=4.5 .+ 0.4 .* sin.(0.05 .* step),
        ncol=ncol,
    )

    layout = regular_layout(nx, ny)
    save = has_flag(args, "no-nc") ? Symbol[] : arg_value(args, "netcdf-vars", "final,history")
    result = Chion.run!(
        domain,
        forcing;
        layout=layout,
        save=save,
        output_dir=arg_value(args, "output-dir", joinpath(@__DIR__, "..", "plots", "synthetic_case")),
        netcdf_path=arg_value(args, "netcdf-path", ""),
        write_outputs=!has_flag(args, "no-output"),
        cycles=parse(Int, arg_value(args, "cycles", "3")),
        backend=arg_value(args, "backend", "threads"),
        history_stride=1,
        io=stdout,
    )
    println("Status: $(result.status)")
    println("Cycles: $(length(result.history))")
    println("Summary: $(result.summary_path)")
    println("History CSV: $(result.history_csv_path)")
    println("NetCDF: $(result.netcdf_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
