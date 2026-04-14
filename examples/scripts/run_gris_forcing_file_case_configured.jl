#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Chion

include("gris_forcing_file_case_backend.jl")

const SM = Chion.SnowpackModel

const DEFAULT_CONFIG = (
    name="GrIS forcing file case",
    forcing_path=DEFAULT_GRIS_FORCING_FILE_PATH,
    output_dir=joinpath(@__DIR__, "..", "plots", "gris_forcing_file_case_configured"),
    netcdf_path="",
    write_outputs=false,
    write_netcdf=false,
    netcdf_variables="all",
    mask_threshold=50.0,
    ntot=20,
    cycles=1000,
    history_stride=1,
    backend=:gpu,
    albedo=:dynamic,
    densification=:bessi,
    fresh_snow_density=:constant,
)

function env_string(name::AbstractString)
    value = get(ENV, String(name), "")
    return isempty(value) ? nothing : value
end

function env_int(name::AbstractString)
    value = env_string(name)
    return isnothing(value) ? nothing : parse(Int, value)
end

function env_float(name::AbstractString)
    value = env_string(name)
    return isnothing(value) ? nothing : parse(Float64, value)
end

function env_bool(name::AbstractString)
    value = env_string(name)
    if isnothing(value)
        return nothing
    end
    lower = lowercase(value)
    lower in ("1", "true", "yes", "on") && return true
    lower in ("0", "false", "no", "off") && return false
    error("Invalid boolean value for $(name): $(value)")
end

function env_symbol(name::AbstractString)
    value = env_string(name)
    return isnothing(value) ? nothing : Symbol(value)
end

function config_from_env(base=DEFAULT_CONFIG)
    overrides = Pair{Symbol, Any}[]
    for item in (
        :forcing_path => env_string("FORCING_PATH"),
        :output_dir => env_string("OUTPUT_DIR"),
        :netcdf_path => env_string("NETCDF_PATH"),
        :backend => env_symbol("BACKEND"),
        :cycles => env_int("CYCLES"),
        :history_stride => env_int("HISTORY_STRIDE"),
        :write_outputs => env_bool("WRITE_OUTPUTS"),
        :write_netcdf => env_bool("WRITE_NETCDF"),
        :netcdf_variables => env_string("NETCDF_VARIABLES"),
        :mask_threshold => env_float("MASK_THRESHOLD"),
        :ntot => env_int("NTOT"),
        :albedo => env_symbol("ALBEDO"),
        :densification => env_symbol("DENSIFICATION"),
        :fresh_snow_density => env_symbol("FRESH_SNOW_DENSITY"),
    )
        key, value = item
        isnothing(value) || push!(overrides, key => value)
    end
    return isempty(overrides) ? base : merge(base, NamedTuple(overrides))
end

function main(config=DEFAULT_CONFIG)
    physics = Chion.physics(
        albedo=config.albedo,
        densification=config.densification,
        fresh_snow_density=config.fresh_snow_density,
    )

    run = Chion.RunConfig(
        name=config.name,
        output_dir=config.output_dir,
        netcdf_path=config.netcdf_path,
        write_outputs=config.write_outputs,
        write_netcdf=config.write_netcdf,
        netcdf_variables=config.netcdf_variables,
        cycles=config.cycles,
        history_stride=config.history_stride,
        backend=config.backend,
    )

    return run_gris_forcing_file_case(
        config.forcing_path;
        io=stdout,
        mask_threshold=config.mask_threshold,
        ntot=config.ntot,
        run=run,
        physics=physics,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(config_from_env())
end
