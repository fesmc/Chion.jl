#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using Printf
using Statistics
using Base.Threads
using CUDA
import Libdl
using NCDatasets
using Chion

const SM = Chion
const DEFAULT_OUT_DIR_EQUIL = joinpath(@__DIR__, "..", "plots", "gris_equilibrium")
const NC_NOERR = 0
const NC_CLOBBER = 0x0000
const NC_NETCDF4 = 0x1000
const NC_GLOBAL = -1
const NC_FLOAT = 5
const NC_DOUBLE = 6
const NC_INT = 4
const NC_UNLIMITED = 0

@inline valid_or(default::Float64, x::Float64) = isfinite(x) ? x : default
@inline mmwe_day_to_kgm2s(x::Float64) = isfinite(x) ? max(x, 0.0) / 86_400.0 : 0.0

function resolve_libnetcdf()
    haskey(ENV, "NETCDF_LIB") && return ENV["NETCDF_LIB"]
    lib = Libdl.find_library(["netcdf", "libnetcdf"])
    isempty(lib) || return lib
    isdefined(NCDatasets, :libnetcdf) && return String(getproperty(NCDatasets, :libnetcdf))
    error("Could not locate libnetcdf. Set NETCDF_LIB to the shared library path.")
end

const LIBNETCDF = resolve_libnetcdf()

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
            snow_mass=0.0,
            liquid_mass=0.0,
            wet_mass=0.0,
            thickness=0.0,
            bulk_density=NaN,
            base_mass=mass_base,
            smb_ice=smb_ice,
            runoff=runoff,
            snow_cover=0.0,
            surface_temperature_c=NaN,
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
        snow_mass=snow_mass,
        liquid_mass=liquid_mass,
        wet_mass=wet_mass,
        thickness=thickness,
        bulk_density=bulk_density,
        base_mass=mass_base,
        smb_ice=smb_ice,
        runoff=runoff,
        snow_cover=snow_cover,
        surface_temperature_c=Float64(temperature[1]) - T0,
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

mutable struct TimingStats
    totals::Dict{Symbol, Float64}
    counts::Dict{Symbol, Int}
end

TimingStats() = TimingStats(Dict{Symbol, Float64}(), Dict{Symbol, Int}())

function add_timing!(stats::TimingStats, key::Symbol, dt_sec::Float64, count::Int=1)
    stats.totals[key] = get(stats.totals, key, 0.0) + dt_sec
    stats.counts[key] = get(stats.counts, key, 0) + count
    return dt_sec
end

function time_block!(stats::TimingStats, key::Symbol, f::F) where {F<:Function}
    t0 = time_ns()
    value = f()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

time_block!(f::F, stats::TimingStats, key::Symbol) where {F<:Function} = time_block!(stats, key, f)

function timing_rows(stats::TimingStats; total_wall_sec::Union{Nothing, Float64}=nothing)
    rows = NamedTuple[]
    total = sum(values(stats.totals))
    share_total = isnothing(total_wall_sec) ? total : total_wall_sec
    for key in keys(stats.totals)
        dt = stats.totals[key]
        count = stats.counts[key]
        push!(rows, (
            key = key,
            total_sec = dt,
            count = count,
            mean_sec = count > 0 ? dt / count : NaN,
            share_pct = share_total > 0.0 ? 100.0 * dt / share_total : 0.0,
        ))
    end
    sort!(rows; by=row -> row.total_sec, rev=true)
    return rows, total
end

function print_timing_summary(io::IO, stats::TimingStats; total_wall_sec::Union{Nothing, Float64}=nothing)
    rows, total = timing_rows(stats; total_wall_sec=total_wall_sec)
    println(io, "Timing summary")
    println(io, @sprintf("  %-24s %12s %9s %12s %10s", "stage", "total [s]", "share", "mean [ms]", "count"))
    for row in rows
        println(
            io,
            @sprintf(
                "  %-24s %12.3f %8.1f%% %12.3f %10d",
                String(row.key),
                row.total_sec,
                row.share_pct,
                row.mean_sec * 1.0e3,
                row.count,
            ),
        )
    end
    if !isnothing(total_wall_sec)
        unaccounted = max(total_wall_sec - total, 0.0)
        println(
            io,
            @sprintf(
                "  %-24s %12.3f %8.1f%% %12s %10s",
                "unaccounted",
                unaccounted,
                total_wall_sec > 0.0 ? 100.0 * unaccounted / total_wall_sec : 0.0,
                "",
                "",
            ),
        )
        println(io, @sprintf("  %-24s %12.3f", "total_accounted", total))
        println(io, @sprintf("  %-24s %12.3f", "run_wall_total", total_wall_sec))
    else
        println(io, @sprintf("  %-24s %12.3f", "total_accounted", total))
    end
    return
end

function print_spinup_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_equilibrium.jl [options]")
    println()
    println("Options:")
    println("  --nc=PATH                    Prepared HDF5/NetCDF forcing file")
    println("  --out-dir=PATH               Output directory (default: examples/plots/gris_equilibrium)")
    println("  --out-nc=PATH                Output NetCDF path (default: OUT_DIR/gris_equilibrium_final_state.nc)")
    println("  --no-nc                      Skip NetCDF output and NetCDF-only step diagnostics for cleaner timing")
    println("  --no-output                  Skip all file output (implies --no-nc) for clean timing runs")
    println("  --ntot=N                     Chion maximum active layers (default: 80)")
    println("  --max-cycles=N               Maximum forcing-cycle repeats (default: 10)")
    println("  --backend=threads|gpu        Execution backend (default: threads)")
    println("  --help                       Show this message")
end

function normalize_backend(name::AbstractString)
    backend = Symbol(lowercase(strip(name)))
    backend in (:threads, :gpu) || error("Unsupported backend '$name'. Use `threads` or `gpu`.")
    return backend
end

@inline summary_backend(backend::Symbol) = backend == :gpu ? :kernelabstractions : :threads

function parse_spinup_config(args::Vector{String})
    nc_path = arg_value(args, "nc", DEFAULT_NC_PATH)
    isempty(nc_path) && error("Pass --nc=PATH or place the prepared forcing file at $(DEFAULT_NC_PATH).")
    write_outputs = !has_flag(args, "no-output")
    return (
        nc_path = nc_path,
        out_dir = arg_value(args, "out-dir", DEFAULT_OUT_DIR_EQUIL),
        out_nc = arg_value(args, "out-nc", ""),
        write_outputs = write_outputs,
        write_netcdf = write_outputs && !has_flag(args, "no-nc"),
        ntot = parse(Int, arg_value(args, "ntot", "20")),
        max_cycles = parse(Int, arg_value(args, "max-cycles", "10")),
        backend = normalize_backend(arg_value(args, "backend", "threads")),
    )
end

function build_annual_output_schedule(time_values::Vector{DateTime})
    ntime = length(time_values)
    write_output = falses(ntime)
    output_slot = zeros(Int, ntime)
    source_indices = Int32[]
    source_codes = Int32[]

    years = unique(year.(time_values))
    for (slot, yr) in enumerate(years)
        last_t = findlast(t -> year(time_values[t]) == yr, eachindex(time_values))
        isnothing(last_t) && error("Could not determine the last timestep for source year $yr.")
        write_output[last_t] = true
        output_slot[last_t] = slot
        push!(source_indices, Int32(last_t))
        ts = time_values[last_t]
        push!(source_codes, Int32(year(ts) * 1000000 + month(ts) * 10000 + day(ts) * 100 + hour(ts)))
    end

    return (
        write_output=write_output,
        output_slot=output_slot,
        source_indices=source_indices,
        source_codes=source_codes,
        years=years,
    )
end

function read_full_timeseries_3d(nc_path::AbstractString, varname::AbstractString, shapes::Dict{String, Vector{Int}})
    data = read_hdf5_full(nc_path, varname, shapes)
    if ndims(data) == 4
        size(data, 2) == 1 || error("Variable '$varname' has an unexpected non-singleton vertical dimension.")
        return dropdims(data; dims=(2,))
    elseif ndims(data) == 3
        return data
    else
        error("Variable '$varname' does not have a supported timeseries layout.")
    end
end

function read_first_available_timeseries_3d(
    nc_path::AbstractString,
    candidate_names::Vector{String},
    shapes::Dict{String, Vector{Int}},
)
    for name in candidate_names
        if haskey(shapes, name)
            return (name=name, data=read_full_timeseries_3d(nc_path, name, shapes))
        end
    end
    return nothing
end

@inline function nc_check(code::Integer)
    if code != NC_NOERR
        msg = unsafe_string(ccall((:nc_strerror, LIBNETCDF), Cstring, (Cint,), code))
        error("NetCDF error: $msg")
    end
    return
end

function nc_create(path::AbstractString)
    ncid = Ref{Cint}()
    nc_check(ccall((:nc_create, LIBNETCDF), Cint, (Cstring, Cint, Ref{Cint}), path, NC_CLOBBER | NC_NETCDF4, ncid))
    return ncid[]
end

function nc_close(ncid::Cint)
    nc_check(ccall((:nc_close, LIBNETCDF), Cint, (Cint,), ncid))
    return
end

function nc_enddef(ncid::Cint)
    nc_check(ccall((:nc_enddef, LIBNETCDF), Cint, (Cint,), ncid))
    return
end

function nc_redef(ncid::Cint)
    nc_check(ccall((:nc_redef, LIBNETCDF), Cint, (Cint,), ncid))
    return
end

function nc_def_dim(ncid::Cint, name::AbstractString, len::Integer)
    dimid = Ref{Cint}()
    nc_check(ccall((:nc_def_dim, LIBNETCDF), Cint, (Cint, Cstring, Csize_t, Ref{Cint}), ncid, name, len, dimid))
    return dimid[]
end

function nc_def_var(ncid::Cint, name::AbstractString, xtype::Integer, dimids::Vector{Cint})
    varid = Ref{Cint}()
    nc_check(ccall((:nc_def_var, LIBNETCDF), Cint, (Cint, Cstring, Cint, Cint, Ptr{Cint}, Ref{Cint}), ncid, name, xtype, length(dimids), dimids, varid))
    return varid[]
end

function nc_put_att_text(ncid::Cint, varid::Integer, name::AbstractString, value::AbstractString)
    nc_check(ccall((:nc_put_att_text, LIBNETCDF), Cint, (Cint, Cint, Cstring, Csize_t, Cstring), ncid, Cint(varid), name, sizeof(value), value))
    return
end

function nc_put_att_float(ncid::Cint, varid::Integer, name::AbstractString, value::Float32)
    ref = Ref{Float32}(value)
    nc_check(ccall((:nc_put_att_float, LIBNETCDF), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ptr{Cfloat}), ncid, Cint(varid), name, NC_FLOAT, 1, ref))
    return
end

function nc_put_var_double(ncid::Cint, varid::Integer, data::Vector{Float64})
    nc_check(ccall((:nc_put_var_double, LIBNETCDF), Cint, (Cint, Cint, Ptr{Cdouble}), ncid, Cint(varid), data))
    return
end

function nc_put_var_int(ncid::Cint, varid::Integer, data::Vector{Int32})
    nc_check(ccall((:nc_put_var_int, LIBNETCDF), Cint, (Cint, Cint, Ptr{Cint}), ncid, Cint(varid), data))
    return
end

function nc_put_var_float_2d(ncid::Cint, varid::Integer, data::AbstractMatrix{<:Real})
    buf = permutedims(Float32.(data), (2, 1))
    nc_check(ccall((:nc_put_var_float, LIBNETCDF), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_int_2d(ncid::Cint, varid::Integer, data::AbstractMatrix{<:Integer})
    buf = permutedims(Int32.(data), (2, 1))
    nc_check(ccall((:nc_put_var_int, LIBNETCDF), Cint, (Cint, Cint, Ptr{Cint}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_3d(ncid::Cint, varid::Integer, data::Array{Float64, 3})
    buf = permutedims(Float32.(data), (3, 2, 1))
    nc_check(ccall((:nc_put_var_float, LIBNETCDF), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_1d(ncid::Cint, varid::Integer, data::Vector{Float64})
    buf = Float32.(data)
    nc_check(ccall((:nc_put_var_float, LIBNETCDF), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_vara_float_4d_step_layer_yx(
    ncid::Cint,
    varid::Integer,
    step_index::Integer,
    data::Array{Float64, 3},
)
    start = Csize_t[Csize_t(step_index - 1), Csize_t(0), Csize_t(0), Csize_t(0)]
    count = Csize_t[Csize_t(1), Csize_t(size(data, 1)), Csize_t(size(data, 2)), Csize_t(size(data, 3))]
    buf = permutedims(Float32.(data), (3, 2, 1))
    nc_check(
        ccall(
            (:nc_put_vara_float, LIBNETCDF),
            Cint,
            (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Cfloat}),
            ncid,
            Cint(varid),
            start,
            count,
            buf,
        ),
    )
    return
end

function nc_put_vara_float_3d_step_yx(
    ncid::Cint,
    varid::Integer,
    step_index::Integer,
    data::AbstractMatrix{<:Real},
)
    start = Csize_t[Csize_t(step_index - 1), Csize_t(0), Csize_t(0)]
    count = Csize_t[Csize_t(1), Csize_t(size(data, 1)), Csize_t(size(data, 2))]
    buf = permutedims(Float32.(data), (2, 1))
    nc_check(
        ccall(
            (:nc_put_vara_float, LIBNETCDF),
            Cint,
            (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Cfloat}),
            ncid,
            Cint(varid),
            start,
            count,
            buf,
        ),
    )
    return
end

function define_nc_output_variable(ncid::Cint, dimids::Vector{Cint}, name::AbstractString, long_name::AbstractString, units::AbstractString)
    varid = nc_def_var(ncid, name, NC_FLOAT, dimids)
    nc_put_att_text(ncid, varid, "long_name", long_name)
    nc_put_att_text(ncid, varid, "units", units)
    nc_put_att_float(ncid, varid, "_FillValue", Float32(NaN))
    return varid
end

function scatter_to_grid(values::Vector{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    out = fill(NaN, grid_shape)
    @inbounds for idx in eachindex(values)
        out[js[idx], is[idx]] = values[idx]
    end
    return out
end

function allocate_summary_buffers(n::Int)
    thickness = Vector{Float64}(undef, n)
    wet_mass = Vector{Float64}(undef, n)
    bulk_density = Vector{Float64}(undef, n)
    base_mass = Vector{Float64}(undef, n)
    smb_ice = Vector{Float64}(undef, n)
    liquid_water = Vector{Float64}(undef, n)
    runoff = Vector{Float64}(undef, n)
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
        smb_ice=smb_ice,
        liquid_water=liquid_water,
        runoff=runoff,
    )
end

function allocate_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int)
    thickness = similar(domain.mass, Float64, n)
    wet_mass = similar(domain.mass, Float64, n)
    bulk_density = similar(domain.mass, Float64, n)
    base_mass = similar(domain.mass, Float64, n)
    smb_ice = similar(domain.mass, Float64, n)
    liquid_water = similar(domain.mass, Float64, n)
    runoff = similar(domain.mass, Float64, n)
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
        smb_ice=smb_ice,
        liquid_water=liquid_water,
        runoff=runoff,
    )
end

function allocate_cycle_summary_buffers(n::Int)
    thickness = Vector{Float64}(undef, n)
    wet_mass = Vector{Float64}(undef, n)
    bulk_density = Vector{Float64}(undef, n)
    base_mass = Vector{Float64}(undef, n)
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
    )
end

function allocate_cycle_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int)
    thickness = similar(domain.mass, Float64, n)
    wet_mass = similar(domain.mass, Float64, n)
    bulk_density = similar(domain.mass, Float64, n)
    base_mass = similar(domain.mass, Float64, n)
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
    )
end

function summarize_columns!(summary, domain::SM.SnowpackDomain; backend::Symbol=:threads, device_summary=nothing)
    if backend == :kernelabstractions
        if isnothing(device_summary)
            device_summary = allocate_summary_buffers(domain, length(summary.thickness))
        end
        SM.summarize_domain_state!(
            device_summary.thickness,
            device_summary.wet_mass,
            device_summary.bulk_density,
            device_summary.base_mass,
            device_summary.smb_ice,
            device_summary.liquid_water,
            device_summary.runoff,
            domain
            )
        copyto!(summary.thickness, device_summary.thickness)
        copyto!(summary.wet_mass, device_summary.wet_mass)
        copyto!(summary.bulk_density, device_summary.bulk_density)
        copyto!(summary.base_mass, device_summary.base_mass)
        copyto!(summary.smb_ice, device_summary.smb_ice)
        copyto!(summary.liquid_water, device_summary.liquid_water)
        copyto!(summary.runoff, device_summary.runoff)
    else
        SM.summarize_domain_state!(
            summary.thickness,
            summary.wet_mass,
            summary.bulk_density,
            summary.base_mass,
            summary.smb_ice,
            summary.liquid_water,
            summary.runoff,
            domain
            )
    end
    return summary
end

function summarize_cycle_columns!(summary, domain::SM.SnowpackDomain; backend::Symbol=:threads, device_summary=nothing)
    if backend == :kernelabstractions
        if isnothing(device_summary)
            device_summary = allocate_cycle_summary_buffers(domain, length(summary.thickness))
        end
        SM.summarize_cycle_state!(
            device_summary.thickness,
            device_summary.wet_mass,
            device_summary.bulk_density,
            device_summary.base_mass,
            domain
            )
        copyto!(summary.thickness, device_summary.thickness)
        copyto!(summary.wet_mass, device_summary.wet_mass)
        copyto!(summary.bulk_density, device_summary.bulk_density)
        copyto!(summary.base_mass, device_summary.base_mass)
    else
        SM.summarize_cycle_state!(
            summary.thickness,
            summary.wet_mass,
            summary.bulk_density,
            summary.base_mass,
            domain
            )
    end
    return summary
end

function cpu_summary(summary)
    return (
        thickness=Array(summary.thickness),
        wet_mass=Array(summary.wet_mass),
        bulk_density=Array(summary.bulk_density),
        base_mass=Array(summary.base_mass),
        smb_ice=Array(summary.smb_ice),
        liquid_water=Array(summary.liquid_water),
        runoff=Array(summary.runoff),
    )
end

function cpu_cycle_summary(summary)
    return (
        thickness=Array(summary.thickness),
        wet_mass=Array(summary.wet_mass),
        bulk_density=Array(summary.bulk_density),
        base_mass=Array(summary.base_mass),
    )
end

function domain_mean_vector(v::Vector{Float64})
    total = 0.0
    n = 0
    for x in v
        if isfinite(x)
            total += x
            n += 1
        end
    end
    return n == 0 ? NaN : total / n
end

function domain_mean_abs_vector(v::Vector{Float64})
    total = 0.0
    n = 0
    for x in v
        if isfinite(x)
            total += abs(x)
            n += 1
        end
    end
    return n == 0 ? NaN : total / n
end

function domain_max_abs_vector(v::Vector{Float64})
    vmax = 0.0
    hasval = false
    for x in v
        if isfinite(x)
            vmax = max(vmax, abs(x))
            hasval = true
        end
    end
    return hasval ? vmax : NaN
end

function make_cycle_record(
    cycle::Int,
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    delta_thickness,
    delta_wet_mass,
    delta_base_mass,
)
    n = length(thickness)
    return (
        cycle=cycle,
        mean_thickness=sum(thickness) / n,
        mean_wet_mass=sum(wet_mass) / n,
        mean_bulk_density=sum(bulk_density) / n,
        mean_base_mass=sum(base_mass) / n,
        mean_signed_delta_thickness=sum(delta_thickness) / n,
        mean_abs_delta_thickness=sum(abs, delta_thickness) / n,
        max_abs_delta_thickness=maximum(abs, delta_thickness),
        mean_signed_delta_wet_mass=sum(delta_wet_mass) / n,
        mean_abs_delta_wet_mass=sum(abs, delta_wet_mass) / n,
        max_abs_delta_wet_mass=maximum(abs, delta_wet_mass),
        mean_signed_delta_base_mass=sum(delta_base_mass) / n,
        mean_abs_delta_base_mass=sum(abs, delta_base_mass) / n,
        max_abs_delta_base_mass=maximum(abs, delta_base_mass),
    )
end

function cycle_log_line(record)
    return @sprintf(
        "cycle=%d mean_th=%.5f m mean_swe=%.5f mmWE mean_base=%.5f mmWE mean_abs_dth=%.5f m mean_abs_dswe=%.5f mmWE mean_abs_dbase=%.5f mmWE",
        record.cycle,
        record.mean_thickness,
        record.mean_wet_mass,
        record.mean_base_mass,
        record.mean_abs_delta_thickness,
        record.mean_abs_delta_wet_mass,
        record.mean_abs_delta_base_mass,
    )
end

function run_spinup_cycles_no_netcdf!(
    timings::TimingStats,
    config,
    domain::SM.SnowpackDomain,
    workspace::SM.ColumnarStepWorkspace,
    step_fields::SM.SnowpackForcing,
    nvalid::Int,
)
    config.backend == :threads || error("CPU fast path requires `--backend=threads`.")
    !config.write_netcdf || error("CPU fast path only applies when NetCDF output is disabled.")

    ntime = size(step_fields.air_temperature, 2)
    ncol = SM.column_count(domain)
    prev = allocate_cycle_summary_buffers(nvalid)
    final = allocate_cycle_summary_buffers(nvalid)
    time_block!(timings, :summarize_columns_initial) do
        summarize_cycle_columns!(prev, domain; backend=:threads)
    end

    previous_cycle_smb_ice_vec = config.write_outputs ? copy(domain.smb_ice) : Float64[]
    history = NamedTuple[]
    status = :max_cycles
    last_delta_thickness_vec = fill(NaN, nvalid)
    last_delta_wet_mass_vec = fill(NaN, nvalid)
    last_delta_base_mass_vec = fill(NaN, nvalid)
    last_delta_ice_sheet_smb_vec = fill(NaN, nvalid)
    simulation_wall_t0 = time_ns()

    for cycle in 1:config.max_cycles
        time_counted_block!(timings, :model_step_wall, ncol * ntime) do
            SM.step!(domain, step_fields, workspace)
        end

        time_block!(timings, :summarize_columns_cycle) do
            summarize_cycle_columns!(final, domain; backend=:threads)
        end
        time_block!(timings, :cycle_state_deltas) do
            last_delta_thickness_vec .= final.thickness .- prev.thickness
            last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
            last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
            if config.write_outputs
                current_smb_ice_vec = copy(domain.smb_ice)
                last_delta_ice_sheet_smb_vec .= current_smb_ice_vec .- previous_cycle_smb_ice_vec
                previous_cycle_smb_ice_vec .= current_smb_ice_vec
            end
        end

        record = time_block!(timings, :cycle_metrics) do
            (
                cycle = cycle,
                mean_thickness = domain_mean_vector(final.thickness),
                mean_wet_mass = domain_mean_vector(final.wet_mass),
                mean_bulk_density = domain_mean_vector(final.bulk_density),
                mean_base_mass = domain_mean_vector(final.base_mass),
                mean_signed_delta_thickness = domain_mean_vector(last_delta_thickness_vec),
                mean_abs_delta_thickness = domain_mean_abs_vector(last_delta_thickness_vec),
                max_abs_delta_thickness = domain_max_abs_vector(last_delta_thickness_vec),
                mean_signed_delta_wet_mass = domain_mean_vector(last_delta_wet_mass_vec),
                mean_abs_delta_wet_mass = domain_mean_abs_vector(last_delta_wet_mass_vec),
                max_abs_delta_wet_mass = domain_max_abs_vector(last_delta_wet_mass_vec),
                mean_signed_delta_base_mass = domain_mean_vector(last_delta_base_mass_vec),
                mean_abs_delta_base_mass = domain_mean_abs_vector(last_delta_base_mass_vec),
                max_abs_delta_base_mass = domain_max_abs_vector(last_delta_base_mass_vec),
            )
        end
        push!(history, record)

        time_block!(timings, :cycle_logging) do
            println(cycle_log_line(record))
        end

        prev, final = final, prev
    end

    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    final_state = if status == :max_cycles
        history[end].cycle == config.max_cycles ? final : prev
    else
        final
    end

    return (
        history=history,
        status=status,
        final_state=final_state,
        simulation_wall_sec=simulation_wall_sec,
        last_delta_thickness_vec=last_delta_thickness_vec,
        last_delta_wet_mass_vec=last_delta_wet_mass_vec,
        last_delta_base_mass_vec=last_delta_base_mass_vec,
        last_delta_ice_sheet_smb_vec=last_delta_ice_sheet_smb_vec,
    )
end

function write_spinup_history_csv(out_path::AbstractString, history::Vector{NamedTuple})
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, "cycle,mean_thickness_m,mean_wet_mass_mmwe,mean_bulk_density_kgm3,mean_base_mass_mmwe,mean_signed_delta_thickness_m,mean_abs_delta_thickness_m,max_abs_delta_thickness_m,mean_signed_delta_wet_mass_mmwe,mean_abs_delta_wet_mass_mmwe,max_abs_delta_wet_mass_mmwe,mean_signed_delta_base_mass_mmwe,mean_abs_delta_base_mass_mmwe,max_abs_delta_base_mass_mmwe")
        for rec in history
            @printf(
                io,
                "%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
                rec.cycle,
                rec.mean_thickness,
                rec.mean_wet_mass,
                rec.mean_bulk_density,
                rec.mean_base_mass,
                rec.mean_signed_delta_thickness,
                rec.mean_abs_delta_thickness,
                rec.max_abs_delta_thickness,
                rec.mean_signed_delta_wet_mass,
                rec.mean_abs_delta_wet_mass,
                rec.max_abs_delta_wet_mass,
                rec.mean_signed_delta_base_mass,
                rec.mean_abs_delta_base_mass,
                rec.max_abs_delta_base_mass,
            )
        end
    end
end

function write_spinup_summary(
    out_path::AbstractString,
    config,
    time_values::Vector{DateTime},
    nvalid::Int,
    ngrid::Int,
    history::Vector{NamedTuple},
    status::Symbol,
    timings::TimingStats,
)
    last_record = history[end]
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, "Chion GrIS equilibrium spin-up")
        println(io, "Forcing file       : ", abspath(config.nc_path))
        println(io, "Forcing start      : ", first(time_values))
        println(io, "Forcing end        : ", last(time_values))
        println(io, "Forcing steps      : ", length(time_values))
        println(io, "GrIS cells         : ", nvalid, " / ", ngrid)
        println(io, "Backend            : ", String(config.backend))
        println(io, "Threads            : ", nthreads())
        println(io, "File output        : ", config.write_outputs ? "enabled" : "disabled (--no-output)")
        println(io, "NetCDF output      : ", config.write_netcdf ? "enabled" : "disabled (--no-nc)")
        println(io, "Status             : ", string(status))
        println(io, "Cycles completed   : ", length(history))
        println(io)
        println(io, "Final domain means over valid GrIS cells")
        println(io, @sprintf("Thickness (m)              : %.6f", last_record.mean_thickness))
        println(io, @sprintf("Wet mass (mmWE)            : %.6f", last_record.mean_wet_mass))
        println(io, @sprintf("Bulk density (kg m-3)      : %.6f", last_record.mean_bulk_density))
        println(io, @sprintf("Firn-to-ice mass (mmWE)    : %.6f", last_record.mean_base_mass))
        println(io)
        println(io, "Last cycle deltas")
        println(io, @sprintf("Mean signed dThickness (m) : %.6f", last_record.mean_signed_delta_thickness))
        println(io, @sprintf("Mean abs dThickness (m)    : %.6f", last_record.mean_abs_delta_thickness))
        println(io, @sprintf("Max abs dThickness (m)     : %.6f", last_record.max_abs_delta_thickness))
        println(io, @sprintf("Mean signed dSWE (mmWE)    : %.6f", last_record.mean_signed_delta_wet_mass))
        println(io, @sprintf("Mean abs dSWE (mmWE)       : %.6f", last_record.mean_abs_delta_wet_mass))
        println(io, @sprintf("Max abs dSWE (mmWE)        : %.6f", last_record.max_abs_delta_wet_mass))
        println(io, @sprintf("Mean signed dBase (mmWE)   : %.6f", last_record.mean_signed_delta_base_mass))
        println(io, @sprintf("Mean abs dBase (mmWE)      : %.6f", last_record.mean_abs_delta_base_mass))
        println(io, @sprintf("Max abs dBase (mmWE)       : %.6f", last_record.max_abs_delta_base_mass))
        if status == :max_cycles
            println(io)
            println(io, "Interpretation     : Requested cycles completed.")
        else
            println(io)
            println(io, "Interpretation     : Run completed.")
        end
        println(io)
        print_timing_summary(io, timings)
    end
end

function render_spinup_history_plot(out_path::AbstractString, history::Vector{NamedTuple}, config)
    P = plots_module()
    cycles = [rec.cycle for rec in history]
    mean_thickness = [rec.mean_thickness for rec in history]
    mean_wet_mass = [rec.mean_wet_mass for rec in history]
    mean_base_mass = [rec.mean_base_mass for rec in history]
    mean_abs_dth = [rec.mean_abs_delta_thickness for rec in history]
    mean_abs_dswe = [rec.mean_abs_delta_wet_mass for rec in history]
    mean_abs_dbase = [rec.mean_abs_delta_base_mass for rec in history]

    p1 = P.plot(cycles, mean_thickness; lw=3, marker=:circle, color=:steelblue, xlabel="Cycle", ylabel="m", title="Mean snow thickness", framestyle=:box)
    p2 = P.plot(cycles, mean_wet_mass; lw=3, marker=:circle, color=:forestgreen, xlabel="Cycle", ylabel="mmWE", title="Mean snow wet mass", framestyle=:box)
    p3 = P.plot(cycles, mean_base_mass; lw=3, marker=:circle, color=:purple, xlabel="Cycle", ylabel="mmWE", title="Mean firn-to-ice mass", framestyle=:box)
    p4 = P.plot(cycles, mean_abs_dth; lw=3, marker=:circle, color=:firebrick, xlabel="Cycle", ylabel="m", title="Mean abs cycle dThickness", framestyle=:box)
    p5 = P.plot(cycles, mean_abs_dswe; lw=3, marker=:circle, color=:darkorange, xlabel="Cycle", ylabel="mmWE", title="Mean abs cycle dSWE", framestyle=:box)
    p6 = P.plot(cycles, mean_abs_dbase; lw=3, marker=:circle, color=:indigo, xlabel="Cycle", ylabel="mmWE", title="Mean abs cycle dBase", framestyle=:box)

    fig = P.plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1700, 950), plot_title="Chion GrIS spin-up cycle history")
    mkpath(dirname(out_path))
    P.savefig(fig, out_path)
end

function render_spinup_fields_plot(
    out_path::AbstractString,
    x::Vector{Float64},
    y::Vector{Float64},
    initial_thickness::Matrix{Float64},
    final_thickness::Matrix{Float64},
    final_bulk_density::Matrix{Float64},
    final_base_mass::Matrix{Float64},
    last_delta_thickness::Matrix{Float64},
    last_delta_wet_mass::Matrix{Float64},
    last_delta_base_mass::Matrix{Float64},
    final_runoff::Matrix{Float64},
    status::Symbol,
    cycles_completed::Int,
)
    title_suffix = "status=$(string(status)), cycles=$(cycles_completed)"
    P = plots_module()
    fig = P.plot(
        heatmap_panel(x, y, final_thickness; title="Final snow thickness", unit="m", color=:ice),
        heatmap_panel(x, y, final_thickness .- initial_thickness; title="Thickness change", unit="m", clim=symmetric_clims(final_thickness .- initial_thickness), color=P.cgrad([:navy, :white, :firebrick])),
        heatmap_panel(x, y, final_bulk_density; title="Final bulk density", unit="kg/m^3", color=:dense),
        heatmap_panel(x, y, final_base_mass; title="Final firn-to-ice mass", unit="mmWE", color=:amp),
        heatmap_panel(x, y, final_runoff; title="Final cumulative runoff", unit="mmWE", color=:rainbow),
        heatmap_panel(x, y, last_delta_thickness; title="Last cycle dThickness", unit="m", clim=symmetric_clims(last_delta_thickness), color=P.cgrad([:navy, :white, :firebrick])),
        heatmap_panel(x, y, last_delta_wet_mass; title="Last cycle dSWE", unit="mmWE", clim=symmetric_clims(last_delta_wet_mass), color=P.cgrad([:navy, :white, :firebrick])),
        heatmap_panel(x, y, last_delta_base_mass; title="Last cycle dBase", unit="mmWE", clim=symmetric_clims(last_delta_base_mass), color=P.cgrad([:navy, :white, :firebrick])),
        layout=(2, 4),
        size=(2200, 1100),
        plot_title="Chion GrIS equilibrium spin-up fields, " * title_suffix,
    )
    mkpath(dirname(out_path))
    P.savefig(fig, out_path)
end

function collect_final_layer_grids(
    domain::SM.SnowpackDomain,
    js::Vector{Int},
    is::Vector{Int},
    grid_shape::Tuple{Int, Int},
    nlayer::Int,
)
    ny, nx = grid_shape
    n_active = fill(Int32(-1), ny, nx)
    layer_density = fill(NaN, nlayer, ny, nx)
    layer_thickness = fill(NaN, nlayer, ny, nx)
    layer_snow_mass = fill(NaN, nlayer, ny, nx)
    layer_liquid_mass = fill(NaN, nlayer, ny, nx)
    layer_temperature_c = fill(NaN, nlayer, ny, nx)

    @inbounds for idx in 1:SM.column_count(domain)
        j = js[idx]
        i = is[idx]
        n_active[j, i] = Int32(domain.N[idx])
        for k in 1:domain.N[idx]
            rho = domain.density[k, idx]
            m = domain.mass[k, idx]
            mw = domain.mass_w[k, idx]
            layer_density[k, j, i] = rho
            layer_snow_mass[k, j, i] = m
            layer_liquid_mass[k, j, i] = mw
            layer_temperature_c[k, j, i] = domain.temperature[k, idx] - domain.c.T0
            if isfinite(rho) && rho > 0.0 && isfinite(m)
                layer_thickness[k, j, i] = m / rho
            end
        end
    end

    return (
        n_active = n_active,
        layer_density = layer_density,
        layer_thickness = layer_thickness,
        layer_snow_mass = layer_snow_mass,
        layer_liquid_mass = layer_liquid_mass,
        layer_temperature_c = layer_temperature_c,
    )
end

function monthly_vectors_to_grids(values::Matrix{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    nmonth, nvalid = size(values)
    ny, nx = grid_shape
    out = fill(NaN, nmonth, ny, nx)
    @inbounds for m in 1:nmonth, idx in 1:nvalid
        out[m, js[idx], is[idx]] = values[m, idx]
    end
    return out
end

struct SpinupNetCDFWriter
    ncid::Cint
    var_final_th::Cint
    var_final_wet::Cint
    var_final_rho::Cint
    var_final_base::Cint
    var_final_smb_ice::Cint
    var_final_runoff::Cint
    var_dth::Cint
    var_dswe::Cint
    var_dbase::Cint
    var_dsmb_ice::Cint
    var_n::Cint
    var_layer_rho::Cint
    var_layer_th::Cint
    var_layer_mass::Cint
    var_layer_mass_w::Cint
    var_layer_temp::Cint
    var_hist_th::Cint
    var_hist_wet::Cint
    var_hist_rho::Cint
    var_hist_base::Cint
    var_hist_dth::Cint
    var_hist_dswe::Cint
    var_hist_dbase::Cint
    var_monthly_th::Cint
    var_monthly_wet::Cint
    var_monthly_rho::Cint
    var_monthly_base::Cint
    var_monthly_smb_ice::Cint
    var_monthly_export::Cint
    var_monthly_net_ice::Cint
    var_monthly_runoff::Cint
    var_step_valid::Cint
    var_step_export::Cint
    var_step_smb_ice::Cint
    max_steps::Int
    max_cycles::Int
end

function init_spinup_netcdf(
    out_nc::AbstractString,
    config,
    time_values::Vector{DateTime},
    x::Vector{Float64},
    y::Vector{Float64},
    mask::AbstractMatrix{<:Real},
    initial_thickness::Matrix{Float64},
    js::Vector{Int},
    is::Vector{Int},
    month_cycle::Vector{Int32},
    month_of_year::Vector{Int32},
    source_month_code::Vector{Int32},
    annual_output_source_indices::Vector{Int32},
    annual_output_source_codes::Vector{Int32},
)
    mkpath(dirname(out_nc))
    isfile(out_nc) && rm(out_nc, force=true)

    ny, nx = size(mask)
    nlayer = config.ntot
    npoint = length(js)
    max_steps = config.max_cycles * length(annual_output_source_indices)

    step_cycle = Vector{Int32}(undef, max_steps)
    step_source_index = Vector{Int32}(undef, max_steps)
    step_source_code = Vector{Int32}(undef, max_steps)
    step_counter = 0
    for cyc in 1:config.max_cycles
        for annual_idx in eachindex(annual_output_source_indices)
            step_counter += 1
            step_cycle[step_counter] = Int32(cyc)
            step_source_index[step_counter] = annual_output_source_indices[annual_idx]
            step_source_code[step_counter] = annual_output_source_codes[annual_idx]
        end
    end

    ncid = nc_create(out_nc)
    dim_y = nc_def_dim(ncid, "y", ny)
    dim_x = nc_def_dim(ncid, "x", nx)
    dim_layer = nc_def_dim(ncid, "layer", nlayer)
    dim_cycle = nc_def_dim(ncid, "cycle", config.max_cycles)
    dim_month = nc_def_dim(ncid, "month", length(month_cycle))
    dim_point = nc_def_dim(ncid, "point", npoint)
    dim_step = nc_def_dim(ncid, "step", max_steps)

    dims_yx = Cint[dim_y, dim_x]
    dims_lyx = Cint[dim_layer, dim_y, dim_x]
    dims_c = Cint[dim_cycle]
    dims_myx = Cint[dim_month, dim_y, dim_x]
    dims_p = Cint[dim_point]
    dims_syx = Cint[dim_step, dim_y, dim_x]

    var_x = nc_def_var(ncid, "x", NC_DOUBLE, Cint[dim_x])
    nc_put_att_text(ncid, var_x, "units", "km")
    nc_put_att_text(ncid, var_x, "axis", "X")

    var_y = nc_def_var(ncid, "y", NC_DOUBLE, Cint[dim_y])
    nc_put_att_text(ncid, var_y, "units", "km")
    nc_put_att_text(ncid, var_y, "axis", "Y")

    var_layer = nc_def_var(ncid, "layer", NC_INT, Cint[dim_layer])
    nc_put_att_text(ncid, var_layer, "long_name", "Chion internal layer index from surface downward")
    var_cycle = nc_def_var(ncid, "cycle", NC_INT, dims_c)
    nc_put_att_text(ncid, var_cycle, "long_name", "Repeated annual forcing cycle index")
    var_month = nc_def_var(ncid, "month", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_month, "long_name", "Sequential monthly output index")
    var_month_cycle = nc_def_var(ncid, "month_cycle", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_month_cycle, "long_name", "Forcing cycle associated with monthly output")
    var_month_of_year = nc_def_var(ncid, "month_of_year", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_month_of_year, "long_name", "Calendar month of the repeated forcing")
    var_source_month_code = nc_def_var(ncid, "source_month_code", NC_INT, Cint[dim_month])
    nc_put_att_text(ncid, var_source_month_code, "long_name", "Source forcing month code YYYYMM")
    var_point = nc_def_var(ncid, "point", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point, "long_name", "Compact valid GrIS cell index")
    var_point_j = nc_def_var(ncid, "point_j", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point_j, "long_name", "1-based y-index for each compact valid GrIS cell")
    var_point_i = nc_def_var(ncid, "point_i", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point_i, "long_name", "1-based x-index for each compact valid GrIS cell")
    var_point_y = nc_def_var(ncid, "point_y_km", NC_DOUBLE, dims_p)
    nc_put_att_text(ncid, var_point_y, "long_name", "Y coordinate for each compact valid GrIS cell")
    nc_put_att_text(ncid, var_point_y, "units", "km")
    var_point_x = nc_def_var(ncid, "point_x_km", NC_DOUBLE, dims_p)
    nc_put_att_text(ncid, var_point_x, "long_name", "X coordinate for each compact valid GrIS cell")
    nc_put_att_text(ncid, var_point_x, "units", "km")
    var_step = nc_def_var(ncid, "step", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step, "long_name", "Sequential yearly output index across repeated annual cycles")
    var_step_cycle = nc_def_var(ncid, "step_cycle", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_cycle, "long_name", "Repeated annual forcing cycle index for each yearly output")
    var_step_source_index = nc_def_var(ncid, "step_source_index", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_index, "long_name", "1-based index of the last forcing step included in each yearly output")
    var_step_source_code = nc_def_var(ncid, "step_source_code", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_code, "long_name", "Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output")
    var_step_valid = nc_def_var(ncid, "step_valid", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_valid, "long_name", "1 where a yearly output record was completed and written, 0 for unused trailing slots")

    var_mask = define_nc_output_variable(ncid, dims_yx, "gris_mask", "Prepared forcing valid-cell mask", "1")
    var_init_th = define_nc_output_variable(ncid, dims_yx, "initial_thickness", "Initial snow thickness", "m")
    var_final_th = define_nc_output_variable(ncid, dims_yx, "final_thickness", "Final snow thickness", "m")
    var_final_wet = define_nc_output_variable(ncid, dims_yx, "final_wet_mass", "Final snow wet mass", "mmWE")
    var_final_rho = define_nc_output_variable(ncid, dims_yx, "final_bulk_density", "Final bulk snow density", "kg m-3")
    var_final_base = define_nc_output_variable(ncid, dims_yx, "final_base_mass", "Cumulative firn mass exported to the ice model", "mmWE")
    var_final_smb_ice = define_nc_output_variable(ncid, dims_yx, "final_ice_sheet_smb", "Cumulative net mass forcing to the ice sheet", "mmWE")
    var_final_runoff = define_nc_output_variable(ncid, dims_yx, "final_runoff", "Final cumulative runoff", "mmWE")
    var_dth = define_nc_output_variable(ncid, dims_yx, "last_cycle_delta_thickness", "Last cycle snow-thickness change", "m")
    var_dswe = define_nc_output_variable(ncid, dims_yx, "last_cycle_delta_wet_mass", "Last cycle wet-mass change", "mmWE")
    var_dbase = define_nc_output_variable(ncid, dims_yx, "last_cycle_delta_base_mass", "Last cycle firn mass exported to the ice model", "mmWE")
    var_dsmb_ice = define_nc_output_variable(ncid, dims_yx, "last_cycle_delta_ice_sheet_smb", "Last cycle net mass forcing to the ice sheet", "mmWE")
    var_n = nc_def_var(ncid, "n_active", NC_INT, dims_yx)
    nc_put_att_text(ncid, var_n, "long_name", "Number of active Chion layers")

    var_layer_rho = define_nc_output_variable(ncid, dims_lyx, "layer_density", "Final Chion layer density", "kg m-3")
    var_layer_th = define_nc_output_variable(ncid, dims_lyx, "layer_thickness", "Final Chion layer thickness", "m")
    var_layer_mass = define_nc_output_variable(ncid, dims_lyx, "layer_snow_mass", "Final Chion layer snow mass", "kg m-2")
    var_layer_mass_w = define_nc_output_variable(ncid, dims_lyx, "layer_liquid_mass", "Final Chion layer liquid-water mass", "kg m-2")
    var_layer_temp = define_nc_output_variable(ncid, dims_lyx, "layer_temperature_c", "Final Chion layer temperature", "C")

    var_hist_th = define_nc_output_variable(ncid, dims_c, "history_mean_thickness", "Cycle-mean snow thickness", "m")
    var_hist_wet = define_nc_output_variable(ncid, dims_c, "history_mean_wet_mass", "Cycle-mean snow wet mass", "mmWE")
    var_hist_rho = define_nc_output_variable(ncid, dims_c, "history_mean_bulk_density", "Cycle-mean bulk snow density", "kg m-3")
    var_hist_base = define_nc_output_variable(ncid, dims_c, "history_mean_base_mass", "Cycle-mean firn mass exported to the ice model", "mmWE")
    var_hist_dth = define_nc_output_variable(ncid, dims_c, "history_mean_abs_delta_thickness", "Cycle mean absolute snow-thickness change", "m")
    var_hist_dswe = define_nc_output_variable(ncid, dims_c, "history_mean_abs_delta_wet_mass", "Cycle mean absolute wet-mass change", "mmWE")
    var_hist_dbase = define_nc_output_variable(ncid, dims_c, "history_mean_abs_delta_base_mass", "Cycle mean absolute firn mass exported to the ice model", "mmWE")
    var_monthly_th = define_nc_output_variable(ncid, dims_myx, "monthly_mean_thickness", "Monthly mean snow thickness", "m")
    var_monthly_wet = define_nc_output_variable(ncid, dims_myx, "monthly_mean_wet_mass", "Monthly mean snow wet mass", "mmWE")
    var_monthly_rho = define_nc_output_variable(ncid, dims_myx, "monthly_mean_bulk_density", "Monthly mean bulk snow density", "kg m-3")
    var_monthly_base = define_nc_output_variable(ncid, dims_myx, "monthly_mean_base_mass", "Monthly mean cumulative firn mass exported to the ice model", "mmWE")
    var_monthly_smb_ice = define_nc_output_variable(ncid, dims_myx, "monthly_mean_ice_sheet_smb", "Monthly net mass forcing to the ice sheet", "mmWE")
    var_monthly_export = define_nc_output_variable(ncid, dims_myx, "monthly_export_to_ice", "Monthly firn mass exported to the ice model", "mmWE")
    var_monthly_net_ice = define_nc_output_variable(ncid, dims_myx, "monthly_net_ice_sheet_forcing", "Monthly net mass forcing to the ice sheet", "mmWE")
    var_monthly_runoff = define_nc_output_variable(ncid, dims_myx, "monthly_runoff", "Monthly runoff production", "mmWE")
    var_step_export = define_nc_output_variable(ncid, dims_syx, "step_export_to_ice", "Annual firn mass exported to the ice model for each written output interval", "mmWE")
    var_step_smb_ice = define_nc_output_variable(ncid, dims_syx, "step_ice_sheet_smb", "Annual net mass forcing to the ice sheet for each written output interval", "mmWE")

    nc_put_att_text(ncid, NC_GLOBAL, "title", "Chion GrIS spin-up final state")
    nc_put_att_text(ncid, NC_GLOBAL, "source_model", "Chion")
    nc_put_att_text(ncid, NC_GLOBAL, "history", "Created by examples/scripts/run_gris_equilibrium.jl")
    nc_put_att_text(ncid, NC_GLOBAL, "forcing_start", string(first(time_values)))
    nc_put_att_text(ncid, NC_GLOBAL, "forcing_end", string(last(time_values)))
    nc_put_att_text(ncid, NC_GLOBAL, "cycles_completed", "pending")
    nc_put_att_text(ncid, NC_GLOBAL, "status", "pending")
    nc_put_att_text(ncid, NC_GLOBAL, "layer_note", "Layer dimension is the internal Chion layer index; inactive layers are stored as NaN.")
    nc_put_att_text(ncid, NC_GLOBAL, "firn_export_note", "final_base_mass is the cumulative firn mass transferred to the ice model proxy; last_cycle_delta_base_mass is the export during the final annual cycle.")
    nc_put_att_text(ncid, NC_GLOBAL, "ice_sheet_smb_note", "final_ice_sheet_smb is the cumulative net mass forcing to the ice sheet: positive firn export minus bare-ice ablation.")
    nc_put_att_text(ncid, NC_GLOBAL, "monthly_note", "Monthly mean state fields are averages over Chion daily states; monthly forcing fields are integrated monthly totals.")
    nc_put_att_text(ncid, NC_GLOBAL, "step_note", "step_export_to_ice and step_ice_sheet_smb are yearly accumulated fields written once per source year in each repeated forcing cycle.")
    nc_put_att_text(ncid, NC_GLOBAL, "created", string(now()))

    nc_enddef(ncid)

    nc_put_var_double(ncid, var_x, x)
    nc_put_var_double(ncid, var_y, y)
    nc_put_var_int(ncid, var_layer, Int32.(collect(1:nlayer)))
    nc_put_var_int(ncid, var_cycle, Int32.(collect(1:config.max_cycles)))
    nc_put_var_int(ncid, var_month, Int32.(collect(1:length(month_cycle))))
    nc_put_var_int(ncid, var_month_cycle, month_cycle)
    nc_put_var_int(ncid, var_month_of_year, month_of_year)
    nc_put_var_int(ncid, var_source_month_code, source_month_code)
    nc_put_var_int(ncid, var_point, Int32.(collect(1:npoint)))
    nc_put_var_int(ncid, var_point_j, Int32.(js))
    nc_put_var_int(ncid, var_point_i, Int32.(is))
    nc_put_var_double(ncid, var_point_y, y[js])
    nc_put_var_double(ncid, var_point_x, x[is])
    nc_put_var_int(ncid, var_step, Int32.(collect(1:max_steps)))
    nc_put_var_int(ncid, var_step_cycle, step_cycle)
    nc_put_var_int(ncid, var_step_source_index, step_source_index)
    nc_put_var_int(ncid, var_step_source_code, step_source_code)
    nc_put_var_float_2d(ncid, var_mask, mask)
    nc_put_var_float_2d(ncid, var_init_th, initial_thickness)

    return SpinupNetCDFWriter(
        ncid,
        var_final_th,
        var_final_wet,
        var_final_rho,
        var_final_base,
        var_final_smb_ice,
        var_final_runoff,
        var_dth,
        var_dswe,
        var_dbase,
        var_dsmb_ice,
        var_n,
        var_layer_rho,
        var_layer_th,
        var_layer_mass,
        var_layer_mass_w,
        var_layer_temp,
        var_hist_th,
        var_hist_wet,
        var_hist_rho,
        var_hist_base,
        var_hist_dth,
        var_hist_dswe,
        var_hist_dbase,
        var_monthly_th,
        var_monthly_wet,
        var_monthly_rho,
        var_monthly_base,
        var_monthly_smb_ice,
        var_monthly_export,
        var_monthly_net_ice,
        var_monthly_runoff,
        var_step_valid,
        var_step_export,
        var_step_smb_ice,
        max_steps,
        config.max_cycles,
    )
end

function write_step_export_to_ice!(writer::SpinupNetCDFWriter, step_index::Int, step_export_to_ice::AbstractMatrix{<:Real})
    nc_put_vara_float_3d_step_yx(writer.ncid, writer.var_step_export, step_index, step_export_to_ice)
    return
end

function write_step_ice_sheet_smb!(writer::SpinupNetCDFWriter, step_index::Int, step_ice_sheet_smb::AbstractMatrix{<:Real})
    nc_put_vara_float_3d_step_yx(writer.ncid, writer.var_step_smb_ice, step_index, step_ice_sheet_smb)
    return
end

function finalize_spinup_netcdf!(
    writer::SpinupNetCDFWriter,
    final_thickness::Matrix{Float64},
    final_wet_mass::Matrix{Float64},
    final_bulk_density::Matrix{Float64},
    final_base_mass::Matrix{Float64},
    final_ice_sheet_smb::Matrix{Float64},
    last_delta_thickness::Matrix{Float64},
    last_delta_wet_mass::Matrix{Float64},
    last_delta_base_mass::Matrix{Float64},
    last_delta_ice_sheet_smb::Matrix{Float64},
    final_runoff::Matrix{Float64},
    layer_grids,
    history::Vector{NamedTuple},
    monthly_mean_thickness::Array{Float64, 3},
    monthly_mean_wet_mass::Array{Float64, 3},
    monthly_mean_bulk_density::Array{Float64, 3},
    monthly_mean_base_mass::Array{Float64, 3},
    monthly_mean_ice_sheet_smb::Array{Float64, 3},
    monthly_export_to_ice::Array{Float64, 3},
    monthly_net_ice_sheet_forcing::Array{Float64, 3},
    monthly_runoff::Array{Float64, 3},
    status::Symbol,
    cycles_completed::Int,
    steps_written::Int,
)
    step_valid = zeros(Int32, writer.max_steps)
    step_valid[1:steps_written] .= 1

    hist_th = fill(NaN, writer.max_cycles)
    hist_wet = fill(NaN, writer.max_cycles)
    hist_rho = fill(NaN, writer.max_cycles)
    hist_base = fill(NaN, writer.max_cycles)
    hist_dth = fill(NaN, writer.max_cycles)
    hist_dswe = fill(NaN, writer.max_cycles)
    hist_dbase = fill(NaN, writer.max_cycles)
    for rec in history
        idx = getfield(rec, :cycle)
        hist_th[idx] = getfield(rec, :mean_thickness)
        hist_wet[idx] = getfield(rec, :mean_wet_mass)
        hist_rho[idx] = getfield(rec, :mean_bulk_density)
        hist_base[idx] = getfield(rec, :mean_base_mass)
        hist_dth[idx] = getfield(rec, :mean_abs_delta_thickness)
        hist_dswe[idx] = getfield(rec, :mean_abs_delta_wet_mass)
        hist_dbase[idx] = getfield(rec, :mean_abs_delta_base_mass)
    end

    nc_put_var_float_2d(writer.ncid, writer.var_final_th, final_thickness)
    nc_put_var_float_2d(writer.ncid, writer.var_final_wet, final_wet_mass)
    nc_put_var_float_2d(writer.ncid, writer.var_final_rho, final_bulk_density)
    nc_put_var_float_2d(writer.ncid, writer.var_final_base, final_base_mass)
    nc_put_var_float_2d(writer.ncid, writer.var_final_smb_ice, final_ice_sheet_smb)
    nc_put_var_float_2d(writer.ncid, writer.var_final_runoff, final_runoff)
    nc_put_var_float_2d(writer.ncid, writer.var_dth, last_delta_thickness)
    nc_put_var_float_2d(writer.ncid, writer.var_dswe, last_delta_wet_mass)
    nc_put_var_float_2d(writer.ncid, writer.var_dbase, last_delta_base_mass)
    nc_put_var_float_2d(writer.ncid, writer.var_dsmb_ice, last_delta_ice_sheet_smb)
    nc_put_var_int_2d(writer.ncid, writer.var_n, layer_grids.n_active)
    nc_put_var_float_3d(writer.ncid, writer.var_layer_rho, layer_grids.layer_density)
    nc_put_var_float_3d(writer.ncid, writer.var_layer_th, layer_grids.layer_thickness)
    nc_put_var_float_3d(writer.ncid, writer.var_layer_mass, layer_grids.layer_snow_mass)
    nc_put_var_float_3d(writer.ncid, writer.var_layer_mass_w, layer_grids.layer_liquid_mass)
    nc_put_var_float_3d(writer.ncid, writer.var_layer_temp, layer_grids.layer_temperature_c)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_th, hist_th)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_wet, hist_wet)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_rho, hist_rho)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_base, hist_base)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_dth, hist_dth)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_dswe, hist_dswe)
    nc_put_var_float_1d(writer.ncid, writer.var_hist_dbase, hist_dbase)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_th, monthly_mean_thickness)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_wet, monthly_mean_wet_mass)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_rho, monthly_mean_bulk_density)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_base, monthly_mean_base_mass)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_smb_ice, monthly_mean_ice_sheet_smb)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_export, monthly_export_to_ice)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_net_ice, monthly_net_ice_sheet_forcing)
    nc_put_var_float_3d(writer.ncid, writer.var_monthly_runoff, monthly_runoff)
    nc_put_var_int(writer.ncid, writer.var_step_valid, step_valid)

    nc_redef(writer.ncid)
    nc_put_att_text(writer.ncid, NC_GLOBAL, "cycles_completed", string(cycles_completed))
    nc_put_att_text(writer.ncid, NC_GLOBAL, "status", string(status))
    nc_put_att_text(writer.ncid, NC_GLOBAL, "steps_written", string(steps_written))
    nc_enddef(writer.ncid)
    nc_close(writer.ncid)
    return
end

function main(args::Vector{String})
    if has_flag(args, "help")
        print_spinup_help()
        return
    end

    ensure_tools!()
    config = parse_spinup_config(args)
    timings = TimingStats()
    run_wall_t0 = time_ns()
    needs_spatial_output = config.write_outputs || config.write_netcdf
    shapes = time_block!(timings, :read_dataset_shapes) do
        read_dataset_shapes(config.nc_path)
    end
    date_codes, time_values = time_block!(timings, :read_forcing_times) do
        read_forcing_times(config.nc_path, shapes)
    end
    dt_days = time_block!(timings, :prepare_timestep_sizes) do
        [infer_dt_days(time_values, t) for t in eachindex(time_values)]
    end

    x = needs_spatial_output ? time_block!(timings, :read_grid) do
        vec(read_hdf5_full(config.nc_path, "x", shapes))
    end : Float64[]
    y = needs_spatial_output ? time_block!(timings, :read_grid) do
        vec(read_hdf5_full(config.nc_path, "y", shapes))
    end : Float64[]
    mask = time_block!(timings, :read_mask) do
        read_hdf5_full(config.nc_path, "MSK", shapes)
    end
    outlay_bounds = time_block!(timings, :read_outlay_bounds) do
        read_hdf5_full(config.nc_path, "OUTLAY_bnds", shapes)
    end

    tt_full = time_block!(timings, :read_forcing_tt) do
        read_full_timeseries_3d(config.nc_path, "TT", shapes)
    end
    sf_full = time_block!(timings, :read_forcing_sf) do
        read_full_timeseries_3d(config.nc_path, "SF", shapes)
    end
    rf_full = time_block!(timings, :read_forcing_rf) do
        read_full_timeseries_3d(config.nc_path, "RF", shapes)
    end
    swd_full = time_block!(timings, :read_forcing_swd) do
        read_full_timeseries_3d(config.nc_path, "SWD", shapes)
    end
    lwd_full = time_block!(timings, :read_forcing_lwd) do
        read_full_timeseries_3d(config.nc_path, "LWD", shapes)
    end
    shf_full = time_block!(timings, :read_forcing_shf) do
        read_full_timeseries_3d(config.nc_path, "SHF", shapes)
    end
    lhf_full = time_block!(timings, :read_forcing_lhf) do
        read_full_timeseries_3d(config.nc_path, "LHF", shapes)
    end
    u_wind_info = time_block!(timings, :read_forcing_wind) do
        read_first_available_timeseries_3d(config.nc_path, ["UU", "U10"], shapes)
    end
    v_wind_info = time_block!(timings, :read_forcing_wind) do
        read_first_available_timeseries_3d(config.nc_path, ["VV", "V10"], shapes)
    end
    wind_full, wind_forcing_message = time_block!(timings, :prepare_wind_forcing) do
        wind_full_local = if !isnothing(u_wind_info) && !isnothing(v_wind_info)
            hypot.(u_wind_info.data, v_wind_info.data)
        else
            nothing
        end
        wind_message = if isnothing(wind_full_local)
            "Wind forcing: file wind components not found; using default 5.0 m s^-1."
        else
            @sprintf(
                "Wind forcing: |V| from components %s and %s.",
                u_wind_info.name,
                v_wind_info.name,
            )
        end
        return wind_full_local, wind_message
    end
    println(wind_forcing_message)

    zn3_init = time_block!(timings, :read_initial_state) do
        read_timeslice_2d(config.nc_path, "ZN3", 1, shapes)
    end
    ro1_init = time_block!(timings, :read_initial_state) do
        read_timeslice_3d(config.nc_path, "RO1", 1, shapes)
    end
    ti1_init = time_block!(timings, :read_initial_state) do
        read_timeslice_3d(config.nc_path, "TI1", 1, shapes)
    end
    wa1_init = time_block!(timings, :read_initial_state) do
        read_timeslice_3d(config.nc_path, "WA1", 1, shapes)
    end

    ny, nx = size(mask)
    valid_mask, valid_indices, nvalid = time_block!(timings, :build_valid_domain) do
        valid_mask_local = falses(ny, nx)
        @inbounds for j in 1:ny, i in 1:nx
            valid_mask_local[j, i] =
                all(isfinite, @view(tt_full[:, j, i])) &&
                all(isfinite, @view(sf_full[:, j, i])) &&
                all(isfinite, @view(rf_full[:, j, i])) &&
                all(isfinite, @view(swd_full[:, j, i]))
        end
        valid_indices_local = findall(valid_mask_local)
        return valid_mask_local, valid_indices_local, length(valid_indices_local)
    end
    ntime = length(time_values)
    annual_output, unique_month_keys, step_month, nmonth_per_cycle, nmonth_total, month_cycle, month_of_year, source_month_code =
        time_block!(timings, :prepare_output_schedule) do
            annual_output_local = config.write_netcdf ? build_annual_output_schedule(time_values) : nothing
            unique_month_keys_local = config.write_netcdf ? unique((year(t), month(t)) for t in time_values) : Tuple{Int, Int}[]
            month_lookup = Dict{Tuple{Int, Int}, Int}()
            if config.write_netcdf
                for (idx, key) in enumerate(unique_month_keys_local)
                    month_lookup[key] = idx
                end
            end
            step_month_local = config.write_netcdf ? [month_lookup[(year(t), month(t))] for t in time_values] : Int[]
            nmonth_per_cycle_local = config.write_netcdf ? length(unique_month_keys_local) : 0
            nmonth_total_local = config.write_netcdf ? config.max_cycles * nmonth_per_cycle_local : 0
            month_cycle_local = Int32[]
            month_of_year_local = Int32[]
            source_month_code_local = Int32[]
            if config.write_netcdf
                for cyc in 1:config.max_cycles, key in unique_month_keys_local
                    push!(month_cycle_local, Int32(cyc))
                    push!(month_of_year_local, Int32(key[2]))
                    push!(source_month_code_local, Int32(key[1] * 100 + key[2]))
                end
            end
            return (
                annual_output_local,
                unique_month_keys_local,
                step_month_local,
                nmonth_per_cycle_local,
                nmonth_total_local,
                month_cycle_local,
                month_of_year_local,
                source_month_code_local,
            )
        end

    js, is, domain, initial_thickness_vec, tair_k, snow_rate, rain_rate, s_boa, q_lw, has_q_lw, q_sh, has_q_sh, q_lh, has_q_lh, wind_speed =
        time_block!(timings, :allocate_simulation_arrays) do
            js_local = needs_spatial_output ? Vector{Int}(undef, nvalid) : Int[]
            is_local = needs_spatial_output ? Vector{Int}(undef, nvalid) : Int[]
            domain_local = SM.SnowpackDomain(ncol=nvalid, Ntot=config.ntot)
            initial_thickness_vec_local = config.write_outputs ? Vector{Float64}(undef, nvalid) : Float64[]

            tair_k_local = Matrix{Float64}(undef, nvalid, ntime)
            snow_rate_local = Matrix{Float64}(undef, nvalid, ntime)
            rain_rate_local = Matrix{Float64}(undef, nvalid, ntime)
            s_boa_local = Matrix{Float64}(undef, nvalid, ntime)
            q_lw_local = Matrix{Float64}(undef, nvalid, ntime)
            has_q_lw_local = fill(false, nvalid, ntime)
            q_sh_local = Matrix{Float64}(undef, nvalid, ntime)
            has_q_sh_local = fill(false, nvalid, ntime)
            q_lh_local = Matrix{Float64}(undef, nvalid, ntime)
            has_q_lh_local = fill(false, nvalid, ntime)
            wind_speed_local = Matrix{Float64}(undef, nvalid, ntime)
            return (
                js_local,
                is_local,
                domain_local,
                initial_thickness_vec_local,
                tair_k_local,
                snow_rate_local,
                rain_rate_local,
                s_boa_local,
                q_lw_local,
                has_q_lw_local,
                q_sh_local,
                has_q_sh_local,
                q_lh_local,
                has_q_lh_local,
                wind_speed_local,
            )
        end

    time_block!(timings, :initialize_columns_and_forcing) do
        @threads :static for idx in eachindex(valid_indices)
            j, i = Tuple(valid_indices[idx])
            if needs_spatial_output
                js[idx] = j
                is[idx] = i
            end

            populate_domain_column_from_forcing_file!(
                domain,
                idx,
                Float64(zn3_init[j, i]),
                @view(ro1_init[:, j, i]),
                @view(ti1_init[:, j, i]),
                @view(wa1_init[:, j, i]),
                outlay_bounds,
            )

            if config.write_outputs
                init_summary = summarize_column(domain, idx)
                initial_thickness_vec[idx] = init_summary.thickness
            end

            for t in 1:ntime
                tair_k[idx, t] = valid_or(-15.0, Float64(tt_full[t, j, i])) + domain.c.T0
                snow_rate[idx, t] = mmwe_day_to_kgm2s(Float64(sf_full[t, j, i]))
                rain_rate[idx, t] = mmwe_day_to_kgm2s(Float64(rf_full[t, j, i]))
                s_boa[idx, t] = valid_or(0.0, Float64(swd_full[t, j, i]))
                q_lw_ij = Float64(lwd_full[t, j, i])
                q_sh_ij = Float64(shf_full[t, j, i])
                q_lh_ij = Float64(lhf_full[t, j, i])
                has_q_lw[idx, t] = isfinite(q_lw_ij)
                has_q_sh[idx, t] = isfinite(q_sh_ij)
                has_q_lh[idx, t] = isfinite(q_lh_ij)
                q_lw[idx, t] = has_q_lw[idx, t] ? q_lw_ij : 0.0
                q_sh[idx, t] = has_q_sh[idx, t] ? q_sh_ij : 0.0
                q_lh[idx, t] = has_q_lh[idx, t] ? q_lh_ij : 0.0
                wind_speed[idx, t] = isnothing(wind_full) ? 5.0 : valid_or(5.0, Float64(wind_full[t, j, i]))
            end
        end
    end

    step_fields = SM.SnowpackForcing(
        dt_days=dt_days,
        air_temperature=tair_k,
        snowfall_rate=snow_rate,
        rainfall_rate=rain_rate,
        shortwave_down=s_boa,
        wind_speed=wind_speed,
        q_lw_down=q_lw,
        has_q_lw_down=has_q_lw,
        q_sh=q_sh,
        has_q_sh=has_q_sh,
        q_lh=q_lh,
        has_q_lh=has_q_lh,
        time_values=time_values,
    )

    if config.backend == :threads && !config.write_outputs && !config.write_netcdf
        workspace = time_block!(timings, :create_workspaces) do
            SM.ColumnarStepWorkspace(domain)
        end
        sim = run_spinup_cycles_no_netcdf!(
            timings,
            config,
            domain,
            workspace,
            step_fields,
            nvalid,
        )
        run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9

        println("GrIS equilibrium spin-up complete.")
        println("Forcing file   : $(abspath(config.nc_path))")
        println("Forcing start  : $(first(time_values))")
        println("Forcing end    : $(last(time_values))")
        println("Backend        : $(String(config.backend))")
        println("Cycles         : $(length(sim.history))")
        println("Status         : $(string(sim.status))")
        println(@sprintf("Simulation wall: %.3f s", sim.simulation_wall_sec))
        println(@sprintf("Run wall total : %.3f s", run_wall_sec))
        if haskey(timings.totals, :model_step_wall)
            println(@sprintf("Model step wall: %.3f s", timings.totals[:model_step_wall]))
        end
        println("Output NetCDF  : skipped (--no-nc)")
        println("File outputs   : skipped (--no-output)")
        print_timing_summary(stdout, timings; total_wall_sec=run_wall_sec)
        return
    end

    if config.backend == :gpu
        SM.cuda_available() || error("`--backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            SM.gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            SM.adapt(CUDA.CuArray, step_fields)
        end
        workspaces = time_block!(timings, :gpu_transfer) do
            SM.ColumnarStepWorkspace(domain)
        end
    else
        workspaces = time_block!(timings, :create_workspaces) do
            SM.ColumnarStepWorkspace(domain)
        end
    end

    gpu_netcdf_diagnostics = config.write_netcdf && config.backend == :gpu
    ncol = SM.column_count(domain)
    prev, final, last_delta_thickness_vec, last_delta_wet_mass_vec, last_delta_base_mass_vec, last_delta_ice_sheet_smb_vec =
        time_block!(timings, :allocate_cycle_buffers) do
            prev_local = config.backend == :gpu ? allocate_cycle_summary_buffers(domain, nvalid) : allocate_cycle_summary_buffers(nvalid)
            final_local = config.backend == :gpu ? allocate_cycle_summary_buffers(domain, nvalid) : allocate_cycle_summary_buffers(nvalid)
            last_delta_thickness_vec_local = config.backend == :gpu ? similar(domain.mass, Float64, nvalid) : fill(NaN, nvalid)
            last_delta_wet_mass_vec_local = config.backend == :gpu ? similar(domain.mass, Float64, nvalid) : fill(NaN, nvalid)
            last_delta_base_mass_vec_local = config.backend == :gpu ? similar(domain.mass, Float64, nvalid) : fill(NaN, nvalid)
            last_delta_ice_sheet_smb_vec_local = config.backend == :gpu ? similar(domain.mass, Float64, nvalid) : fill(NaN, nvalid)
            return (
                prev_local,
                final_local,
                last_delta_thickness_vec_local,
                last_delta_wet_mass_vec_local,
                last_delta_base_mass_vec_local,
                last_delta_ice_sheet_smb_vec_local,
            )
        end
    time_block!(timings, :summarize_columns_initial) do
        if config.backend == :gpu
            SM.summarize_cycle_state!(
                prev.thickness,
                prev.wet_mass,
                prev.bulk_density,
                prev.base_mass,
                domain
            )
        else
            summarize_cycle_columns!(prev, domain; backend=:threads)
        end
    end
    previous_base_mass_vec, previous_smb_ice_vec, previous_runoff_vec, previous_cycle_smb_ice_vec =
        time_block!(timings, :initialize_cycle_tracking) do
            previous_base_mass_vec_local = if gpu_netcdf_diagnostics
                copy(prev.base_mass)
            elseif config.write_netcdf
                copy(prev.base_mass)
            else
                Float64[]
            end
            previous_smb_ice_vec_local = if gpu_netcdf_diagnostics
                copy(domain.smb_ice)
            elseif config.write_netcdf
                copy(domain.smb_ice)
            else
                Float64[]
            end
            previous_runoff_vec_local = if gpu_netcdf_diagnostics
                copy(domain.runoff)
            elseif config.write_netcdf
                copy(domain.runoff)
            else
                Float64[]
            end
            previous_cycle_smb_ice_vec_local = if config.write_outputs || config.write_netcdf
                copy(domain.smb_ice)
            else
                Float64[]
            end
            return (
                previous_base_mass_vec_local,
                previous_smb_ice_vec_local,
                previous_runoff_vec_local,
                previous_cycle_smb_ice_vec_local,
            )
        end
    initial_thickness = config.write_outputs ? time_block!(timings, :prepare_initial_output_fields) do
        scatter_to_grid(initial_thickness_vec, js, is, (ny, nx))
    end : Matrix{Float64}(undef, 0, 0)
    nc_path = isempty(config.out_nc) ? joinpath(config.out_dir, "gris_equilibrium_final_state.nc") : config.out_nc
    writer = config.write_netcdf ? time_block!(timings, :init_netcdf) do
        init_spinup_netcdf(
            nc_path,
            config,
            time_values,
            x,
            y,
            mask,
            initial_thickness,
            js,
            is,
            month_cycle,
            month_of_year,
            source_month_code,
            annual_output.source_indices,
            annual_output.source_codes,
        )
    end : nothing
    step_summary, device_step_summary, annual_export_to_ice, annual_ice_sheet_smb, steps_written, history, status,
    monthly_sum_thickness, monthly_sum_wet_mass, monthly_sum_bulk_density, monthly_sum_base_mass,
    monthly_sum_ice_sheet_smb, monthly_sum_export, monthly_sum_net_ice_sheet_forcing, monthly_sum_runoff, monthly_count =
        time_block!(timings, :allocate_output_buffers) do
            step_summary_local = config.write_netcdf && !gpu_netcdf_diagnostics ? allocate_summary_buffers(nvalid) : nothing
            device_step_summary_local = gpu_netcdf_diagnostics ? allocate_summary_buffers(domain, nvalid) : nothing
            annual_export_to_ice_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nvalid)
            else
                Float64[]
            end
            annual_ice_sheet_smb_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nvalid)
            else
                Float64[]
            end
            history_local = NamedTuple[]
            status_local = :max_cycles
            monthly_sum_thickness_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_wet_mass_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_bulk_density_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_base_mass_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_ice_sheet_smb_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_export_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_net_ice_sheet_forcing_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_runoff_local = if gpu_netcdf_diagnostics
                CUDA.zeros(Float64, nmonth_total, nvalid)
            elseif config.write_netcdf
                zeros(Float64, nmonth_total, nvalid)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_count_local = config.write_netcdf ? zeros(Int32, nmonth_total) : Int32[]
            return (
                step_summary_local,
                device_step_summary_local,
                annual_export_to_ice_local,
                annual_ice_sheet_smb_local,
                0,
                history_local,
                status_local,
                monthly_sum_thickness_local,
                monthly_sum_wet_mass_local,
                monthly_sum_bulk_density_local,
                monthly_sum_base_mass_local,
                monthly_sum_ice_sheet_smb_local,
                monthly_sum_export_local,
                monthly_sum_net_ice_sheet_forcing_local,
                monthly_sum_runoff_local,
                monthly_count_local,
            )
        end
    simulation_wall_t0 = time_ns()

    for cycle in 1:config.max_cycles
        for t in 1:ntime
            month_idx = config.write_netcdf ? (cycle - 1) * nmonth_per_cycle + step_month[t] : 0
            if config.backend == :gpu
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)
                if config.write_netcdf
                    t1 = time_ns()
                    if gpu_netcdf_diagnostics
                        SM.summarize_domain_state!(
                            device_step_summary.thickness,
                            device_step_summary.wet_mass,
                            device_step_summary.bulk_density,
                            device_step_summary.base_mass,
                            device_step_summary.smb_ice,
                            device_step_summary.liquid_water,
                            device_step_summary.runoff,
                            domain
            )
                        current_base = device_step_summary.base_mass
                        current_smb_ice = device_step_summary.smb_ice
                        current_runoff = device_step_summary.runoff
                        monthly_sum_thickness[month_idx, :] .+= device_step_summary.thickness
                        monthly_sum_wet_mass[month_idx, :] .+= device_step_summary.wet_mass
                        monthly_sum_bulk_density[month_idx, :] .+= device_step_summary.bulk_density
                        monthly_sum_base_mass[month_idx, :] .+= current_base
                        monthly_sum_ice_sheet_smb[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                        monthly_sum_export[month_idx, :] .+= current_base .- previous_base_mass_vec
                        monthly_sum_net_ice_sheet_forcing[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                        monthly_sum_runoff[month_idx, :] .+= current_runoff .- previous_runoff_vec
                        annual_export_to_ice .+= current_base .- previous_base_mass_vec
                        annual_ice_sheet_smb .+= current_smb_ice .- previous_smb_ice_vec
                        previous_base_mass_vec .= current_base
                        previous_smb_ice_vec .= current_smb_ice
                        previous_runoff_vec .= current_runoff
                    else
                        summarize_columns!(step_summary, domain; backend=:kernelabstractions, device_summary=device_step_summary)
                        current_base = step_summary.base_mass
                        current_smb_ice = step_summary.smb_ice
                        current_runoff = step_summary.runoff
                        monthly_sum_thickness[month_idx, :] .+= step_summary.thickness
                        monthly_sum_wet_mass[month_idx, :] .+= step_summary.wet_mass
                        monthly_sum_bulk_density[month_idx, :] .+= step_summary.bulk_density
                        monthly_sum_base_mass[month_idx, :] .+= current_base
                        monthly_sum_ice_sheet_smb[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                        monthly_sum_export[month_idx, :] .+= current_base .- previous_base_mass_vec
                        monthly_sum_net_ice_sheet_forcing[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                        monthly_sum_runoff[month_idx, :] .+= current_runoff .- previous_runoff_vec
                        annual_export_to_ice .+= current_base .- previous_base_mass_vec
                        annual_ice_sheet_smb .+= current_smb_ice .- previous_smb_ice_vec
                        previous_base_mass_vec .= current_base
                        previous_smb_ice_vec .= current_smb_ice
                        previous_runoff_vec .= current_runoff
                    end
                    add_timing!(timings, :step_diagnostics, (time_ns() - t1) * 1.0e-9, ncol)
                end
            elseif config.write_netcdf
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                @threads :static for idx in 1:ncol
                    summary = summarize_column(domain, idx)
                    monthly_sum_thickness[month_idx, idx] += summary.thickness
                    monthly_sum_wet_mass[month_idx, idx] += summary.wet_mass
                    monthly_sum_bulk_density[month_idx, idx] += summary.bulk_density
                    monthly_sum_base_mass[month_idx, idx] += domain.mass_base[idx]
                    monthly_sum_ice_sheet_smb[month_idx, idx] += domain.smb_ice[idx] - previous_smb_ice_vec[idx]
                    monthly_sum_export[month_idx, idx] += domain.mass_base[idx] - previous_base_mass_vec[idx]
                    monthly_sum_net_ice_sheet_forcing[month_idx, idx] += domain.smb_ice[idx] - previous_smb_ice_vec[idx]
                    monthly_sum_runoff[month_idx, idx] += domain.runoff[idx] - previous_runoff_vec[idx]
                    annual_export_to_ice[idx] += domain.mass_base[idx] - previous_base_mass_vec[idx]
                    annual_ice_sheet_smb[idx] += domain.smb_ice[idx] - previous_smb_ice_vec[idx]
                    previous_base_mass_vec[idx] = domain.mass_base[idx]
                    previous_smb_ice_vec[idx] = domain.smb_ice[idx]
                    previous_runoff_vec[idx] = domain.runoff[idx]
                end
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)
            else
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)
            end
            if config.write_netcdf
                monthly_count[month_idx] += 1
                if annual_output.write_output[t]
                    steps_written += 1
                    step_export_to_ice_vec, step_ice_sheet_smb_vec, step_export_to_ice, step_ice_sheet_smb =
                        time_block!(timings, :step_output_prepare) do
                            step_export_to_ice_vec_local = gpu_netcdf_diagnostics ? Array(annual_export_to_ice) : annual_export_to_ice
                            step_ice_sheet_smb_vec_local = gpu_netcdf_diagnostics ? Array(annual_ice_sheet_smb) : annual_ice_sheet_smb
                            step_export_to_ice_local = scatter_to_grid(step_export_to_ice_vec_local, js, is, (ny, nx))
                            step_ice_sheet_smb_local = scatter_to_grid(step_ice_sheet_smb_vec_local, js, is, (ny, nx))
                            return (
                                step_export_to_ice_vec_local,
                                step_ice_sheet_smb_vec_local,
                                step_export_to_ice_local,
                                step_ice_sheet_smb_local,
                            )
                        end
                    time_block!(timings, :step_output_write) do
                        write_step_export_to_ice!(writer, steps_written, step_export_to_ice)
                        write_step_ice_sheet_smb!(writer, steps_written, step_ice_sheet_smb)
                    end
                    fill!(annual_export_to_ice, 0.0)
                    fill!(annual_ice_sheet_smb, 0.0)
                end
            end
        end

        record = if config.backend == :gpu
            time_block!(timings, :summarize_columns_cycle) do
                SM.summarize_cycle_state!(
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    domain
            )
            end
            time_block!(timings, :cycle_state_deltas) do
                last_delta_thickness_vec .= final.thickness .- prev.thickness
                last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
                last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
                if config.write_outputs || config.write_netcdf
                    last_delta_ice_sheet_smb_vec .= domain.smb_ice .- previous_cycle_smb_ice_vec
                    previous_cycle_smb_ice_vec .= domain.smb_ice
                end
            end
            time_block!(timings, :cycle_metrics) do
                make_cycle_record(
                    cycle,
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    last_delta_thickness_vec,
                    last_delta_wet_mass_vec,
                    last_delta_base_mass_vec,
                )
            end
        else
            time_block!(timings, :summarize_columns_cycle) do
                summarize_cycle_columns!(final, domain; backend=:threads)
            end
            time_block!(timings, :cycle_state_deltas) do
                last_delta_thickness_vec .= final.thickness .- prev.thickness
                last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
                last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
                if config.write_outputs || config.write_netcdf
                    current_smb_ice_vec = copy(domain.smb_ice)
                    last_delta_ice_sheet_smb_vec .= current_smb_ice_vec .- previous_cycle_smb_ice_vec
                    previous_cycle_smb_ice_vec .= current_smb_ice_vec
                end
            end
            time_block!(timings, :cycle_metrics) do
                (
                    cycle = cycle,
                    mean_thickness = domain_mean_vector(final.thickness),
                    mean_wet_mass = domain_mean_vector(final.wet_mass),
                    mean_bulk_density = domain_mean_vector(final.bulk_density),
                    mean_base_mass = domain_mean_vector(final.base_mass),
                    mean_signed_delta_thickness = domain_mean_vector(last_delta_thickness_vec),
                    mean_abs_delta_thickness = domain_mean_abs_vector(last_delta_thickness_vec),
                    max_abs_delta_thickness = domain_max_abs_vector(last_delta_thickness_vec),
                    mean_signed_delta_wet_mass = domain_mean_vector(last_delta_wet_mass_vec),
                    mean_abs_delta_wet_mass = domain_mean_abs_vector(last_delta_wet_mass_vec),
                    max_abs_delta_wet_mass = domain_max_abs_vector(last_delta_wet_mass_vec),
                    mean_signed_delta_base_mass = domain_mean_vector(last_delta_base_mass_vec),
                    mean_abs_delta_base_mass = domain_mean_abs_vector(last_delta_base_mass_vec),
                    max_abs_delta_base_mass = domain_max_abs_vector(last_delta_base_mass_vec),
                )
            end
        end
        push!(history, record)

        time_block!(timings, :cycle_logging) do
            println(cycle_log_line(record))
        end

        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    final_state = if status == :max_cycles
        history[end].cycle == config.max_cycles ? final : prev
    else
        final
    end
    if config.write_outputs
        final_state_host, final_smb_ice_vec, final_runoff_vec, last_delta_thickness_host, last_delta_wet_mass_host,
        last_delta_base_mass_host, last_delta_ice_sheet_smb_host = time_block!(timings, :finalize_state_transfer) do
            final_state_host_local = config.backend == :gpu ? cpu_cycle_summary(final_state) : final_state
            final_smb_ice_vec_local = config.backend == :gpu ? Array(domain.smb_ice) : copy(domain.smb_ice)
            final_runoff_vec_local = config.backend == :gpu ? Array(domain.runoff) : copy(domain.runoff)
            last_delta_thickness_host_local = config.backend == :gpu ? Array(last_delta_thickness_vec) : last_delta_thickness_vec
            last_delta_wet_mass_host_local = config.backend == :gpu ? Array(last_delta_wet_mass_vec) : last_delta_wet_mass_vec
            last_delta_base_mass_host_local = config.backend == :gpu ? Array(last_delta_base_mass_vec) : last_delta_base_mass_vec
            last_delta_ice_sheet_smb_host_local = config.backend == :gpu ? Array(last_delta_ice_sheet_smb_vec) : last_delta_ice_sheet_smb_vec
            return (
                final_state_host_local,
                final_smb_ice_vec_local,
                final_runoff_vec_local,
                last_delta_thickness_host_local,
                last_delta_wet_mass_host_local,
                last_delta_base_mass_host_local,
                last_delta_ice_sheet_smb_host_local,
            )
        end
        final_thickness, final_wet_mass, final_bulk_density, final_base_mass, final_ice_sheet_smb, final_runoff,
        last_delta_thickness, last_delta_wet_mass, last_delta_base_mass, last_delta_ice_sheet_smb =
            time_block!(timings, :scatter_final_outputs) do
                return (
                    scatter_to_grid(final_state_host.thickness, js, is, (ny, nx)),
                    scatter_to_grid(final_state_host.wet_mass, js, is, (ny, nx)),
                    scatter_to_grid(final_state_host.bulk_density, js, is, (ny, nx)),
                    scatter_to_grid(final_state_host.base_mass, js, is, (ny, nx)),
                    scatter_to_grid(final_smb_ice_vec, js, is, (ny, nx)),
                    scatter_to_grid(final_runoff_vec, js, is, (ny, nx)),
                    scatter_to_grid(last_delta_thickness_host, js, is, (ny, nx)),
                    scatter_to_grid(last_delta_wet_mass_host, js, is, (ny, nx)),
                    scatter_to_grid(last_delta_base_mass_host, js, is, (ny, nx)),
                    scatter_to_grid(last_delta_ice_sheet_smb_host, js, is, (ny, nx)),
                )
            end
    end
    if config.write_netcdf
        final_domain = config.backend == :gpu ? time_block!(timings, :gpu_transfer) do
            SM.cpu_domain(domain)
        end : domain
        layer_grids = time_block!(timings, :collect_final_layer_grids) do
            collect_final_layer_grids(final_domain, js, is, (ny, nx), config.ntot)
        end
        monthly_sum_thickness_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_thickness)
        end : monthly_sum_thickness
        monthly_sum_wet_mass_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_wet_mass)
        end : monthly_sum_wet_mass
        monthly_sum_bulk_density_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_bulk_density)
        end : monthly_sum_bulk_density
        monthly_sum_base_mass_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_base_mass)
        end : monthly_sum_base_mass
        monthly_sum_ice_sheet_smb_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_ice_sheet_smb)
        end : monthly_sum_ice_sheet_smb
        monthly_sum_export_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_export)
        end : monthly_sum_export
        monthly_sum_net_ice_sheet_forcing_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_net_ice_sheet_forcing)
        end : monthly_sum_net_ice_sheet_forcing
        monthly_sum_runoff_host = gpu_netcdf_diagnostics ? time_block!(timings, :gpu_transfer) do
            Array(monthly_sum_runoff)
        end : monthly_sum_runoff
        monthly_mean_thickness = similar(monthly_sum_thickness_host)
        monthly_mean_wet_mass = similar(monthly_sum_wet_mass_host)
        monthly_mean_bulk_density = similar(monthly_sum_bulk_density_host)
        monthly_mean_base_mass = similar(monthly_sum_base_mass_host)
        monthly_mean_ice_sheet_smb = similar(monthly_sum_ice_sheet_smb_host)
        monthly_export_to_ice = similar(monthly_sum_export_host)
        monthly_net_ice_sheet_forcing = similar(monthly_sum_net_ice_sheet_forcing_host)
        monthly_runoff = similar(monthly_sum_runoff_host)
        monthly_mean_thickness_grid, monthly_mean_wet_mass_grid, monthly_mean_bulk_density_grid,
        monthly_mean_base_mass_grid, monthly_mean_ice_sheet_smb_grid, monthly_export_to_ice_grid,
        monthly_net_ice_sheet_forcing_grid, monthly_runoff_grid = time_block!(timings, :aggregate_monthly_outputs) do
            @inbounds for m in 1:nmonth_total
                c = max(monthly_count[m], 1)
                monthly_mean_thickness[m, :] .= monthly_sum_thickness_host[m, :] ./ c
                monthly_mean_wet_mass[m, :] .= monthly_sum_wet_mass_host[m, :] ./ c
                monthly_mean_bulk_density[m, :] .= monthly_sum_bulk_density_host[m, :] ./ c
                monthly_mean_base_mass[m, :] .= monthly_sum_base_mass_host[m, :] ./ c
                monthly_mean_ice_sheet_smb[m, :] .= monthly_sum_ice_sheet_smb_host[m, :]
                monthly_export_to_ice[m, :] .= monthly_sum_export_host[m, :]
                monthly_net_ice_sheet_forcing[m, :] .= monthly_sum_net_ice_sheet_forcing_host[m, :]
                monthly_runoff[m, :] .= monthly_sum_runoff_host[m, :]
            end
            return (
                monthly_vectors_to_grids(monthly_mean_thickness, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_mean_wet_mass, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_mean_bulk_density, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_mean_base_mass, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_mean_ice_sheet_smb, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_export_to_ice, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_net_ice_sheet_forcing, js, is, (ny, nx)),
                monthly_vectors_to_grids(monthly_runoff, js, is, (ny, nx)),
            )
        end
    end

    if config.write_outputs
        mkpath(config.out_dir)
        summary_path = joinpath(config.out_dir, "gris_equilibrium_summary.txt")
        history_csv_path = joinpath(config.out_dir, "gris_equilibrium_history.csv")
        history_plot_path = joinpath(config.out_dir, "gris_equilibrium_history.png")
        fields_plot_path = joinpath(config.out_dir, "gris_equilibrium_fields.png")
        time_block!(timings, :write_summary_text) do
            write_spinup_summary(summary_path, config, time_values, nvalid, length(mask), history, status, timings)
        end
        time_block!(timings, :write_history_csv) do
            write_spinup_history_csv(history_csv_path, history)
        end
        time_block!(timings, :render_history_plot) do
            render_spinup_history_plot(history_plot_path, history, config)
        end
        time_block!(timings, :render_fields_plot) do
            render_spinup_fields_plot(
                fields_plot_path,
                x,
                y,
                initial_thickness,
                final_thickness,
                final_bulk_density,
                final_base_mass,
                last_delta_thickness,
                last_delta_wet_mass,
                last_delta_base_mass,
                final_runoff,
                status,
                length(history),
            )
        end
    end
    if config.write_netcdf
        time_block!(timings, :write_netcdf) do
            finalize_spinup_netcdf!(
                writer,
                final_thickness,
                final_wet_mass,
                final_bulk_density,
                final_base_mass,
                final_ice_sheet_smb,
                last_delta_thickness,
                last_delta_wet_mass,
                last_delta_base_mass,
                last_delta_ice_sheet_smb,
                final_runoff,
                layer_grids,
                history,
                monthly_mean_thickness_grid,
                monthly_mean_wet_mass_grid,
                monthly_mean_bulk_density_grid,
                monthly_mean_base_mass_grid,
                monthly_mean_ice_sheet_smb_grid,
                monthly_export_to_ice_grid,
                monthly_net_ice_sheet_forcing_grid,
                monthly_runoff_grid,
                status,
                length(history),
                steps_written,
            )
        end
    end
    run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9

    println("GrIS equilibrium spin-up complete.")
    println("Forcing file   : $(abspath(config.nc_path))")
    println("Forcing start  : $(first(time_values))")
    println("Forcing end    : $(last(time_values))")
    println("Backend        : $(String(config.backend))")
    println("Cycles         : $(length(history))")
    println("Status         : $(string(status))")
    println(@sprintf("Simulation wall: %.3f s", simulation_wall_sec))
    println(@sprintf("Run wall total : %.3f s", run_wall_sec))
    if haskey(timings.totals, :model_step_wall)
        println(@sprintf("Model step wall: %.3f s", timings.totals[:model_step_wall]))
    end
    if config.write_netcdf
        println("Output NetCDF  : $(abspath(nc_path))")
    else
        println("Output NetCDF  : skipped (--no-nc)")
    end
    if config.write_outputs
        println("History CSV    : $(abspath(history_csv_path))")
        println("History plot   : $(abspath(history_plot_path))")
        println("Fields plot    : $(abspath(fields_plot_path))")
        println("Summary        : $(abspath(summary_path))")
    else
        println("File outputs   : skipped (--no-output)")
    end
    print_timing_summary(stdout, timings; total_wall_sec=run_wall_sec)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
