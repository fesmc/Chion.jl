#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Dates
using Printf
using Statistics
using Base.Threads
using Chion

include("run_gris_one_step.jl")

const DEFAULT_OUT_DIR_EQUIL = joinpath(@__DIR__, "..", "plots", "gris_equilibrium")
const LIBNETCDF = "/opt/homebrew/lib/libnetcdf.dylib"
const NC_NOERR = 0
const NC_CLOBBER = 0x0000
const NC_NETCDF4 = 0x1000
const NC_GLOBAL = -1
const NC_FLOAT = 5
const NC_DOUBLE = 6
const NC_INT = 4
const NC_UNLIMITED = 0

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

function timing_rows(stats::TimingStats)
    rows = NamedTuple[]
    total = sum(values(stats.totals))
    for key in keys(stats.totals)
        dt = stats.totals[key]
        count = stats.counts[key]
        push!(rows, (
            key = key,
            total_sec = dt,
            count = count,
            mean_sec = count > 0 ? dt / count : NaN,
            share_pct = total > 0.0 ? 100.0 * dt / total : 0.0,
        ))
    end
    sort!(rows; by=row -> row.total_sec, rev=true)
    return rows, total
end

function print_timing_summary(io::IO, stats::TimingStats)
    rows, total = timing_rows(stats)
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
    println(io, @sprintf("  %-24s %12.3f", "total_timed", total))
    if haskey(stats.totals, :model_step_wall)
        println(io)
        println(io, "  note: use model_step_wall for cross-run step-performance comparisons.")
        println(io, "  note: model_step sums per-thread elapsed time, so it is scheduler-sensitive.")
    end
    return
end

function print_spinup_help()
    println("Usage:")
    println("  julia --project=. examples/scripts/run_gris_equilibrium.jl [options]")
    println()
    println("Options:")
    println("  --nc=PATH                    MAR NetCDF/HDF5 file")
    println("  --out-dir=PATH               Output directory (default: examples/plots/gris_equilibrium)")
    println("  --out-nc=PATH                Output NetCDF path (default: OUT_DIR/gris_equilibrium_final_state.nc)")
    println("  --no-nc                      Skip NetCDF output and NetCDF-only step diagnostics for cleaner timing")
    println("  --no-output                  Skip all file output (implies --no-nc) for clean timing runs")
    println("  --mask-threshold=VALUE       Minimum MSK value for GrIS cells (default: 50)")
    println("  --ntot=N                     Chion maximum active layers (default: 80)")
    println("  --max-cycles=N               Maximum forcing-cycle repeats (default: 10)")
    println("  --tol-thickness=VALUE        Mean abs cycle delta-thickness tolerance in m (default: 1e-3)")
    println("  --tol-swe=VALUE              Mean abs cycle delta-SWE tolerance in mmWE (default: 0.1)")
    println("  --drift-window=N             Stop early when last N cycles show persistent drift (default: 3)")
    println("  --flip-turbulent-fluxes      Multiply SHF and LHF by -1 before forcing Chion")
    println("  --help                       Show this message")
end

function parse_spinup_config(args::Vector{String})
    nc_path = arg_value(args, "nc", DEFAULT_NC_PATH)
    isempty(nc_path) && error("Pass --nc=PATH or place the MAR file at $(DEFAULT_NC_PATH).")
    write_outputs = !has_flag(args, "no-output")
    return (
        nc_path = nc_path,
        out_dir = arg_value(args, "out-dir", DEFAULT_OUT_DIR_EQUIL),
        out_nc = arg_value(args, "out-nc", ""),
        write_outputs = write_outputs,
        write_netcdf = write_outputs && !has_flag(args, "no-nc"),
        mask_threshold = parse(Float64, arg_value(args, "mask-threshold", "50.0")),
        ntot = parse(Int, arg_value(args, "ntot", "20")),
        max_cycles = parse(Int, arg_value(args, "max-cycles", "10")),
        tol_thickness = parse(Float64, arg_value(args, "tol-thickness", "1.0e-3")),
        tol_swe = parse(Float64, arg_value(args, "tol-swe", "0.1")),
        drift_window = parse(Int, arg_value(args, "drift-window", "3")),
        turbulent_flux_sign = has_flag(args, "flip-turbulent-fluxes") ? -1.0 : 1.0,
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

function summarize_columns(columns::Vector{SM.SnowpackColumn})
    n = length(columns)
    thickness = Vector{Float64}(undef, n)
    wet_mass = Vector{Float64}(undef, n)
    bulk_density = Vector{Float64}(undef, n)
    base_mass = Vector{Float64}(undef, n)
    smb_ice = Vector{Float64}(undef, n)
    liquid_water = Vector{Float64}(undef, n)
    runoff = Vector{Float64}(undef, n)

    @threads :static for idx in eachindex(columns)
        summary = summarize_column(columns[idx])
        thickness[idx] = summary.thickness
        wet_mass[idx] = summary.wet_mass
        bulk_density[idx] = summary.bulk_density
        base_mass[idx] = columns[idx].mass_base
        smb_ice[idx] = columns[idx].smb_ice
        liquid_water[idx] = summary.liquid_mass
        runoff[idx] = summary.runoff
    end

    return (
        thickness = thickness,
        wet_mass = wet_mass,
        bulk_density = bulk_density,
        base_mass = base_mass,
        smb_ice = smb_ice,
        liquid_water = liquid_water,
        runoff = runoff,
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

function persistent_drift(history::Vector{NamedTuple}, window::Int, tol_thickness::Float64, tol_swe::Float64)
    length(history) >= window || return false
    recent = history[(end - window + 1):end]
    dth = [rec.mean_signed_delta_thickness for rec in recent]
    dswe = [rec.mean_signed_delta_wet_mass for rec in recent]

    same_sign(values) = all(v -> v > 0.0, values) || all(v -> v < 0.0, values)
    nearly_constant(values) = abs(values[end]) >= 0.8 * abs(values[1])

    return same_sign(dth) &&
           same_sign(dswe) &&
           minimum(abs.(dth)) > tol_thickness &&
           minimum(abs.(dswe)) > tol_swe &&
           nearly_constant(dth) &&
           nearly_constant(dswe)
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
        println(io, "MAR file           : ", abspath(config.nc_path))
        println(io, "Forcing start      : ", first(time_values))
        println(io, "Forcing end        : ", last(time_values))
        println(io, "Forcing steps      : ", length(time_values))
        println(io, "GrIS cells         : ", nvalid, " / ", ngrid)
        println(io, @sprintf("Mask threshold     : %.2f", config.mask_threshold))
        println(io, "Threads            : ", nthreads())
        println(io, "File output        : ", config.write_outputs ? "enabled" : "disabled (--no-output)")
        println(io, "NetCDF output      : ", config.write_netcdf ? "enabled" : "disabled (--no-nc)")
        println(io, "Status             : ", string(status))
        println(io, "Cycles completed   : ", length(history))
        println(io, @sprintf("Tol thickness (m)  : %.6g", config.tol_thickness))
        println(io, @sprintf("Tol SWE (mmWE)     : %.6g", config.tol_swe))
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
        if status == :drifting
            println(io)
            println(io, "Interpretation     : Forcing cycle shows persistent drift, so a snow equilibrium was not reached.")
        elseif status == :max_cycles
            println(io)
            println(io, "Interpretation     : Max cycles reached before convergence.")
        else
            println(io)
            println(io, "Interpretation     : Convergence thresholds were met.")
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
    P.hline!(p4, [config.tol_thickness]; color=:black, linestyle=:dash, label="tolerance")
    p5 = P.plot(cycles, mean_abs_dswe; lw=3, marker=:circle, color=:darkorange, xlabel="Cycle", ylabel="mmWE", title="Mean abs cycle dSWE", framestyle=:box)
    P.hline!(p5, [config.tol_swe]; color=:black, linestyle=:dash, label="tolerance")
    p6 = P.plot(cycles, mean_abs_dbase; lw=3, marker=:circle, color=:indigo, xlabel="Cycle", ylabel="mmWE", title="Mean abs cycle dBase", framestyle=:box)

    fig = P.plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(1700, 950), plot_title="Chion GrIS spin-up convergence history")
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
    columns::Vector{SM.SnowpackColumn},
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

    @inbounds for idx in eachindex(columns)
        col = columns[idx]
        j = js[idx]
        i = is[idx]
        n_active[j, i] = Int32(col.N)
        for k in 1:col.N
            rho = col.density[k]
            m = col.mass[k]
            mw = col.mass_w[k]
            layer_density[k, j, i] = rho
            layer_snow_mass[k, j, i] = m
            layer_liquid_mass[k, j, i] = mw
            layer_temperature_c[k, j, i] = col.temperature[k] - col.c.T0
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
    var_step_layer_temp::Cint
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
)
    mkpath(dirname(out_nc))
    isfile(out_nc) && rm(out_nc, force=true)

    ny, nx = size(mask)
    nlayer = config.ntot
    npoint = length(js)
    ntime = length(time_values)
    max_steps = config.max_cycles * ntime

    step_cycle = Vector{Int32}(undef, max_steps)
    step_source_index = Vector{Int32}(undef, max_steps)
    step_source_code = Vector{Int32}(undef, max_steps)
    step_counter = 0
    for cyc in 1:config.max_cycles
        for t in 1:ntime
            step_counter += 1
            step_cycle[step_counter] = Int32(cyc)
            step_source_index[step_counter] = Int32(t)
            ts = time_values[t]
            step_source_code[step_counter] = Int32(year(ts) * 1000000 + month(ts) * 10000 + day(ts) * 100 + hour(ts))
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
    dims_slyx = Cint[dim_step, dim_layer, dim_y, dim_x]

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
    nc_put_att_text(ncid, var_step, "long_name", "Sequential model step index across repeated annual cycles")
    var_step_cycle = nc_def_var(ncid, "step_cycle", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_cycle, "long_name", "Repeated annual forcing cycle index for each model step")
    var_step_source_index = nc_def_var(ncid, "step_source_index", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_index, "long_name", "1-based index into the original forcing year for each model step")
    var_step_source_code = nc_def_var(ncid, "step_source_code", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_code, "long_name", "Source forcing timestamp code YYYYMMDDHH for each model step")
    var_step_valid = nc_def_var(ncid, "step_valid", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_valid, "long_name", "1 where a model step was completed and written, 0 for unused trailing slots")

    var_mask = define_nc_output_variable(ncid, dims_yx, "gris_mask", "MAR ice-sheet mask", "1")
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
    var_monthly_smb_ice = define_nc_output_variable(ncid, dims_myx, "monthly_mean_ice_sheet_smb", "Monthly mean cumulative net mass forcing to the ice sheet", "mmWE")
    var_monthly_export = define_nc_output_variable(ncid, dims_myx, "monthly_export_to_ice", "Monthly firn mass exported to the ice model", "mmWE")
    var_monthly_net_ice = define_nc_output_variable(ncid, dims_myx, "monthly_net_ice_sheet_forcing", "Monthly net mass forcing to the ice sheet", "mmWE")
    var_monthly_runoff = define_nc_output_variable(ncid, dims_myx, "monthly_runoff", "Monthly runoff production", "mmWE")
    var_step_export = define_nc_output_variable(ncid, dims_syx, "step_export_to_ice", "Firn mass exported to the ice model for each model step", "mmWE")
    var_step_smb_ice = define_nc_output_variable(ncid, dims_syx, "step_ice_sheet_smb", "Net mass forcing to the ice sheet for each model step", "mmWE")
    var_step_layer_temp = define_nc_output_variable(
        ncid,
        dims_slyx,
        "step_layer_temperature_c",
        "Chion layer temperature for every model step",
        "C",
    )

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
    nc_put_att_text(ncid, NC_GLOBAL, "monthly_note", "Monthly fields are averages or sums over all Chion daily states within each repeated-forcing month.")
    nc_put_att_text(ncid, NC_GLOBAL, "step_note", "step_layer_temperature_c uses dimensions (step, layer, y, x); inactive layers and non-GrIS cells are stored as NaN.")
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
        var_step_layer_temp,
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

function write_step_layer_temperature!(writer::SpinupNetCDFWriter, step_index::Int, step_layer_temperature_c::Array{Float64, 3})
    nc_put_vara_float_4d_step_layer_yx(writer.ncid, writer.var_step_layer_temp, step_index, step_layer_temperature_c)
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
    shapes = time_block!(timings, :read_dataset_shapes) do
        read_dataset_shapes(config.nc_path)
    end
    date_codes, time_values = time_block!(timings, :read_mar_times) do
        read_mar_times(config.nc_path, shapes)
    end
    dt_days = [infer_dt_days(time_values, t) for t in eachindex(time_values)]

    x = time_block!(timings, :read_grid) do
        vec(read_hdf5_full(config.nc_path, "x", shapes))
    end
    y = time_block!(timings, :read_grid) do
        vec(read_hdf5_full(config.nc_path, "y", shapes))
    end
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
        read_full_timeseries_3d(config.nc_path, "SHF", shapes) .* config.turbulent_flux_sign
    end
    lhf_full = time_block!(timings, :read_forcing_lhf) do
        read_full_timeseries_3d(config.nc_path, "LHF", shapes) .* config.turbulent_flux_sign
    end
    u_wind_info = time_block!(timings, :read_forcing_wind) do
        read_first_available_timeseries_3d(config.nc_path, ["UU", "U10"], shapes)
    end
    v_wind_info = time_block!(timings, :read_forcing_wind) do
        read_first_available_timeseries_3d(config.nc_path, ["VV", "V10"], shapes)
    end
    wind_full = if !isnothing(u_wind_info) && !isnothing(v_wind_info)
        hypot.(u_wind_info.data, v_wind_info.data)
    else
        nothing
    end
    if isnothing(wind_full)
        println("Wind forcing: MAR wind components not found; using default 5.0 m s^-1.")
    else
        println(@sprintf(
            "Wind forcing: |V| from MAR components %s and %s.",
            u_wind_info.name,
            v_wind_info.name,
        ))
    end

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
    valid_mask = falses(ny, nx)
    @inbounds for j in 1:ny, i in 1:nx
        valid_mask[j, i] = isfinite(mask[j, i]) && mask[j, i] >= config.mask_threshold && isfinite(tt_full[1, j, i])
    end
    valid_indices = findall(valid_mask)
    nvalid = length(valid_indices)
    ntime = length(time_values)
    unique_month_keys = unique((year(t), month(t)) for t in time_values)
    month_lookup = Dict{Tuple{Int, Int}, Int}()
    for (idx, key) in enumerate(unique_month_keys)
        month_lookup[key] = idx
    end
    step_month = [month_lookup[(year(t), month(t))] for t in time_values]
    nmonth_per_cycle = length(unique_month_keys)
    nmonth_total = config.max_cycles * nmonth_per_cycle
    month_cycle = Int32[]
    month_of_year = Int32[]
    source_month_code = Int32[]
    for cyc in 1:config.max_cycles, key in unique_month_keys
        push!(month_cycle, Int32(cyc))
        push!(month_of_year, Int32(key[2]))
        push!(source_month_code, Int32(key[1] * 100 + key[2]))
    end

    js = Vector{Int}(undef, nvalid)
    is = Vector{Int}(undef, nvalid)
    columns = Vector{SM.SnowpackColumn}(undef, nvalid)
    initial_thickness_vec = Vector{Float64}(undef, nvalid)
    initial_wet_mass_vec = Vector{Float64}(undef, nvalid)

    tair_k = Matrix{Float64}(undef, nvalid, ntime)
    snow_rate = Matrix{Float64}(undef, nvalid, ntime)
    rain_rate = Matrix{Float64}(undef, nvalid, ntime)
    s_boa = Matrix{Float64}(undef, nvalid, ntime)
    q_lw = Matrix{Float64}(undef, nvalid, ntime)
    q_sh = Matrix{Float64}(undef, nvalid, ntime)
    q_lh = Matrix{Float64}(undef, nvalid, ntime)
    wind_speed = Matrix{Float64}(undef, nvalid, ntime)

    time_block!(timings, :initialize_columns_and_forcing) do
        @threads :static for idx in eachindex(valid_indices)
            j, i = Tuple(valid_indices[idx])
            js[idx] = j
            is[idx] = i

            column = build_column_from_mar(
                Float64(zn3_init[j, i]),
                @view(ro1_init[:, j, i]),
                @view(ti1_init[:, j, i]),
                @view(wa1_init[:, j, i]),
                outlay_bounds;
                ntot=config.ntot,
            )
            columns[idx] = column

            init_summary = summarize_column(column)
            initial_thickness_vec[idx] = init_summary.thickness
            initial_wet_mass_vec[idx] = init_summary.wet_mass

            for t in 1:ntime
                tair_k[idx, t] = valid_or(-15.0, Float64(tt_full[t, j, i])) + column.c.T0
                snow_rate[idx, t] = mmwe_day_to_kgm2s(Float64(sf_full[t, j, i]))
                rain_rate[idx, t] = mmwe_day_to_kgm2s(Float64(rf_full[t, j, i]))
                s_boa[idx, t] = valid_or(0.0, Float64(swd_full[t, j, i]))
                q_lw[idx, t] = Float64(lwd_full[t, j, i])
                q_sh[idx, t] = Float64(shf_full[t, j, i])
                q_lh[idx, t] = Float64(lhf_full[t, j, i])
                wind_speed[idx, t] = isnothing(wind_full) ? 5.0 : valid_or(5.0, Float64(wind_full[t, j, i]))
            end
        end
    end

    prev = time_block!(timings, :summarize_columns_initial) do
        summarize_columns(columns)
    end
    initial_thickness = config.write_outputs ? scatter_to_grid(initial_thickness_vec, js, is, (ny, nx)) : Matrix{Float64}(undef, 0, 0)
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
        )
    end : nothing
    step_export_to_ice = config.write_netcdf ? fill(NaN, ny, nx) : Matrix{Float64}(undef, 0, 0)
    step_ice_sheet_smb = config.write_netcdf ? fill(NaN, ny, nx) : Matrix{Float64}(undef, 0, 0)
    step_layer_temperature_c = config.write_netcdf ? fill(NaN, config.ntot, ny, nx) : Array{Float64, 3}(undef, 0, 0, 0)
    steps_written = 0
    history = NamedTuple[]
    status = :max_cycles
    last_delta_thickness_vec = fill(NaN, nvalid)
    last_delta_wet_mass_vec = fill(NaN, nvalid)
    last_delta_base_mass_vec = fill(NaN, nvalid)
    last_delta_ice_sheet_smb_vec = fill(NaN, nvalid)
    monthly_sum_thickness = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_wet_mass = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_bulk_density = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_base_mass = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_ice_sheet_smb = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_export = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_net_ice_sheet_forcing = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_sum_runoff = config.write_netcdf ? zeros(Float64, nmonth_total, nvalid) : Matrix{Float64}(undef, 0, 0)
    monthly_count = config.write_netcdf ? zeros(Int32, nmonth_total) : Int32[]
    final = prev
    simulation_wall_t0 = time_ns()

    for cycle in 1:config.max_cycles
        for t in 1:ntime
            dt = dt_days[t]
            month_idx = (cycle - 1) * nmonth_per_cycle + step_month[t]
            if config.write_netcdf
                fill!(step_export_to_ice, NaN)
                fill!(step_ice_sheet_smb, NaN)
                fill!(step_layer_temperature_c, NaN)
            end
            step_thread_sec = zeros(Float64, Threads.maxthreadid())
            diag_thread_sec = zeros(Float64, Threads.maxthreadid())
            step_wall_t0 = time_ns()
            @threads :static for idx in eachindex(columns)
                col = columns[idx]
                P_snow = snow_rate[idx, t]
                P_rain = rain_rate[idx, t]
                q_lw_down = isfinite(q_lw[idx, t]) ? q_lw[idx, t] : nothing
                q_sh_now = isfinite(q_sh[idx, t]) ? q_sh[idx, t] : nothing
                q_lh_now = isfinite(q_lh[idx, t]) ? q_lh[idx, t] : nothing
                wind_now = wind_speed[idx, t]
                base_before = col.mass_base
                smb_ice_before = col.smb_ice
                runoff_before = col.runoff
                tid = threadid()
                t0 = time_ns()
                SM.step!(
                    col,
                    tair_k[idx, t],
                    P_snow + P_rain,
                    dt;
                    p_snow=P_snow,
                    p_rain=P_rain,
                    s_boa=s_boa[idx, t],
                    wind_speed=wind_now,
                    q_lw_down=q_lw_down,
                    q_sh=q_sh_now,
                    q_lh=q_lh_now,
                )
                step_thread_sec[tid] += (time_ns() - t0) * 1.0e-9
                if config.write_netcdf
                    t1 = time_ns()
                    summary = summarize_column(col)
                    monthly_sum_thickness[month_idx, idx] += summary.thickness
                    monthly_sum_wet_mass[month_idx, idx] += summary.wet_mass
                    monthly_sum_bulk_density[month_idx, idx] += summary.bulk_density
                    monthly_sum_base_mass[month_idx, idx] += col.mass_base
                    monthly_sum_ice_sheet_smb[month_idx, idx] += col.smb_ice
                    monthly_sum_export[month_idx, idx] += col.mass_base - base_before
                    monthly_sum_net_ice_sheet_forcing[month_idx, idx] += col.smb_ice - smb_ice_before
                    monthly_sum_runoff[month_idx, idx] += col.runoff - runoff_before
                    step_export_to_ice[js[idx], is[idx]] = col.mass_base - base_before
                    step_ice_sheet_smb[js[idx], is[idx]] = col.smb_ice - smb_ice_before
                    if col.N > 0
                        @views step_layer_temperature_c[1:col.N, js[idx], is[idx]] .= col.temperature[1:col.N] .- col.c.T0
                    end
                    diag_thread_sec[tid] += (time_ns() - t1) * 1.0e-9
                end
            end
            step_wall_sec = (time_ns() - step_wall_t0) * 1.0e-9
            add_timing!(timings, :model_step, sum(step_thread_sec), length(columns))
            add_timing!(timings, :model_step_wall, step_wall_sec, length(columns))
            if config.write_netcdf
                add_timing!(timings, :step_diagnostics, sum(diag_thread_sec), length(columns))
                monthly_count[month_idx] += 1
                steps_written += 1
                time_block!(timings, :step_output_write) do
                    write_step_export_to_ice!(writer, steps_written, step_export_to_ice)
                    write_step_ice_sheet_smb!(writer, steps_written, step_ice_sheet_smb)
                    write_step_layer_temperature!(writer, steps_written, step_layer_temperature_c)
                end
            end
        end

        final = time_block!(timings, :summarize_columns_cycle) do
            summarize_columns(columns)
        end
        last_delta_thickness_vec .= final.thickness .- prev.thickness
        last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
        last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
        last_delta_ice_sheet_smb_vec .= final.smb_ice .- prev.smb_ice

        record = (
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
        push!(history, record)

        println(
            @sprintf(
                "cycle=%d mean_th=%.5f m mean_swe=%.5f mmWE mean_base=%.5f mmWE mean_abs_dth=%.5f m mean_abs_dswe=%.5f mmWE mean_abs_dbase=%.5f mmWE",
                record.cycle,
                record.mean_thickness,
                record.mean_wet_mass,
                record.mean_base_mass,
                record.mean_abs_delta_thickness,
                record.mean_abs_delta_wet_mass,
                record.mean_abs_delta_base_mass,
            ),
        )

        if record.mean_abs_delta_thickness <= config.tol_thickness &&
           record.mean_abs_delta_wet_mass <= config.tol_swe
            status = :converged
            break
        elseif persistent_drift(history, config.drift_window, config.tol_thickness, config.tol_swe)
            status = :drifting
            break
        end

        prev = final
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9

    if config.write_outputs
        final_thickness = scatter_to_grid(final.thickness, js, is, (ny, nx))
        final_wet_mass = scatter_to_grid(final.wet_mass, js, is, (ny, nx))
        final_bulk_density = scatter_to_grid(final.bulk_density, js, is, (ny, nx))
        final_base_mass = scatter_to_grid(final.base_mass, js, is, (ny, nx))
        final_ice_sheet_smb = scatter_to_grid(final.smb_ice, js, is, (ny, nx))
        final_runoff = scatter_to_grid(final.runoff, js, is, (ny, nx))
        last_delta_thickness = scatter_to_grid(last_delta_thickness_vec, js, is, (ny, nx))
        last_delta_wet_mass = scatter_to_grid(last_delta_wet_mass_vec, js, is, (ny, nx))
        last_delta_base_mass = scatter_to_grid(last_delta_base_mass_vec, js, is, (ny, nx))
        last_delta_ice_sheet_smb = scatter_to_grid(last_delta_ice_sheet_smb_vec, js, is, (ny, nx))
    end
    if config.write_netcdf
        layer_grids = time_block!(timings, :collect_final_layer_grids) do
            collect_final_layer_grids(columns, js, is, (ny, nx), config.ntot)
        end
        monthly_mean_thickness = similar(monthly_sum_thickness)
        monthly_mean_wet_mass = similar(monthly_sum_wet_mass)
        monthly_mean_bulk_density = similar(monthly_sum_bulk_density)
        monthly_mean_base_mass = similar(monthly_sum_base_mass)
        monthly_mean_ice_sheet_smb = similar(monthly_sum_ice_sheet_smb)
        monthly_export_to_ice = similar(monthly_sum_export)
        monthly_net_ice_sheet_forcing = similar(monthly_sum_net_ice_sheet_forcing)
        monthly_runoff = similar(monthly_sum_runoff)
        monthly_mean_thickness_grid, monthly_mean_wet_mass_grid, monthly_mean_bulk_density_grid,
        monthly_mean_base_mass_grid, monthly_mean_ice_sheet_smb_grid, monthly_export_to_ice_grid,
        monthly_net_ice_sheet_forcing_grid, monthly_runoff_grid = time_block!(timings, :aggregate_monthly_outputs) do
            @inbounds for m in 1:nmonth_total
                c = max(monthly_count[m], 1)
                monthly_mean_thickness[m, :] .= monthly_sum_thickness[m, :] ./ c
                monthly_mean_wet_mass[m, :] .= monthly_sum_wet_mass[m, :] ./ c
                monthly_mean_bulk_density[m, :] .= monthly_sum_bulk_density[m, :] ./ c
                monthly_mean_base_mass[m, :] .= monthly_sum_base_mass[m, :] ./ c
                monthly_mean_ice_sheet_smb[m, :] .= monthly_sum_ice_sheet_smb[m, :] ./ c
                monthly_export_to_ice[m, :] .= monthly_sum_export[m, :]
                monthly_net_ice_sheet_forcing[m, :] .= monthly_sum_net_ice_sheet_forcing[m, :]
                monthly_runoff[m, :] .= monthly_sum_runoff[m, :]
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
    println("MAR file       : $(abspath(config.nc_path))")
    println("Forcing start  : $(first(time_values))")
    println("Forcing end    : $(last(time_values))")
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
    print_timing_summary(stdout, timings)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
