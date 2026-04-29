# Year history, run reports, and timing wrappers for Simulation runs.

using Base.Threads: nthreads
using TimerOutputs: TimerOutputs
using Statistics: mean
using CSV

@inline completed_year_count(history::Vector{NamedTuple}, status::Symbol, years::Int) =
    status === :complete ? years : isempty(history) ? 0 : min(history[end].year, years)

@inline year_metrics_schedule_label(stride::Int) =
    stride == 0 ? "final year only" : stride == 1 ? "every year" : "every $(stride) years + final"

const HISTORY_CSV_SPECS = (
    (key=:year, label="year", integer=true),
    (key=:mean_thickness, label="mean_thickness_m", integer=false),
    (key=:mean_wet_mass, label="mean_wet_mass_mmwe", integer=false),
    (key=:mean_bulk_density, label="mean_bulk_density_kgm3", integer=false),
    (key=:mean_base_mass, label="mean_base_mass_mmwe", integer=false),
    (key=:mean_signed_delta_thickness, label="mean_signed_delta_thickness_m", integer=false),
    (key=:mean_abs_delta_thickness, label="mean_abs_delta_thickness_m", integer=false),
    (key=:max_abs_delta_thickness, label="max_abs_delta_thickness_m", integer=false),
    (key=:mean_signed_delta_wet_mass, label="mean_signed_delta_wet_mass_mmwe", integer=false),
    (key=:mean_abs_delta_wet_mass, label="mean_abs_delta_wet_mass_mmwe", integer=false),
    (key=:max_abs_delta_wet_mass, label="max_abs_delta_wet_mass_mmwe", integer=false),
    (key=:mean_signed_delta_base_mass, label="mean_signed_delta_base_mass_mmwe", integer=false),
    (key=:mean_abs_delta_base_mass, label="mean_abs_delta_base_mass_mmwe", integer=false),
    (key=:max_abs_delta_base_mass, label="max_abs_delta_base_mass_mmwe", integer=false),
)

const SUMMARY_REPORT_SPECS = (
    (title="Final domain means", fields=(
        ("Thickness (m)", :mean_thickness),
        ("Wet mass (mmWE)", :mean_wet_mass),
        ("Bulk density (kg m-3)", :mean_bulk_density),
        ("Firn-to-ice mass (mmWE)", :mean_base_mass),
    )),
    (title="Last year deltas", fields=(
        ("Mean signed dThickness (m)", :mean_signed_delta_thickness),
        ("Mean abs dThickness (m)", :mean_abs_delta_thickness),
        ("Max abs dThickness (m)", :max_abs_delta_thickness),
        ("Mean signed dSWE (mmWE)", :mean_signed_delta_wet_mass),
        ("Mean abs dSWE (mmWE)", :mean_abs_delta_wet_mass),
        ("Max abs dSWE (mmWE)", :max_abs_delta_wet_mass),
        ("Mean signed dBase (mmWE)", :mean_signed_delta_base_mass),
        ("Mean abs dBase (mmWE)", :mean_abs_delta_base_mass),
        ("Max abs dBase (mmWE)", :max_abs_delta_base_mass),
    )),
)

@inline function _finite_mean(data)
    finite = filter(isfinite, data)
    return isempty(finite) ? NaN : mean(finite)
end

function _delta_stats(data)
    finite = filter(isfinite, data)
    isempty(finite) && return (mean_signed=NaN, mean_abs=NaN, max_abs=NaN)
    abs_vals = abs.(finite)
    return (
        mean_signed=mean(finite),
        mean_abs=mean(abs_vals),
        max_abs=maximum(abs_vals),
    )
end

function make_year_record_and_deltas!(
    year::Int,
    delta_thickness,
    delta_wet_mass,
    delta_base_mass,
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    prev_thickness,
    prev_wet_mass,
    prev_base_mass,
)
    delta_thickness .= thickness .- prev_thickness
    delta_wet_mass .= wet_mass .- prev_wet_mass
    delta_base_mass .= base_mass .- prev_base_mass
    dth = _delta_stats(_host_vector(delta_thickness))
    dwet = _delta_stats(_host_vector(delta_wet_mass))
    dbase = _delta_stats(_host_vector(delta_base_mass))
    return (
        year=year,
        mean_thickness=_finite_mean(_host_vector(thickness)),
        mean_wet_mass=_finite_mean(_host_vector(wet_mass)),
        mean_bulk_density=_finite_mean(_host_vector(bulk_density)),
        mean_base_mass=_finite_mean(_host_vector(base_mass)),
        mean_signed_delta_thickness=dth.mean_signed,
        mean_abs_delta_thickness=dth.mean_abs,
        max_abs_delta_thickness=dth.max_abs,
        mean_signed_delta_wet_mass=dwet.mean_signed,
        mean_abs_delta_wet_mass=dwet.mean_abs,
        max_abs_delta_wet_mass=dwet.max_abs,
        mean_signed_delta_base_mass=dbase.mean_signed,
        mean_abs_delta_base_mass=dbase.mean_abs,
        max_abs_delta_base_mass=dbase.max_abs,
    )
end

function _copy_summary_fields!(summary, device_summary, fields)
    for field in fields
        copyto!(getfield(summary, field), getfield(device_summary, field))
    end
    return summary
end

function year_log_line(record)
    return @sprintf(
        "year=%d mean_th=%.5f m mean_swe=%.5f mmWE mean_base=%.5f mmWE mean_abs_dth=%.5f m mean_abs_dswe=%.5f mmWE mean_abs_dbase=%.5f mmWE",
        record.year,
        record.mean_thickness,
        record.mean_wet_mass,
        record.mean_base_mass,
        record.mean_abs_delta_thickness,
        record.mean_abs_delta_wet_mass,
        record.mean_abs_delta_base_mass,
    )
end

function write_run_history_csv(out_path::AbstractString, history::Vector{NamedTuple})
    mkpath(dirname(out_path))
    # Build a NamedTuple with the published column labels from the history records.
    cols = NamedTuple{Tuple(Symbol(spec.label) for spec in HISTORY_CSV_SPECS)}(
        Tuple(getfield.(history, spec.key) for spec in HISTORY_CSV_SPECS)
    )
    CSV.write(out_path, cols)
end

function write_run_summary(
    out_path::AbstractString,
    options::RunOptions,
    time_values::Vector{DateTime},
    ncol::Int,
    history::Vector{NamedTuple},
    status::Symbol,
    timings::StepTimingStats,
)
    last_record = isempty(history) ? nothing : history[end]
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, options.name)
        println(io, "Input label        : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
        println(io, "Forcing start      : ", first(time_values))
        println(io, "Forcing end        : ", last(time_values))
        println(io, "Forcing steps      : ", length(time_values))
        println(io, "Columns            : ", ncol)
        println(io, "Backend            : ", String(options.backend))
        println(io, "Threads            : ", nthreads())
        println(io, "File output        : ", options.write_outputs ? "enabled" : "disabled (--no-output)")
        println(io, "NetCDF output      : ", options.write_netcdf ? "enabled" : "disabled (--no-nc)")
        println(io, "Year metrics       : ", year_metrics_schedule_label(options.history_year_stride))
        println(io, "Status             : ", string(status))
        println(io, "Years completed    : ", completed_year_count(history, status, options.years))
        if last_record !== nothing
            for section in SUMMARY_REPORT_SPECS
                println(io)
                println(io, section.title)
                for (label, key) in section.fields
                    println(io, @sprintf("%-28s : %.6f", label, getfield(last_record, key)))
                end
            end
        end
        println(io)
        println(io, "Interpretation     : ", status === :complete ? "Requested years completed." : "Run finalized before all requested years completed.")
        println(io)
        print_timing_summary(io, timings)
    end
end

function print_run_report(
    io::IO,
    options::RunOptions,
    time_values::Vector{DateTime},
    history::Vector{NamedTuple},
    status::Symbol,
    simulation_wall_sec::Float64,
    run_wall_sec::Float64,
    timings::StepTimingStats;
    nc_path::AbstractString="",
    summary_path::AbstractString="",
    history_csv_path::AbstractString="",
)
    println(io, "$(options.name) complete.")
    println(io, "Input label     : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
    println(io, "Forcing start   : $(first(time_values))")
    println(io, "Forcing end     : $(last(time_values))")
    println(io, "Backend         : ", String(options.backend))
    println(io, "Years           : ", completed_year_count(history, status, options.years))
    println(io, "Status          : ", string(status))
    println(io, "Year metrics    : ", year_metrics_schedule_label(options.history_year_stride))
    println(io, @sprintf("Simulation wall : %.3f s", simulation_wall_sec))
    println(io, @sprintf("Run wall total  : %.3f s", run_wall_sec))
    haskey(timings.to.inner_timers, "model_step_wall") && println(io, @sprintf("Model step wall : %.3f s", TimerOutputs.time(timings.to.inner_timers["model_step_wall"]) * 1e-9))
    println(io, "Output NetCDF   : ", options.write_netcdf ? abspath(nc_path) : "skipped (--no-nc)")
    if options.write_outputs
        println(io, "History CSV     : $(abspath(history_csv_path))")
        println(io, "Summary         : $(abspath(summary_path))")
    else
        println(io, "File outputs    : skipped (--no-output)")
    end
    print_timing_summary(io, timings; total_wall_sec=run_wall_sec)
end
function time_block!(stats, key::Symbol, f; synchronize=nothing)
    synchronize === nothing || synchronize()
    t0 = time_ns()
    value = f()
    synchronize === nothing || synchronize()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

time_block!(f, stats, key::Symbol; kwargs...) = time_block!(stats, key, f; kwargs...)

function time_counted_block!(stats, key::Symbol, count::Int, f; synchronize=nothing)
    synchronize === nothing || synchronize()
    t0 = time_ns()
    value = f()
    synchronize === nothing || synchronize()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9, count)
    return value
end

time_counted_block!(f, stats, key::Symbol, count::Int; kwargs...) =
    time_counted_block!(stats, key, count, f; kwargs...)
