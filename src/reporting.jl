# Year history, run reports, and timing wrappers for Simulation runs.

using TimerOutputs: TimerOutputs

@inline completed_year_count(history::YearHistory, status::Symbol, years::Int) =
    status === :complete ? years : isempty(history) ? 0 : min(history[end].year, years)

@inline year_metrics_schedule_label(stride::Int, enabled::Bool=true) =
    !enabled ? "disabled" : stride == 0 ? "final year only" : stride == 1 ? "every year" : "every $(stride) years + final"

@inline function _finite_mean(data)
    total = 0.0
    count = 0
    for value in data
        if isfinite(value)
            total += value
            count += 1
        end
    end
    return count == 0 ? NaN : total / count
end

function _delta_stats(data)
    signed_total = 0.0
    abs_total = 0.0
    max_abs = 0.0
    count = 0
    for value in data
        if isfinite(value)
            abs_value = abs(value)
            signed_total += value
            abs_total += abs_value
            max_abs = max(max_abs, abs_value)
            count += 1
        end
    end
    count == 0 && return (mean_signed=NaN, mean_abs=NaN, max_abs=NaN)
    return (
        mean_signed=signed_total / count,
        mean_abs=abs_total / count,
        max_abs=max_abs,
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

function print_run_report(
    io::IO,
    options::RunOptions,
    time_values::Vector{DateTime},
    history::YearHistory,
    status::Symbol,
    run_wall_sec::Float64,
    timings::StepTimingStats;
    nc_path::AbstractString="",
)
    println(io, "$(options.name) complete.")
    println(io, "Input label     : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
    println(io, "Forcing start   : $(first(time_values))")
    println(io, "Forcing end     : $(last(time_values))")
    println(io, "Backend         : ", String(options.backend))
    println(io, "Years           : ", completed_year_count(history, status, options.years))
    println(io, "Status          : ", string(status))
    println(io, "Year metrics    : ", year_metrics_schedule_label(options.history_year_stride, options.compute_year_metrics))
    println(io, @sprintf("Run wall total  : %.3f s", run_wall_sec))
    haskey(timings.to.inner_timers, "model_step_wall") && println(io, @sprintf("Model step wall : %.3f s", TimerOutputs.time(timings.to.inner_timers["model_step_wall"]) * 1e-9))
    println(io, "Output NetCDF   : ", isempty(nc_path) ? "skipped" : abspath(nc_path))
    print_timing_summary(io, timings; total_wall_sec=run_wall_sec)
end
