"""
Timing helpers for profiling the major stages of `step!`.
"""

"""
    StepTimingStats

Accumulator for profiling the major internal stages of [`step!`](@ref).
"""
mutable struct StepTimingStats
    totals::Dict{Symbol, Float64}
    counts::Dict{Symbol, Int}
end

"""
    StepTimingStats()

Create an empty accumulator for `step!` stage timings. The returned object is
mutated in-place by `add_timing!`, `_time_block!`, and `_time_call!`.
"""
StepTimingStats() = StepTimingStats(Dict{Symbol, Float64}(), Dict{Symbol, Int}())

"""
    add_timing!(stats, key, dt_sec, count=1)

Accumulate `dt_sec` seconds under `key` inside `stats` and increment the call
count by `count`. Mutates `stats` and returns `dt_sec`.
"""
function add_timing!(stats::StepTimingStats, key::Symbol, dt_sec::Float64, count::Int=1)
    stats.totals[key] = get(stats.totals, key, 0.0) + dt_sec
    stats.counts[key] = get(stats.counts, key, 0) + count
    return dt_sec
end

"""
    _time_block!(stats, key, f)

Execute the zero-argument function `f` and optionally record its wall-clock
runtime under `key`. When `stats === nothing`, the call is forwarded without
measurement.
"""
@inline function _time_block!(
    stats::Nothing,
    key::Symbol,
    f::F,
) where {F<:Function}
    return f()
end

@inline function _time_block!(
    stats::StepTimingStats,
    key::Symbol,
    f::F,
) where {F<:Function}
    t0 = time_ns()
    value = f()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

"""
    _time_block!(f, stats, key)

Compatibility argument order for `_time_block!`. Executes `f` and records
timing information when `stats` is a `StepTimingStats`.
"""
@inline _time_block!(f::F, stats, key::Symbol) where {F<:Function} =
    _time_block!(stats, key, f)

"""
    _time_call!(stats, key, f, args...)

Call `f(args...)` and optionally accumulate its runtime under `key`. The
function returns the value produced by `f` and never allocates timing state on
its own.
"""
@inline function _time_call!(
    stats::Nothing,
    key::Symbol,
    f,
    args...,
)
    return f(args...)
end

@inline function _time_call!(
    stats::StepTimingStats,
    key::Symbol,
    f,
    args...,
)
    t0 = time_ns()
    value = f(args...)
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

"""
    timing_rows(stats)

Return `(rows, total)` summary data for a `StepTimingStats`. Each row contains
the total time, invocation count, mean time, and percentage share for one
timed stage.
"""
function timing_rows(stats::StepTimingStats)
    rows = NamedTuple[]
    total = sum(values(stats.totals))
    for key in keys(stats.totals)
        dt = stats.totals[key]
        count = stats.counts[key]
        push!(rows, (
            key=key,
            total_sec=dt,
            count=count,
            mean_sec=count > 0 ? dt / count : NaN,
            share_pct=total > 0.0 ? 100.0 * dt / total : 0.0,
        ))
    end
    sort!(rows; by=row -> row.total_sec, rev=true)
    return rows, total
end

"""
    print_timing_summary(io, stats; total_wall_sec=nothing)

Write a human-readable table of accumulated `step!` stage timings to `io`.
This is purely diagnostic and does not modify model state.
"""
function print_timing_summary(
    io::IO,
    stats::StepTimingStats;
    total_wall_sec::Union{Nothing, Float64}=nothing,
)
    rows, total = timing_rows(stats)
    share_total = isnothing(total_wall_sec) ? total : total_wall_sec
    println(io, "Timing summary")
    println(io, @sprintf("  %-24s %12s %9s %12s %10s", "stage", "total [s]", "share", "mean [ms]", "count"))
    for row in rows
        share_pct = share_total > 0.0 ? 100.0 * row.total_sec / share_total : 0.0
        println(io, @sprintf("  %-24s %12.3f %8.1f%% %12.3f %10d", String(row.key), row.total_sec, share_pct, row.mean_sec * 1.0e3, row.count))
    end
    if isnothing(total_wall_sec)
        println(io, @sprintf("  %-24s %12.3f", "total_accounted", total))
        return
    end
    unaccounted = max(total_wall_sec - total, 0.0)
    println(io, @sprintf("  %-24s %12.3f %8.1f%% %12s %10s", "unaccounted", unaccounted, share_total > 0.0 ? 100.0 * unaccounted / share_total : 0.0, "", ""))
    println(io, @sprintf("  %-24s %12.3f", "total_accounted", total))
    println(io, @sprintf("  %-24s %12.3f", "run_wall_total", total_wall_sec))
end
