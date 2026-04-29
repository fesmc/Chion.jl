"""
Timing helpers for profiling the major stages of simulation runs.
Backed by TimerOutputs.jl.
"""

using TimerOutputs: TimerOutput, TimerOutputs

"""
    StepTimingStats

Accumulator for profiling major internal stages. Wraps a
`TimerOutputs.TimerOutput` and a per-key item count separate from call count.
"""
mutable struct StepTimingStats
    to::TimerOutput
    item_counts::Dict{Symbol, Int}
end

"""
    StepTimingStats()

Create an empty timing accumulator.
"""
StepTimingStats() = StepTimingStats(TimerOutput(), Dict{Symbol, Int}())

"""
    add_timing!(stats, key, dt_sec, count=1)

Record `dt_sec` seconds under `key` and increment the item count by `count`.
"""
function add_timing!(stats::StepTimingStats, key::Symbol, dt_sec::Float64, count::Int=1)
    name = String(key)
    dt_ns = round(Int64, dt_sec * 1e9)
    child = get!(stats.to.inner_timers, name) do
        to = TimerOutput()
        to.name = name
        to
    end
    data = child.accumulated_data
    child.accumulated_data = TimerOutputs.TimeData(
        data.ncalls + 1,
        data.time + dt_ns,
        data.allocs,
        data.firstexec == 0 ? time_ns() : data.firstexec,
    )
    stats.item_counts[key] = get(stats.item_counts, key, 0) + count
    return dt_sec
end

"""
    timing_rows(stats)

Return `(rows, total)`, where `rows` is a vector of per-stage `NamedTuple`s and
`total` is the sum of all recorded times in seconds.
"""
function timing_rows(stats::StepTimingStats)
    rows = NamedTuple[]
    total = sum(TimerOutputs.time(t) for t in values(stats.to.inner_timers); init=Int64(0)) * 1e-9
    for (name, timer) in stats.to.inner_timers
        dt = TimerOutputs.time(timer) * 1e-9
        ncalls = TimerOutputs.ncalls(timer)
        key = Symbol(name)
        item_count = get(stats.item_counts, key, ncalls)
        share_pct = total > 0.0 ? 100.0 * dt / total : 0.0
        push!(rows, (
            key=key,
            total_sec=dt,
            count=item_count,
            mean_sec=ncalls > 0 ? dt / ncalls : NaN,
            share_pct=share_pct,
        ))
    end
    sort!(rows; by=row -> row.total_sec, rev=true)
    return rows, total
end

@inline _time_block!(::Nothing, ::Symbol, f::F) where {F <: Function} = f()

@inline function _time_block!(stats::StepTimingStats, key::Symbol, f::F) where {F <: Function}
    t0 = time_ns()
    value = f()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

@inline _time_block!(f::F, stats, key::Symbol) where {F <: Function} = _time_block!(stats, key, f)

@inline _time_call!(::Nothing, ::Symbol, f, args...) = f(args...)

@inline function _time_call!(stats::StepTimingStats, key::Symbol, f, args...)
    t0 = time_ns()
    value = f(args...)
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

"""
    time_block!(stats, key, f; synchronize=nothing)

Execute `f()` and record elapsed wall time under `key`.
"""
function time_block!(stats, key::Symbol, f; synchronize=nothing)
    synchronize === nothing || synchronize()
    t0 = time_ns()
    value = f()
    synchronize === nothing || synchronize()
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

time_block!(f, stats, key::Symbol; kwargs...) = time_block!(stats, key, f; kwargs...)

"""
    time_counted_block!(stats, key, count, f; synchronize=nothing)

Execute `f()` and record elapsed wall time under `key`, adding `count` to the
item count for that key.
"""
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

"""
    print_timing_summary(io, stats; total_wall_sec=nothing)

Write a human-readable table of accumulated stage timings to `io`.
"""
function print_timing_summary(
    io::IO,
    stats::StepTimingStats;
    total_wall_sec::Union{Nothing, Float64}=nothing,
)
    show(io, stats.to; sortby=:time)
    println(io)
    if !isnothing(total_wall_sec)
        _, total = timing_rows(stats)
        unaccounted = max(total_wall_sec - total, 0.0)
        println(io, @sprintf("  %-24s %12.3f", "unaccounted", unaccounted))
        println(io, @sprintf("  %-24s %12.3f", "run_wall_total", total_wall_sec))
    end
end
