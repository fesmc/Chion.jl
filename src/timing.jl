"""
Timing helpers for profiling the major stages of `step!`.
"""

mutable struct StepTimingStats
    totals::Dict{Symbol, Float64}
    counts::Dict{Symbol, Int}
end

StepTimingStats() = StepTimingStats(Dict{Symbol, Float64}(), Dict{Symbol, Int}())

function add_timing!(stats::StepTimingStats, key::Symbol, dt_sec::Float64, count::Int=1)
    stats.totals[key] = get(stats.totals, key, 0.0) + dt_sec
    stats.counts[key] = get(stats.counts, key, 0) + count
    return dt_sec
end

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

@inline _time_block!(f::F, stats, key::Symbol) where {F<:Function} =
    _time_block!(stats, key, f)

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

function print_timing_summary(io::IO, stats::StepTimingStats)
    rows, total = timing_rows(stats)
    println(io, "Step timing summary")
    println(
        io,
        @sprintf("  %-24s %12s %9s %12s %10s", "stage", "total [s]", "share", "mean [micro s]", "count"),
    )
    for row in rows
        println(
            io,
            @sprintf(
                "  %-24s %12.3f %8.1f%% %12.3f %10d",
                String(row.key),
                row.total_sec,
                row.share_pct,
                row.mean_sec * 1.0e6,
                row.count,
            ),
        )
    end
    println(io, @sprintf("  %-24s %12.3f", "total_timed", total))
    return
end
