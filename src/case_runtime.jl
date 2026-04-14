const SM = SnowpackModel
const CASE_OUTPUT_GROUPS = (
    final=(
        :final_thickness,
        :final_wet_mass,
        :final_bulk_density,
        :final_base_mass,
        :final_ice_sheet_smb,
        :final_runoff,
        :last_cycle_delta_thickness,
        :last_cycle_delta_wet_mass,
        :last_cycle_delta_base_mass,
        :last_cycle_delta_ice_sheet_smb,
    ),
    layers=(
        :n_active,
        :layer_density,
        :layer_thickness,
        :layer_snow_mass,
        :layer_liquid_mass,
        :layer_temperature_c,
    ),
    history=(
        :history_mean_thickness,
        :history_mean_wet_mass,
        :history_mean_bulk_density,
        :history_mean_base_mass,
        :history_mean_abs_delta_thickness,
        :history_mean_abs_delta_wet_mass,
        :history_mean_abs_delta_base_mass,
    ),
    monthly=(
        :monthly_mean_thickness,
        :monthly_mean_wet_mass,
        :monthly_mean_bulk_density,
        :monthly_mean_base_mass,
        :monthly_mean_ice_sheet_smb,
        :monthly_export_to_ice,
        :monthly_net_ice_sheet_forcing,
        :monthly_runoff,
    ),
    step=(:step_export_to_ice, :step_ice_sheet_smb),
)
const CASE_NETCDF_VARIABLE_GROUPS = Dict(key => collect(values) for (key, values) in pairs(CASE_OUTPUT_GROUPS))
const CASE_NETCDF_VARIABLES = unique(Symbol[var for group in values(CASE_OUTPUT_GROUPS) for var in group])
const FINAL_GRID_KEYS = CASE_OUTPUT_GROUPS.final
const LAYER_GRID_KEYS = CASE_OUTPUT_GROUPS.layers
const MONTHLY_GRID_KEYS = CASE_OUTPUT_GROUPS.monthly
const HISTORY_OUTPUT_SPECS = (
    (output=:history_mean_thickness, record=:mean_thickness),
    (output=:history_mean_wet_mass, record=:mean_wet_mass),
    (output=:history_mean_bulk_density, record=:mean_bulk_density),
    (output=:history_mean_base_mass, record=:mean_base_mass),
    (output=:history_mean_abs_delta_thickness, record=:mean_abs_delta_thickness),
    (output=:history_mean_abs_delta_wet_mass, record=:mean_abs_delta_wet_mass),
    (output=:history_mean_abs_delta_base_mass, record=:mean_abs_delta_base_mass),
)
const TimingStats = SM.StepTimingStats
const add_timing! = SM.add_timing!
const timing_rows = SM.timing_rows

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

function print_timing_summary(io::IO, stats::StepTimingStats; total_wall_sec::Union{Nothing, Float64}=nothing)
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

using Dates
using Base.Threads: @threads, nthreads
using NCDatasets
import CUDA
import Libdl

include("cases/runtime_core.jl")
include("cases/netcdf.jl")
include("cases/runtime_execute.jl")
