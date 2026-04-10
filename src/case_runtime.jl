const SM = SnowpackModel

using Dates
using Base.Threads: @threads, nthreads
import CUDA
import Libdl

const NC_NOERR = 0
const NC_CLOBBER = 0x0000
const NC_NETCDF4 = 0x1000
const NC_GLOBAL = -1
const NC_FLOAT = 5
const NC_DOUBLE = 6
const NC_INT = 4

"""
    CASE_NETCDF_VARIABLE_GROUPS

Named groups of NetCDF output variables accepted by [`RunConfig`](@ref) and
the case-running helpers.
"""
const CASE_NETCDF_VARIABLE_GROUPS = Dict(
    :final => [
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
    ],
    :layers => [
        :n_active,
        :layer_density,
        :layer_thickness,
        :layer_snow_mass,
        :layer_liquid_mass,
        :layer_temperature_c,
    ],
    :history => [
        :history_mean_thickness,
        :history_mean_wet_mass,
        :history_mean_bulk_density,
        :history_mean_base_mass,
        :history_mean_abs_delta_thickness,
        :history_mean_abs_delta_wet_mass,
        :history_mean_abs_delta_base_mass,
    ],
    :monthly => [
        :monthly_mean_thickness,
        :monthly_mean_wet_mass,
        :monthly_mean_bulk_density,
        :monthly_mean_base_mass,
        :monthly_mean_ice_sheet_smb,
        :monthly_export_to_ice,
        :monthly_net_ice_sheet_forcing,
        :monthly_runoff,
    ],
    :step => [
        :step_export_to_ice,
        :step_ice_sheet_smb,
    ],
)

"""
    CASE_NETCDF_VARIABLES

Flat list of all supported NetCDF output variable names.
"""
const CASE_NETCDF_VARIABLES = unique(vcat(values(CASE_NETCDF_VARIABLE_GROUPS)...))
const _LIBNETCDF_CACHE = Ref{Union{Nothing, String}}(nothing)

"""
    TimingStats

Run-level timing accumulator used by the high-level case runtime.
"""
mutable struct TimingStats
    totals::Dict{Symbol, Float64}
    counts::Dict{Symbol, Int}
end

TimingStats() = TimingStats(Dict{Symbol, Float64}(), Dict{Symbol, Int}())

"""
    ForcingData

Normalized forcing bundle used by case execution. All meteorological fields
share the same `(ncol, ntime)` shape and `dt_days` stores one duration per time
step.
"""
struct ForcingData
    time_values::Vector{DateTime}
    dt_days::Vector{Float64}
    air_temperature
    snowfall_rate
    rainfall_rate
    shortwave_down
    wind_speed
    q_lw_down
    has_q_lw_down
    q_sh
    has_q_sh
    q_lh
    has_q_lh
end

"""
    GridLayout

Mapping between column indices and a regular output grid for gridded case
inputs and NetCDF export.
"""
struct GridLayout
    x::Vector{Float64}
    y::Vector{Float64}
    js::Vector{Int}
    is::Vector{Int}
    mask::Matrix{Float64}
end

"""
    SnowpackStateFields

Plain Julia container for initializing a [`SnowpackDomain`](@ref) from explicit
state arrays and domain parameters.
"""
struct SnowpackStateFields
    N::Vector{Int}
    mass::Matrix{Float64}
    mass_w::Matrix{Float64}
    density::Matrix{Float64}
    temperature::Matrix{Float64}
    mass_base::Vector{Float64}
    smb_ice::Vector{Float64}
    runoff::Vector{Float64}
    Tsrf::Vector{Float64}
    snow_cover::Vector{Float64}
    albedo_dynamic::Vector{Float64}
    physics::SM.SnowpackPhysicalConstants{Float64}
    mass_max::Float64
    mass_split::Float64
    mass_min::Float64
    rho_max::Float64
end

"""
    RunConfig

High-level execution settings for [`build_case`](@ref) and [`run_case`](@ref),
including backend selection and output options.
"""
struct RunConfig
    name::String
    input_label::String
    output_dir::String
    netcdf_path::String
    write_outputs::Bool
    write_netcdf::Bool
    netcdf_variables::Vector{Symbol}
    cycles::Int
    backend::Symbol
    history_stride::Int
end

struct CaseNetCDFWriter
    ncid::Cint
    vars::Dict{Symbol, Cint}
    max_steps::Int
    cycles::Int
end

"""
    RunResult

Collected outputs and diagnostics returned by [`run_case`](@ref).
"""
struct RunResult
    history::Vector{NamedTuple}
    status::Symbol
    timings::TimingStats
    simulation_wall_sec::Float64
    run_wall_sec::Float64
    netcdf_path::String
    summary_path::String
    history_csv_path::String
    domain
    run::RunConfig
end

function add_timing!(stats::TimingStats, key::Symbol, dt_sec::Float64, count::Int=1)
    stats.totals[key] = get(stats.totals, key, 0.0) + dt_sec
    stats.counts[key] = get(stats.counts, key, 0) + count
    return dt_sec
end

@inline _maybe_synchronize!(::Nothing) = nothing
@inline _maybe_synchronize!(sync::F) where {F <: Function} = (sync(); nothing)

function time_block!(stats::TimingStats, key::Symbol, f::F; synchronize=nothing) where {F <: Function}
    _maybe_synchronize!(synchronize)
    t0 = time_ns()
    value = f()
    _maybe_synchronize!(synchronize)
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9)
    return value
end

time_block!(f::F, stats::TimingStats, key::Symbol; kwargs...) where {F <: Function} = time_block!(stats, key, f; kwargs...)

function time_counted_block!(stats::TimingStats, key::Symbol, count::Int, f::F; synchronize=nothing) where {F <: Function}
    _maybe_synchronize!(synchronize)
    t0 = time_ns()
    value = f()
    _maybe_synchronize!(synchronize)
    add_timing!(stats, key, (time_ns() - t0) * 1.0e-9, count)
    return value
end

time_counted_block!(f::F, stats::TimingStats, key::Symbol, count::Int; kwargs...) where {F <: Function} =
    time_counted_block!(stats, key, count, f; kwargs...)

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

@inline function normalize_case_backend(backend)
    value = lowercase(strip(String(backend)))
    value == "cpu" && return :threads
    value in ("threads", "gpu") || error("Unsupported backend '$backend'. Use `threads`, `cpu`, or `gpu`.")
    return Symbol(value)
end

@inline function normalize_history_stride(stride::Integer)
    value = Int(stride)
    value >= 0 || error("`history_stride` must be >= 0.")
    return value
end

@inline function should_record_cycle_metrics(cycle::Int, cycles::Int, stride::Int)
    cycle == cycles && return true
    stride == 0 && return false
    return mod(cycle, stride) == 0
end

@inline function completed_cycle_count(history::Vector{NamedTuple}, status::Symbol, cycles::Int)
    isempty(history) && return 0
    return status == :cycles ? cycles : history[end].cycle
end

@inline function cycle_metrics_schedule_label(stride::Int)
    stride == 0 && return "final cycle only"
    stride == 1 && return "every cycle"
    return "every $(stride) cycles + final"
end

@inline function _looks_like_directory_path(path::AbstractString)
    isempty(path) && return false
    return endswith(path, '/') || endswith(path, '\\')
end

@inline function resolve_case_netcdf_path(options::RunConfig)
    default_name = "$(options.name)_final_state.nc"
    isempty(options.netcdf_path) && return joinpath(options.output_dir, default_name)
    return (isdir(options.netcdf_path) || _looks_like_directory_path(options.netcdf_path)) ?
        joinpath(options.netcdf_path, default_name) :
        options.netcdf_path
end

function normalize_case_netcdf_variables(spec)
    if spec isa AbstractVector
        tokens = String[string(x) for x in spec]
    else
        text = lowercase(strip(String(spec)))
        isempty(text) && return copy(CASE_NETCDF_VARIABLES)
        tokens = split(text, ',')
    end

    selected = Symbol[]
    allowed_groups = sort!(String.(collect(keys(CASE_NETCDF_VARIABLE_GROUPS))))
    for token in tokens
        stripped = strip(token)
        isempty(stripped) && continue
        key = Symbol(lowercase(stripped))
        if key == :all
            append!(selected, CASE_NETCDF_VARIABLES)
        elseif key == :none
            continue
        elseif haskey(CASE_NETCDF_VARIABLE_GROUPS, key)
            append!(selected, CASE_NETCDF_VARIABLE_GROUPS[key])
        elseif key in CASE_NETCDF_VARIABLES
            push!(selected, key)
        else
            error(
                "Unsupported NetCDF variable selector '$token'. " *
                "Use `all`, `none`, a group ($(join(allowed_groups, ", "))), or an explicit variable name.",
            )
        end
    end
    return unique(selected)
end

function RunConfig(;
    name::AbstractString="snowpack_case",
    input_label::AbstractString="",
    output_dir::AbstractString="",
    netcdf_path::AbstractString="",
    write_outputs::Bool=true,
    write_netcdf::Bool=true,
    netcdf_variables=copy(CASE_NETCDF_VARIABLES),
    cycles::Integer=10,
    backend=:threads,
    history_stride::Integer=1,
)
    return RunConfig(
        String(name),
        String(input_label),
        String(output_dir),
        String(netcdf_path),
        write_outputs,
        write_netcdf,
        normalize_case_netcdf_variables(netcdf_variables),
        Int(cycles),
        normalize_case_backend(backend),
        normalize_history_stride(history_stride),
    )
end

function _synthesized_time_values(dt_days::Vector{Float64})
    base = DateTime(2000, 1, 1, 12)
    out = Vector{DateTime}(undef, length(dt_days))
    elapsed_ms = 0
    for idx in eachindex(dt_days)
        out[idx] = base + Dates.Millisecond(elapsed_ms)
        elapsed_ms += round(Int, dt_days[idx] * 86_400_000)
    end
    return out
end

function _ensure_matching_field_sizes(reference::Tuple{Int, Int}, name::AbstractString, field)
    size(field) == reference || error("`$name` must have shape $(reference), got $(size(field)).")
    return
end

function _forcing_column_count(field, ntime::Int)
    field isa Number && return 1
    data = collect(field)
    ndims(data) == 1 && return 1
    ndims(data) == 2 || error("Forcing fields must be scalars, vectors, or matrices.")
    size(data, 2) == ntime ||
        error("Matrix forcing fields must have $ntime columns, got $(size(data, 2)).")
    return size(data, 1)
end

function _forcing_numeric_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Number
        return fill(Float64(field), ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return repeat(reshape(Float64.(data), 1, ntime), ncol, 1)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Float64}(data)
    end
    error("`$name` must be a scalar, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

function _forcing_bool_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Bool
        return fill(field, ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return repeat(reshape(Bool.(data), 1, ntime), ncol, 1)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Bool.(data)
    end
    error("`$name` must be a Bool, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

function ForcingData(;
    dt_days,
    air_temperature=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    air_temperature_c=nothing,
    snowfall_mm_day=nothing,
    rainfall_mm_day=nothing,
    shortwave_down,
    ncol::Union{Nothing, Integer}=nothing,
    wind_speed=nothing,
    q_lw_down=nothing,
    has_q_lw_down=nothing,
    q_sh=nothing,
    has_q_sh=nothing,
    q_lh=nothing,
    has_q_lh=nothing,
    time_values=nothing,
)
    has_native_inputs = !isnothing(air_temperature) || !isnothing(snowfall_rate) || !isnothing(rainfall_rate)
    has_user_inputs = !isnothing(air_temperature_c) || !isnothing(snowfall_mm_day) || !isnothing(rainfall_mm_day)
    has_native_inputs && has_user_inputs &&
        error("Pass either model-native forcing fields or user-facing forcing fields, not both.")

    dt_days_v = Float64.(collect(dt_days))
    isempty(dt_days_v) && error("`dt_days` must not be empty.")
    ntime = length(dt_days_v)
    all(>(0.0), dt_days_v) || error("All `dt_days` entries must be positive.")

    column_count = 1
    if isnothing(ncol)
        for candidate in (
            air_temperature,
            snowfall_rate,
            rainfall_rate,
            air_temperature_c,
            snowfall_mm_day,
            rainfall_mm_day,
            shortwave_down,
            wind_speed,
            q_lw_down,
            q_sh,
            q_lh,
        )
            if !isnothing(candidate)
                column_count = _forcing_column_count(candidate, ntime)
                break
            end
        end
    else
        column_count = Int(ncol)
    end
    column_count > 0 || error("`ncol` must be positive.")

    if has_user_inputs
        isnothing(air_temperature_c) && error("`air_temperature_c` is required.")
        isnothing(snowfall_mm_day) && error("`snowfall_mm_day` is required.")
        isnothing(rainfall_mm_day) && error("`rainfall_mm_day` is required.")
        air_temperature = _forcing_numeric_matrix(air_temperature_c, column_count, ntime, "air_temperature_c") .+ 273.15
        snowfall_rate = _forcing_numeric_matrix(snowfall_mm_day, column_count, ntime, "snowfall_mm_day") ./ 86_400.0
        rainfall_rate = _forcing_numeric_matrix(rainfall_mm_day, column_count, ntime, "rainfall_mm_day") ./ 86_400.0
    else
        isnothing(air_temperature) && error("`air_temperature` or `air_temperature_c` is required.")
        isnothing(snowfall_rate) && error("`snowfall_rate` or `snowfall_mm_day` is required.")
        isnothing(rainfall_rate) && error("`rainfall_rate` or `rainfall_mm_day` is required.")
        air_temperature = _forcing_numeric_matrix(air_temperature, column_count, ntime, "air_temperature")
        snowfall_rate = _forcing_numeric_matrix(snowfall_rate, column_count, ntime, "snowfall_rate")
        rainfall_rate = _forcing_numeric_matrix(rainfall_rate, column_count, ntime, "rainfall_rate")
    end

    shortwave_down_m = _forcing_numeric_matrix(shortwave_down, column_count, ntime, "shortwave_down")
    dims = size(air_temperature)
    wind_speed_m = isnothing(wind_speed) ? fill(5.0, dims) : _forcing_numeric_matrix(wind_speed, column_count, ntime, "wind_speed")
    q_lw_down_m = isnothing(q_lw_down) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_lw_down, column_count, ntime, "q_lw_down")
    has_q_lw_down_m = isnothing(q_lw_down) ?
        fill(false, dims) :
        (isnothing(has_q_lw_down) ? fill(true, dims) : _forcing_bool_matrix(has_q_lw_down, column_count, ntime, "has_q_lw_down"))
    q_sh_m = isnothing(q_sh) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_sh, column_count, ntime, "q_sh")
    has_q_sh_m = isnothing(q_sh) ?
        fill(false, dims) :
        (isnothing(has_q_sh) ? fill(true, dims) : _forcing_bool_matrix(has_q_sh, column_count, ntime, "has_q_sh"))
    q_lh_m = isnothing(q_lh) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_lh, column_count, ntime, "q_lh")
    has_q_lh_m = isnothing(q_lh) ?
        fill(false, dims) :
        (isnothing(has_q_lh) ? fill(true, dims) : _forcing_bool_matrix(has_q_lh, column_count, ntime, "has_q_lh"))

    for (name, field) in (
        ("snowfall_rate", snowfall_rate),
        ("rainfall_rate", rainfall_rate),
        ("shortwave_down", shortwave_down_m),
        ("wind_speed", wind_speed_m),
        ("q_lw_down", q_lw_down_m),
        ("has_q_lw_down", has_q_lw_down_m),
        ("q_sh", q_sh_m),
        ("has_q_sh", has_q_sh_m),
        ("q_lh", q_lh_m),
        ("has_q_lh", has_q_lh_m),
    )
        _ensure_matching_field_sizes(dims, name, field)
    end
    time_values_v = if isnothing(time_values)
        _synthesized_time_values(dt_days_v)
    else
        DateTime.(collect(time_values))
    end
    length(time_values_v) == dims[2] || error("`time_values` must have one entry per forcing timestep.")
    return ForcingData(
        time_values_v,
        dt_days_v,
        air_temperature,
        snowfall_rate,
        rainfall_rate,
        shortwave_down_m,
        wind_speed_m,
        q_lw_down_m,
        has_q_lw_down_m,
        q_sh_m,
        has_q_sh_m,
        q_lh_m,
        has_q_lh_m,
    )
end

function GridLayout(
    x::AbstractVector,
    y::AbstractVector,
    js::AbstractVector{<:Integer},
    is::AbstractVector{<:Integer},
    mask,
)
    x_v = Float64.(x)
    y_v = Float64.(y)
    js_v = Int.(js)
    is_v = Int.(is)
    mask_m = Matrix{Float64}(mask)
    length(js_v) == length(is_v) || error("`js` and `is` must have the same length.")
    size(mask_m, 1) == length(y_v) || error("`mask` y dimension must match `y`.")
    size(mask_m, 2) == length(x_v) || error("`mask` x dimension must match `x`.")
    return GridLayout(x_v, y_v, js_v, is_v, mask_m)
end

function SnowpackStateFields(
    N::AbstractVector{<:Integer},
    mass,
    mass_w,
    density,
    temperature;
    mass_base=zeros(Float64, length(N)),
    smb_ice=zeros(Float64, length(N)),
    runoff=zeros(Float64, length(N)),
    Tsrf=fill(SM.SnowpackPhysicalConstants().T0, length(N)),
    snow_cover=zeros(Float64, length(N)),
    albedo_dynamic=fill(SM.SnowpackPhysicalConstants().alpha_dry, length(N)),
    physics::SM.SnowpackPhysicalConstants{Float64}=SM.SnowpackPhysicalConstants(),
    mass_max::Real=SM.DEFAULT_MASS_MAX,
    mass_split::Real=SM.DEFAULT_MASS_SPLIT,
    mass_min::Real=SM.DEFAULT_MASS_MIN,
    rho_max::Real=SM.DEFAULT_RHO_MAX,
)
    N_v = Int.(N)
    mass_m = Matrix{Float64}(mass)
    mass_w_m = Matrix{Float64}(mass_w)
    density_m = Matrix{Float64}(density)
    temperature_m = Matrix{Float64}(temperature)
    ncol = length(N_v)
    size(mass_m, 2) == ncol || error("`mass` must have one column per entry of `N`.")
    size(mass_w_m) == size(mass_m) || error("`mass_w` must match `mass`.")
    size(density_m) == size(mass_m) || error("`density` must match `mass`.")
    size(temperature_m) == size(mass_m) || error("`temperature` must match `mass`.")
    return SnowpackStateFields(
        N_v,
        mass_m,
        mass_w_m,
        density_m,
        temperature_m,
        Float64.(mass_base),
        Float64.(smb_ice),
        Float64.(runoff),
        Float64.(Tsrf),
        Float64.(snow_cover),
        Float64.(albedo_dynamic),
        physics,
        Float64(mass_max),
        Float64(mass_split),
        Float64(mass_min),
        Float64(rho_max),
    )
end

function SM.SnowpackDomain(state::SnowpackStateFields)
    return SM.SnowpackDomain(
        state.N,
        state.mass,
        state.mass_w,
        state.density,
        state.temperature,
        state.mass_base,
        state.smb_ice,
        state.runoff,
        state.Tsrf,
        state.snow_cover,
        state.albedo_dynamic;
        c=state.physics,
        mass_max=state.mass_max,
        mass_split=state.mass_split,
        mass_min=state.mass_min,
        rho_max=state.rho_max,
    )
end

@inline case_output_enabled(options::RunConfig) = options.write_outputs || options.write_netcdf

@inline function case_selected(options::RunConfig, group::Symbol)
    group_vars = Set(CASE_NETCDF_VARIABLE_GROUPS[group])
    return any(var -> var in group_vars, options.netcdf_variables)
end

@inline function _grid_shape(layout::GridLayout)
    return size(layout.mask)
end

function _netcdf_search_dirs()
    dirs = String[]
    for key in ("NETCDF_DIR", "NETCDF_HOME", "HOMEBREW_PREFIX")
        root = strip(get(ENV, key, ""))
        isempty(root) || push!(dirs, joinpath(expanduser(root), "lib"))
    end
    for key in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH")
        value = strip(get(ENV, key, ""))
        isempty(value) && continue
        append!(dirs, filter(!isempty, split(value, ':')))
    end
    append!(dirs, Base.DL_LOAD_PATH)
    append!(dirs, ("/opt/homebrew/lib", "/usr/local/lib", "/opt/local/lib", "/usr/lib", "/lib"))
    return unique(filter(isdir, map(expanduser, dirs)))
end

function _netcdf_library_candidates()
    candidates = String[]
    env_lib = strip(get(ENV, "NETCDF_LIB", ""))
    if !isempty(env_lib)
        push!(candidates, expanduser(env_lib))
    end
    append!(candidates, ("libnetcdf", "libnetcdf.so", "libnetcdf.dylib"))
    for dir in _netcdf_search_dirs()
        append!(candidates, (
            joinpath(dir, "libnetcdf"),
            joinpath(dir, "libnetcdf.so"),
            joinpath(dir, "libnetcdf.dylib"),
        ))
    end
    return unique(candidates)
end

function resolve_libnetcdf()
    for candidate in _netcdf_library_candidates()
        try
            handle = Libdl.dlopen(candidate)
            Libdl.dlclose(handle)
            return candidate
        catch
        end
    end
    searched_dirs = join(_netcdf_search_dirs(), ", ")
    suggestion = Sys.isapple() && isfile("/opt/homebrew/lib/libnetcdf.dylib") ?
        " Try `export NETCDF_LIB=/opt/homebrew/lib/libnetcdf.dylib`." : ""
    error(
        "Could not load NetCDF library. Set NETCDF_LIB to the shared library path or load a NetCDF module." *
        (isempty(searched_dirs) ? "" : " Searched: " * searched_dirs * ".") *
        suggestion,
    )
end

function _libnetcdf()
    if isnothing(_LIBNETCDF_CACHE[])
        _LIBNETCDF_CACHE[] = resolve_libnetcdf()
    end
    return _LIBNETCDF_CACHE[]::String
end

@inline function nc_check(code::Integer)
    if code != NC_NOERR
        msg = unsafe_string(ccall((:nc_strerror, _libnetcdf()), Cstring, (Cint,), code))
        error("NetCDF error: $msg")
    end
    return
end

function nc_create(path::AbstractString)
    isdir(path) && error("NetCDF output path '$(abspath(path))' is a directory; pass a file path ending in `.nc`.")
    ncid = Ref{Cint}()
    code = ccall((:nc_create, _libnetcdf()), Cint, (Cstring, Cint, Ref{Cint}), path, NC_CLOBBER | NC_NETCDF4, ncid)
    if code != NC_NOERR
        msg = unsafe_string(ccall((:nc_strerror, _libnetcdf()), Cstring, (Cint,), code))
        error("NetCDF error while creating '$(abspath(path))': $msg")
    end
    return ncid[]
end

function nc_close(ncid::Cint)
    nc_check(ccall((:nc_close, _libnetcdf()), Cint, (Cint,), ncid))
    return
end

function nc_enddef(ncid::Cint)
    nc_check(ccall((:nc_enddef, _libnetcdf()), Cint, (Cint,), ncid))
    return
end

function nc_redef(ncid::Cint)
    nc_check(ccall((:nc_redef, _libnetcdf()), Cint, (Cint,), ncid))
    return
end

function nc_def_dim(ncid::Cint, name::AbstractString, len::Integer)
    dimid = Ref{Cint}()
    nc_check(ccall((:nc_def_dim, _libnetcdf()), Cint, (Cint, Cstring, Csize_t, Ref{Cint}), ncid, name, len, dimid))
    return dimid[]
end

function nc_def_var(ncid::Cint, name::AbstractString, xtype::Integer, dimids::Vector{Cint})
    varid = Ref{Cint}()
    nc_check(ccall((:nc_def_var, _libnetcdf()), Cint, (Cint, Cstring, Cint, Cint, Ptr{Cint}, Ref{Cint}), ncid, name, xtype, length(dimids), dimids, varid))
    return varid[]
end

function nc_put_att_text(ncid::Cint, varid::Integer, name::AbstractString, value::AbstractString)
    nc_check(ccall((:nc_put_att_text, _libnetcdf()), Cint, (Cint, Cint, Cstring, Csize_t, Cstring), ncid, Cint(varid), name, sizeof(value), value))
    return
end

function nc_put_att_float(ncid::Cint, varid::Integer, name::AbstractString, value::Float32)
    buf = Ref{Float32}(value)
    nc_check(ccall((:nc_put_att_float, _libnetcdf()), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ref{Float32}), ncid, Cint(varid), name, NC_FLOAT, 1, buf))
    return
end

function nc_put_var_double(ncid::Cint, varid::Integer, data::Vector{Float64})
    nc_check(ccall((:nc_put_var_double, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cdouble}), ncid, Cint(varid), data))
    return
end

function nc_put_var_int(ncid::Cint, varid::Integer, data::Vector{Int32})
    nc_check(ccall((:nc_put_var_int, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cint}), ncid, Cint(varid), data))
    return
end

function nc_put_var_int_2d(ncid::Cint, varid::Integer, data::AbstractMatrix{Int32})
    buf = permutedims(data, (2, 1))
    nc_check(ccall((:nc_put_var_int, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cint}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_2d(ncid::Cint, varid::Integer, data::AbstractMatrix{<:Real})
    buf = permutedims(Float32.(data), (2, 1))
    nc_check(ccall((:nc_put_var_float, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_3d(ncid::Cint, varid::Integer, data::Array{Float64, 3})
    buf = permutedims(Float32.(data), (3, 2, 1))
    nc_check(ccall((:nc_put_var_float, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
    return
end

function nc_put_var_float_1d(ncid::Cint, varid::Integer, data::Vector{Float64})
    buf = Float32.(data)
    nc_check(ccall((:nc_put_var_float, _libnetcdf()), Cint, (Cint, Cint, Ptr{Cfloat}), ncid, Cint(varid), buf))
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
            (:nc_put_vara_float, _libnetcdf()),
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

function monthly_vectors_to_grids(values::Matrix{Float64}, js::Vector{Int}, is::Vector{Int}, grid_shape::Tuple{Int, Int})
    nmonth, nvalid = size(values)
    ny, nx = grid_shape
    out = fill(NaN, nmonth, ny, nx)
    @inbounds for m in 1:nmonth, idx in 1:nvalid
        out[m, js[idx], is[idx]] = values[m, idx]
    end
    return out
end

function allocate_summary_buffers(n::Int)
    return (
        thickness = Vector{Float64}(undef, n),
        wet_mass = Vector{Float64}(undef, n),
        bulk_density = Vector{Float64}(undef, n),
        base_mass = Vector{Float64}(undef, n),
        smb_ice = Vector{Float64}(undef, n),
        liquid_water = Vector{Float64}(undef, n),
        runoff = Vector{Float64}(undef, n),
    )
end

function allocate_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int)
    return (
        thickness = similar(domain.mass, Float64, n),
        wet_mass = similar(domain.mass, Float64, n),
        bulk_density = similar(domain.mass, Float64, n),
        base_mass = similar(domain.mass, Float64, n),
        smb_ice = similar(domain.mass, Float64, n),
        liquid_water = similar(domain.mass, Float64, n),
        runoff = similar(domain.mass, Float64, n),
    )
end

function allocate_cycle_summary_buffers(n::Int)
    return (
        thickness = Vector{Float64}(undef, n),
        wet_mass = Vector{Float64}(undef, n),
        bulk_density = Vector{Float64}(undef, n),
        base_mass = Vector{Float64}(undef, n),
    )
end

function allocate_cycle_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int)
    return (
        thickness = similar(domain.mass, Float64, n),
        wet_mass = similar(domain.mass, Float64, n),
        bulk_density = similar(domain.mass, Float64, n),
        base_mass = similar(domain.mass, Float64, n),
    )
end

const CYCLE_METRIC_BUFFER_LENGTH = 13

struct CycleMetricsWorkspace{A <: AbstractVector{Float64}}
    device_buffer::A
    host_buffer::Vector{Float64}
end

function CycleMetricsWorkspace(domain::SM.AbstractSnowpackDomain)
    return CycleMetricsWorkspace(
        similar(domain.mass, Float64, CYCLE_METRIC_BUFFER_LENGTH),
        zeros(Float64, CYCLE_METRIC_BUFFER_LENGTH),
    )
end

function _cycle_metrics_kernel!(
    metrics,
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
    idx = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
    if idx <= length(thickness)
        thickness_val = thickness[idx]
        wet_mass_val = wet_mass[idx]
        bulk_density_val = bulk_density[idx]
        base_mass_val = base_mass[idx]
        delta_thickness_val = thickness_val - prev_thickness[idx]
        delta_wet_mass_val = wet_mass_val - prev_wet_mass[idx]
        delta_base_mass_val = base_mass_val - prev_base_mass[idx]
        delta_thickness[idx] = delta_thickness_val
        delta_wet_mass[idx] = delta_wet_mass_val
        delta_base_mass[idx] = delta_base_mass_val
        abs_delta_thickness_val = abs(delta_thickness_val)
        abs_delta_wet_mass_val = abs(delta_wet_mass_val)
        abs_delta_base_mass_val = abs(delta_base_mass_val)
        CUDA.@atomic metrics[1] += thickness_val
        CUDA.@atomic metrics[2] += wet_mass_val
        CUDA.@atomic metrics[3] += bulk_density_val
        CUDA.@atomic metrics[4] += base_mass_val
        CUDA.@atomic metrics[5] += delta_thickness_val
        CUDA.@atomic metrics[6] += abs_delta_thickness_val
        CUDA.@atomic metrics[7] = max(metrics[7], abs_delta_thickness_val)
        CUDA.@atomic metrics[8] += delta_wet_mass_val
        CUDA.@atomic metrics[9] += abs_delta_wet_mass_val
        CUDA.@atomic metrics[10] = max(metrics[10], abs_delta_wet_mass_val)
        CUDA.@atomic metrics[11] += delta_base_mass_val
        CUDA.@atomic metrics[12] += abs_delta_base_mass_val
        CUDA.@atomic metrics[13] = max(metrics[13], abs_delta_base_mass_val)
    end
    return nothing
end

@inline function _metric_mean(total::Float64, count::Int)
    return count == 0 ? NaN : total / count
end

@inline function _metric_max(value::Float64, has_value::Bool)
    return has_value ? value : NaN
end

function cycle_record_from_metrics(cycle::Int, metrics::AbstractVector{Float64}, n::Int)
    scale = inv(Float64(n))
    return (
        cycle = cycle,
        mean_thickness = metrics[1] * scale,
        mean_wet_mass = metrics[2] * scale,
        mean_bulk_density = metrics[3] * scale,
        mean_base_mass = metrics[4] * scale,
        mean_signed_delta_thickness = metrics[5] * scale,
        mean_abs_delta_thickness = metrics[6] * scale,
        max_abs_delta_thickness = metrics[7],
        mean_signed_delta_wet_mass = metrics[8] * scale,
        mean_abs_delta_wet_mass = metrics[9] * scale,
        max_abs_delta_wet_mass = metrics[10],
        mean_signed_delta_base_mass = metrics[11] * scale,
        mean_abs_delta_base_mass = metrics[12] * scale,
        max_abs_delta_base_mass = metrics[13],
    )
end

function make_cycle_record_and_deltas!(
    cycle::Int,
    delta_thickness::Vector{Float64},
    delta_wet_mass::Vector{Float64},
    delta_base_mass::Vector{Float64},
    thickness::Vector{Float64},
    wet_mass::Vector{Float64},
    bulk_density::Vector{Float64},
    base_mass::Vector{Float64},
    prev_thickness::Vector{Float64},
    prev_wet_mass::Vector{Float64},
    prev_base_mass::Vector{Float64},
)
    sum_thickness = 0.0
    sum_wet_mass = 0.0
    sum_bulk_density = 0.0
    sum_base_mass = 0.0
    count_thickness = 0
    count_wet_mass = 0
    count_bulk_density = 0
    count_base_mass = 0
    sum_delta_thickness = 0.0
    sum_abs_delta_thickness = 0.0
    max_abs_delta_thickness = 0.0
    count_delta_thickness = 0
    has_delta_thickness = false
    sum_delta_wet_mass = 0.0
    sum_abs_delta_wet_mass = 0.0
    max_abs_delta_wet_mass = 0.0
    count_delta_wet_mass = 0
    has_delta_wet_mass = false
    sum_delta_base_mass = 0.0
    sum_abs_delta_base_mass = 0.0
    max_abs_delta_base_mass = 0.0
    count_delta_base_mass = 0
    has_delta_base_mass = false
    @inbounds for idx in eachindex(thickness)
        thickness_val = thickness[idx]
        wet_mass_val = wet_mass[idx]
        bulk_density_val = bulk_density[idx]
        base_mass_val = base_mass[idx]
        delta_thickness_val = thickness_val - prev_thickness[idx]
        delta_wet_mass_val = wet_mass_val - prev_wet_mass[idx]
        delta_base_mass_val = base_mass_val - prev_base_mass[idx]
        delta_thickness[idx] = delta_thickness_val
        delta_wet_mass[idx] = delta_wet_mass_val
        delta_base_mass[idx] = delta_base_mass_val
        if isfinite(thickness_val)
            sum_thickness += thickness_val
            count_thickness += 1
        end
        if isfinite(wet_mass_val)
            sum_wet_mass += wet_mass_val
            count_wet_mass += 1
        end
        if isfinite(bulk_density_val)
            sum_bulk_density += bulk_density_val
            count_bulk_density += 1
        end
        if isfinite(base_mass_val)
            sum_base_mass += base_mass_val
            count_base_mass += 1
        end
        if isfinite(delta_thickness_val)
            abs_delta_thickness_val = abs(delta_thickness_val)
            sum_delta_thickness += delta_thickness_val
            sum_abs_delta_thickness += abs_delta_thickness_val
            max_abs_delta_thickness = max(max_abs_delta_thickness, abs_delta_thickness_val)
            count_delta_thickness += 1
            has_delta_thickness = true
        end
        if isfinite(delta_wet_mass_val)
            abs_delta_wet_mass_val = abs(delta_wet_mass_val)
            sum_delta_wet_mass += delta_wet_mass_val
            sum_abs_delta_wet_mass += abs_delta_wet_mass_val
            max_abs_delta_wet_mass = max(max_abs_delta_wet_mass, abs_delta_wet_mass_val)
            count_delta_wet_mass += 1
            has_delta_wet_mass = true
        end
        if isfinite(delta_base_mass_val)
            abs_delta_base_mass_val = abs(delta_base_mass_val)
            sum_delta_base_mass += delta_base_mass_val
            sum_abs_delta_base_mass += abs_delta_base_mass_val
            max_abs_delta_base_mass = max(max_abs_delta_base_mass, abs_delta_base_mass_val)
            count_delta_base_mass += 1
            has_delta_base_mass = true
        end
    end
    return (
        cycle = cycle,
        mean_thickness = _metric_mean(sum_thickness, count_thickness),
        mean_wet_mass = _metric_mean(sum_wet_mass, count_wet_mass),
        mean_bulk_density = _metric_mean(sum_bulk_density, count_bulk_density),
        mean_base_mass = _metric_mean(sum_base_mass, count_base_mass),
        mean_signed_delta_thickness = _metric_mean(sum_delta_thickness, count_delta_thickness),
        mean_abs_delta_thickness = _metric_mean(sum_abs_delta_thickness, count_delta_thickness),
        max_abs_delta_thickness = _metric_max(max_abs_delta_thickness, has_delta_thickness),
        mean_signed_delta_wet_mass = _metric_mean(sum_delta_wet_mass, count_delta_wet_mass),
        mean_abs_delta_wet_mass = _metric_mean(sum_abs_delta_wet_mass, count_delta_wet_mass),
        max_abs_delta_wet_mass = _metric_max(max_abs_delta_wet_mass, has_delta_wet_mass),
        mean_signed_delta_base_mass = _metric_mean(sum_delta_base_mass, count_delta_base_mass),
        mean_abs_delta_base_mass = _metric_mean(sum_abs_delta_base_mass, count_delta_base_mass),
        max_abs_delta_base_mass = _metric_max(max_abs_delta_base_mass, has_delta_base_mass),
    )
end

function make_cycle_record_and_deltas!(
    workspace::CycleMetricsWorkspace,
    cycle::Int,
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
    fill!(workspace.device_buffer, 0.0)
    n = length(thickness)
    threads = min(256, n)
    blocks = cld(n, threads)
    CUDA.@cuda threads=threads blocks=blocks _cycle_metrics_kernel!(
        workspace.device_buffer,
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
    copyto!(workspace.host_buffer, workspace.device_buffer)
    return cycle_record_from_metrics(cycle, workspace.host_buffer, n)
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
            domain;
            backend=backend,
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
            domain;
            backend=backend,
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
            domain;
            backend=backend,
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
            domain;
            backend=backend,
        )
    end
    return summary
end

function cpu_cycle_summary(summary)
    return (
        thickness = Array(summary.thickness),
        wet_mass = Array(summary.wet_mass),
        bulk_density = Array(summary.bulk_density),
        base_mass = Array(summary.base_mass),
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
        write_output = write_output,
        output_slot = output_slot,
        source_indices = source_indices,
        source_codes = source_codes,
        years = years,
    )
end

function summarize_column_state(domain::SM.AbstractSnowpackDomain, idx::Int)
    n = domain.N[idx]
    thickness = 0.0
    wet_mass = 0.0
    solid_mass = 0.0
    liquid_water = 0.0
    @inbounds for layer_index in 1:n
        solid = domain.mass[layer_index, idx]
        liquid = domain.mass_w[layer_index, idx]
        rho = domain.density[layer_index, idx]
        solid_mass += solid
        wet_mass += solid + liquid
        liquid_water += liquid
        if solid > 0.0 && rho > SM.EPS_TINY
            thickness += solid / rho
        end
    end
    return (
        thickness = thickness,
        wet_mass = wet_mass,
        bulk_density = thickness > SM.EPS_TINY ? solid_mass / thickness : 0.0,
        base_mass = domain.mass_base[idx],
        smb_ice = domain.smb_ice[idx],
        liquid_water = liquid_water,
        runoff = domain.runoff[idx],
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

function write_case_history_csv(out_path::AbstractString, history::Vector{NamedTuple})
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
    return
end

function write_case_summary(
    out_path::AbstractString,
    options::RunConfig,
    time_values::Vector{DateTime},
    ncol::Int,
    history::Vector{NamedTuple},
    status::Symbol,
    timings::TimingStats,
)
    last_record = history[end]
    cycles_completed = completed_cycle_count(history, status, options.cycles)
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
        println(io, "Cycle metrics      : ", cycle_metrics_schedule_label(options.history_stride))
        println(io, "Status             : ", string(status))
        println(io, "Cycles completed   : ", cycles_completed)
        println(io)
        println(io, "Final domain means")
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
        if status == :cycles
            println(io)
            println(io, "Interpretation     : Requested cycles completed.")
        else
            println(io)
            println(io, "Interpretation     : Run completed.")
        end
        println(io)
        print_timing_summary(io, timings)
    end
    return
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

function maybe_define_nc_output_variable!(
    vars::Dict{Symbol, Cint},
    selected::Set{Symbol},
    ncid::Cint,
    dims::Vector{Cint},
    key::Symbol,
    name::AbstractString,
    long_name::AbstractString,
    units::AbstractString,
)
    if key in selected
        vars[key] = define_nc_output_variable(ncid, dims, name, long_name, units)
    end
    return
end

function maybe_define_nc_int_variable!(
    vars::Dict{Symbol, Cint},
    selected::Set{Symbol},
    ncid::Cint,
    dims::Vector{Cint},
    key::Symbol,
    name::AbstractString,
    long_name::AbstractString,
)
    if key in selected
        vars[key] = nc_def_var(ncid, name, NC_INT, dims)
        nc_put_att_text(ncid, vars[key], "long_name", long_name)
    end
    return
end

function init_case_netcdf(
    netcdf_path::AbstractString,
    options::RunConfig,
    time_values::Vector{DateTime},
    nlayer::Int,
    layout::GridLayout,
    initial_thickness::Matrix{Float64},
    month_cycle::Vector{Int32},
    month_of_year::Vector{Int32},
    source_month_code::Vector{Int32},
    annual_output_source_indices::Vector{Int32},
    annual_output_source_codes::Vector{Int32},
)
    mkpath(dirname(netcdf_path))
    isfile(netcdf_path) && rm(netcdf_path, force=true)
    ny, nx = _grid_shape(layout)
    max_steps = options.cycles * length(annual_output_source_indices)
    selected = Set(options.netcdf_variables)
    step_cycle = Vector{Int32}(undef, max_steps)
    step_source_index = Vector{Int32}(undef, max_steps)
    step_source_code = Vector{Int32}(undef, max_steps)
    step_counter = 0
    for cyc in 1:options.cycles
        for annual_idx in eachindex(annual_output_source_indices)
            step_counter += 1
            step_cycle[step_counter] = Int32(cyc)
            step_source_index[step_counter] = annual_output_source_indices[annual_idx]
            step_source_code[step_counter] = annual_output_source_codes[annual_idx]
        end
    end

    ncid = nc_create(netcdf_path)
    dim_y = nc_def_dim(ncid, "y", ny)
    dim_x = nc_def_dim(ncid, "x", nx)
    dim_layer = nc_def_dim(ncid, "layer", max(nlayer, 1))
    dim_cycle = nc_def_dim(ncid, "cycle", options.cycles)
    dim_month = nc_def_dim(ncid, "month", length(month_cycle))
    dim_point = nc_def_dim(ncid, "point", length(layout.js))
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
    nc_put_att_text(ncid, var_point, "long_name", "Compact valid cell index")
    var_point_j = nc_def_var(ncid, "point_j", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point_j, "long_name", "1-based y-index for each compact valid cell")
    var_point_i = nc_def_var(ncid, "point_i", NC_INT, dims_p)
    nc_put_att_text(ncid, var_point_i, "long_name", "1-based x-index for each compact valid cell")
    var_point_y = nc_def_var(ncid, "point_y_km", NC_DOUBLE, dims_p)
    nc_put_att_text(ncid, var_point_y, "long_name", "Y coordinate for each compact valid cell")
    nc_put_att_text(ncid, var_point_y, "units", "km")
    var_point_x = nc_def_var(ncid, "point_x_km", NC_DOUBLE, dims_p)
    nc_put_att_text(ncid, var_point_x, "long_name", "X coordinate for each compact valid cell")
    nc_put_att_text(ncid, var_point_x, "units", "km")
    var_step = nc_def_var(ncid, "step", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step, "long_name", "Sequential yearly output index across repeated annual cycles")
    var_step_cycle = nc_def_var(ncid, "step_cycle", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_cycle, "long_name", "Repeated annual forcing cycle index for each yearly output")
    var_step_source_index = nc_def_var(ncid, "step_source_index", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_index, "long_name", "1-based index of the last forcing step included in each yearly output")
    var_step_source_code = nc_def_var(ncid, "step_source_code", NC_INT, Cint[dim_step])
    nc_put_att_text(ncid, var_step_source_code, "long_name", "Source forcing timestamp code YYYYMMDDHH for the final step included in each yearly output")

    vars = Dict{Symbol, Cint}()
    var_mask = define_nc_output_variable(ncid, dims_yx, "domain_mask", "Domain mask", "1")
    var_init_th = define_nc_output_variable(ncid, dims_yx, "initial_thickness", "Initial snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_thickness, "final_thickness", "Final snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_wet_mass, "final_wet_mass", "Final snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_bulk_density, "final_bulk_density", "Final bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_base_mass, "final_base_mass", "Cumulative firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_ice_sheet_smb, "final_ice_sheet_smb", "Cumulative net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :final_runoff, "final_runoff", "Final cumulative runoff", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_thickness, "last_cycle_delta_thickness", "Last cycle snow-thickness change", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_wet_mass, "last_cycle_delta_wet_mass", "Last cycle wet-mass change", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_base_mass, "last_cycle_delta_base_mass", "Last cycle firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_yx, :last_cycle_delta_ice_sheet_smb, "last_cycle_delta_ice_sheet_smb", "Last cycle net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_int_variable!(vars, selected, ncid, dims_yx, :n_active, "n_active", "Number of active Chion layers")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_density, "layer_density", "Final Chion layer density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_thickness, "layer_thickness", "Final Chion layer thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_snow_mass, "layer_snow_mass", "Final Chion layer snow mass", "kg m-2")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_liquid_mass, "layer_liquid_mass", "Final Chion layer liquid-water mass", "kg m-2")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_lyx, :layer_temperature_c, "layer_temperature_c", "Final Chion layer temperature", "C")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_thickness, "history_mean_thickness", "Cycle-mean snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_wet_mass, "history_mean_wet_mass", "Cycle-mean snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_bulk_density, "history_mean_bulk_density", "Cycle-mean bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_base_mass, "history_mean_base_mass", "Cycle-mean firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_abs_delta_thickness, "history_mean_abs_delta_thickness", "Cycle mean absolute snow-thickness change", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_abs_delta_wet_mass, "history_mean_abs_delta_wet_mass", "Cycle mean absolute wet-mass change", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_c, :history_mean_abs_delta_base_mass, "history_mean_abs_delta_base_mass", "Cycle mean absolute firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_thickness, "monthly_mean_thickness", "Monthly mean snow thickness", "m")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_wet_mass, "monthly_mean_wet_mass", "Monthly mean snow wet mass", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_bulk_density, "monthly_mean_bulk_density", "Monthly mean bulk snow density", "kg m-3")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_base_mass, "monthly_mean_base_mass", "Monthly mean cumulative firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_mean_ice_sheet_smb, "monthly_mean_ice_sheet_smb", "Monthly net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_export_to_ice, "monthly_export_to_ice", "Monthly firn mass exported to the ice model", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_net_ice_sheet_forcing, "monthly_net_ice_sheet_forcing", "Monthly net mass forcing to the ice sheet", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_myx, :monthly_runoff, "monthly_runoff", "Monthly runoff production", "mmWE")
    if any(v -> v in selected, (:step_export_to_ice, :step_ice_sheet_smb))
        vars[:step_valid] = nc_def_var(ncid, "step_valid", NC_INT, Cint[dim_step])
        nc_put_att_text(ncid, vars[:step_valid], "long_name", "1 where a yearly output record was completed and written, 0 for unused trailing slots")
    end
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_syx, :step_export_to_ice, "step_export_to_ice", "Annual firn mass exported to the ice model for each written output interval", "mmWE")
    maybe_define_nc_output_variable!(vars, selected, ncid, dims_syx, :step_ice_sheet_smb, "step_ice_sheet_smb", "Annual net mass forcing to the ice sheet for each written output interval", "mmWE")

    nc_put_att_text(ncid, NC_GLOBAL, "title", options.name)
    nc_put_att_text(ncid, NC_GLOBAL, "source_model", "Chion")
    nc_put_att_text(ncid, NC_GLOBAL, "input_label", isempty(options.input_label) ? "not provided" : options.input_label)
    nc_put_att_text(ncid, NC_GLOBAL, "forcing_start", string(first(time_values)))
    nc_put_att_text(ncid, NC_GLOBAL, "forcing_end", string(last(time_values)))
    nc_put_att_text(ncid, NC_GLOBAL, "cycles_completed", "pending")
    nc_put_att_text(ncid, NC_GLOBAL, "status", "pending")
    nc_put_att_text(ncid, NC_GLOBAL, "created", string(now()))

    nc_enddef(ncid)
    nc_put_var_double(ncid, var_x, layout.x)
    nc_put_var_double(ncid, var_y, layout.y)
    nc_put_var_int(ncid, var_layer, Int32.(collect(1:max(nlayer, 1))))
    nc_put_var_int(ncid, var_cycle, Int32.(collect(1:options.cycles)))
    nc_put_var_int(ncid, var_month, Int32.(collect(1:length(month_cycle))))
    nc_put_var_int(ncid, var_month_cycle, month_cycle)
    nc_put_var_int(ncid, var_month_of_year, month_of_year)
    nc_put_var_int(ncid, var_source_month_code, source_month_code)
    nc_put_var_int(ncid, var_point, Int32.(collect(1:length(layout.js))))
    nc_put_var_int(ncid, var_point_j, Int32.(layout.js))
    nc_put_var_int(ncid, var_point_i, Int32.(layout.is))
    nc_put_var_double(ncid, var_point_y, layout.y[layout.js])
    nc_put_var_double(ncid, var_point_x, layout.x[layout.is])
    nc_put_var_int(ncid, var_step, Int32.(collect(1:max_steps)))
    nc_put_var_int(ncid, var_step_cycle, step_cycle)
    nc_put_var_int(ncid, var_step_source_index, step_source_index)
    nc_put_var_int(ncid, var_step_source_code, step_source_code)
    nc_put_var_float_2d(ncid, var_mask, layout.mask)
    nc_put_var_float_2d(ncid, var_init_th, initial_thickness)
    return CaseNetCDFWriter(ncid, vars, max_steps, options.cycles)
end

function maybe_write_step_output!(writer::CaseNetCDFWriter, step_index::Int, key::Symbol, data::AbstractMatrix{<:Real})
    if haskey(writer.vars, key)
        nc_put_vara_float_3d_step_yx(writer.ncid, writer.vars[key], step_index, data)
    end
    return
end

function maybe_write_output!(writer::CaseNetCDFWriter, key::Symbol, data, writer_fn)
    if haskey(writer.vars, key)
        writer_fn(writer.ncid, writer.vars[key], data)
    end
    return
end

function finalize_case_netcdf!(
    writer::CaseNetCDFWriter,
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
    hist_th = fill(NaN, writer.cycles)
    hist_wet = fill(NaN, writer.cycles)
    hist_rho = fill(NaN, writer.cycles)
    hist_base = fill(NaN, writer.cycles)
    hist_dth = fill(NaN, writer.cycles)
    hist_dswe = fill(NaN, writer.cycles)
    hist_dbase = fill(NaN, writer.cycles)
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
    maybe_write_output!(writer, :final_thickness, final_thickness, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_wet_mass, final_wet_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_bulk_density, final_bulk_density, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_base_mass, final_base_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_ice_sheet_smb, final_ice_sheet_smb, nc_put_var_float_2d)
    maybe_write_output!(writer, :final_runoff, final_runoff, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_thickness, last_delta_thickness, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_wet_mass, last_delta_wet_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_base_mass, last_delta_base_mass, nc_put_var_float_2d)
    maybe_write_output!(writer, :last_cycle_delta_ice_sheet_smb, last_delta_ice_sheet_smb, nc_put_var_float_2d)
    maybe_write_output!(writer, :n_active, layer_grids.n_active, nc_put_var_int_2d)
    maybe_write_output!(writer, :layer_density, layer_grids.layer_density, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_thickness, layer_grids.layer_thickness, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_snow_mass, layer_grids.layer_snow_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_liquid_mass, layer_grids.layer_liquid_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :layer_temperature_c, layer_grids.layer_temperature_c, nc_put_var_float_3d)
    maybe_write_output!(writer, :history_mean_thickness, hist_th, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_wet_mass, hist_wet, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_bulk_density, hist_rho, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_base_mass, hist_base, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_abs_delta_thickness, hist_dth, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_abs_delta_wet_mass, hist_dswe, nc_put_var_float_1d)
    maybe_write_output!(writer, :history_mean_abs_delta_base_mass, hist_dbase, nc_put_var_float_1d)
    maybe_write_output!(writer, :monthly_mean_thickness, monthly_mean_thickness, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_wet_mass, monthly_mean_wet_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_bulk_density, monthly_mean_bulk_density, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_base_mass, monthly_mean_base_mass, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_mean_ice_sheet_smb, monthly_mean_ice_sheet_smb, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_export_to_ice, monthly_export_to_ice, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_net_ice_sheet_forcing, monthly_net_ice_sheet_forcing, nc_put_var_float_3d)
    maybe_write_output!(writer, :monthly_runoff, monthly_runoff, nc_put_var_float_3d)
    if haskey(writer.vars, :step_valid)
        step_valid = zeros(Int32, writer.max_steps)
        step_valid[1:steps_written] .= 1
        nc_put_var_int(writer.ncid, writer.vars[:step_valid], step_valid)
    end
    nc_redef(writer.ncid)
    nc_put_att_text(writer.ncid, NC_GLOBAL, "cycles_completed", string(cycles_completed))
    nc_put_att_text(writer.ncid, NC_GLOBAL, "status", string(status))
    nc_put_att_text(writer.ncid, NC_GLOBAL, "steps_written", string(steps_written))
    nc_enddef(writer.ncid)
    nc_close(writer.ncid)
    return
end

function run_case_cycles_threads_no_netcdf!(
    timings::TimingStats,
    options::RunConfig,
    domain::SM.SnowpackDomain,
    workspaces::AbstractVector{<:SM.StepWorkspace},
    step_fields::SM.SnowpackStepFields,
)
    options.backend == :threads || error("Threaded fast path requires `backend=:threads`.")
    !options.write_netcdf || error("Threaded fast path only applies when NetCDF output is disabled.")
    ntime = size(step_fields.air_temperature, 2)
    ncol = SM.column_count(domain)
    prev = allocate_cycle_summary_buffers(ncol)
    final = allocate_cycle_summary_buffers(ncol)
    time_block!(timings, :summarize_columns_initial) do
        summarize_cycle_columns!(prev, domain; backend=:threads)
    end
    history = NamedTuple[]
    status = :cycles
    last_delta_thickness_vec = fill(NaN, ncol)
    last_delta_wet_mass_vec = fill(NaN, ncol)
    last_delta_base_mass_vec = fill(NaN, ncol)
    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        time_counted_block!(timings, :model_step_wall, ncol * ntime) do
            SM.step!(domain, step_fields, workspaces)
        end
        time_block!(timings, :summarize_columns_cycle) do
            summarize_cycle_columns!(final, domain; backend=:threads)
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics) do
                make_cycle_record_and_deltas!(
                    cycle,
                    last_delta_thickness_vec,
                    last_delta_wet_mass_vec,
                    last_delta_base_mass_vec,
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    prev.thickness,
                    prev.wet_mass,
                    prev.base_mass,
                )
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(cycle_log_line(record))
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    return (history=history, status=status, simulation_wall_sec=simulation_wall_sec)
end

function run_case_cycles_gpu_no_netcdf!(
    timings::TimingStats,
    options::RunConfig,
    domain::SM.SnowpackDomain,
    workspace::SM.ColumnarStepWorkspace,
    step_fields::SM.SnowpackStepFields,
)
    options.backend == :gpu || error("GPU fast path requires `backend=:gpu`.")
    !options.write_netcdf || error("GPU fast path only applies when NetCDF output is disabled.")
    ntime = size(step_fields.air_temperature, 2)
    ncol = SM.column_count(domain)
    prev = allocate_cycle_summary_buffers(domain, ncol)
    final = allocate_cycle_summary_buffers(domain, ncol)
    time_block!(timings, :summarize_columns_initial; synchronize=CUDA.synchronize) do
        SM.summarize_cycle_state!(
            prev.thickness,
            prev.wet_mass,
            prev.bulk_density,
            prev.base_mass,
            domain;
            backend=:kernelabstractions,
        )
    end
    history = NamedTuple[]
    status = :cycles
    last_delta_thickness_vec = similar(domain.mass, Float64, ncol)
    last_delta_wet_mass_vec = similar(domain.mass, Float64, ncol)
    last_delta_base_mass_vec = similar(domain.mass, Float64, ncol)
    cycle_metrics_workspace = CycleMetricsWorkspace(domain)
    simulation_wall_t0 = time_ns()
    for cycle in 1:options.cycles
        time_counted_block!(timings, :model_step_wall, ncol * ntime; synchronize=CUDA.synchronize) do
            SM.step!(domain, step_fields, workspace; update_snow_cover=false)
        end
        time_block!(timings, :summarize_columns_cycle; synchronize=CUDA.synchronize) do
            SM.summarize_cycle_state!(
                final.thickness,
                final.wet_mass,
                final.bulk_density,
                final.base_mass,
                domain;
                backend=:kernelabstractions,
            )
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = time_block!(timings, :cycle_metrics; synchronize=CUDA.synchronize) do
                make_cycle_record_and_deltas!(
                    cycle_metrics_workspace,
                    cycle,
                    last_delta_thickness_vec,
                    last_delta_wet_mass_vec,
                    last_delta_base_mass_vec,
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    prev.thickness,
                    prev.wet_mass,
                    prev.base_mass,
                )
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(cycle_log_line(record))
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    return (history=history, status=status, simulation_wall_sec=simulation_wall_sec)
end

@inline function _empty_layer_grids()
    return (
        n_active = Matrix{Int32}(undef, 0, 0),
        layer_density = Array{Float64}(undef, 0, 0, 0),
        layer_thickness = Array{Float64}(undef, 0, 0, 0),
        layer_snow_mass = Array{Float64}(undef, 0, 0, 0),
        layer_liquid_mass = Array{Float64}(undef, 0, 0, 0),
        layer_temperature_c = Array{Float64}(undef, 0, 0, 0),
    )
end

@inline _empty_monthly_grid() = Array{Float64}(undef, 0, 0, 0)

function _print_run_report(
    io::IO,
    options::RunConfig,
    time_values::Vector{DateTime},
    history::Vector{NamedTuple},
    status::Symbol,
    simulation_wall_sec::Float64,
    run_wall_sec::Float64,
    timings::TimingStats;
    nc_path::AbstractString="",
    summary_path::AbstractString="",
    history_csv_path::AbstractString="",
)
    cycles_completed = completed_cycle_count(history, status, options.cycles)
    println(io, "$(options.name) complete.")
    println(io, "Input label     : ", isempty(options.input_label) ? "(not provided)" : options.input_label)
    println(io, "Forcing start   : $(first(time_values))")
    println(io, "Forcing end     : $(last(time_values))")
    println(io, "Backend         : $(String(options.backend))")
    println(io, "Cycles          : $(cycles_completed)")
    println(io, "Status          : $(string(status))")
    println(io, "Cycle metrics   : $(cycle_metrics_schedule_label(options.history_stride))")
    println(io, @sprintf("Simulation wall : %.3f s", simulation_wall_sec))
    println(io, @sprintf("Run wall total  : %.3f s", run_wall_sec))
    if haskey(timings.totals, :model_step_wall)
        println(io, @sprintf("Model step wall : %.3f s", timings.totals[:model_step_wall]))
    end
    if options.write_netcdf
        println(io, "Output NetCDF   : $(abspath(nc_path))")
    else
        println(io, "Output NetCDF   : skipped (--no-nc)")
    end
    if options.write_outputs
        println(io, "History CSV     : $(abspath(history_csv_path))")
        println(io, "Summary         : $(abspath(summary_path))")
    else
        println(io, "File outputs    : skipped (--no-output)")
    end
    print_timing_summary(io, timings; total_wall_sec=run_wall_sec)
    return
end

function execute_case!(
    domain::SM.SnowpackDomain,
    forcing::ForcingData;
    layout::Union{Nothing, GridLayout}=nothing,
    options::RunConfig=RunConfig(),
    io::IO=stdout,
    timings::TimingStats=TimingStats(),
    run_wall_t0::Integer=time_ns(),
)
    ncol = SM.column_count(domain)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the domain column count.")
    options.write_netcdf && isnothing(layout) && error("NetCDF output requires a grid layout.")
    !isnothing(layout) && length(layout.js) == ncol || isnothing(layout) || error("Grid-layout point count must match the domain column count.")

    write_final_fields = options.write_netcdf && !isnothing(layout)
    selected = Set(options.netcdf_variables)
    need_step_outputs = options.write_netcdf && case_selected(options, :step)
    need_monthly_outputs = options.write_netcdf && case_selected(options, :monthly)
    need_layer_outputs = options.write_netcdf && case_selected(options, :layers)
    need_last_cycle_smb_delta = options.write_netcdf && :last_cycle_delta_ice_sheet_smb in selected
    need_netcdf_step_diagnostics = need_step_outputs || need_monthly_outputs

    annual_output, unique_month_keys, step_month, nmonth_per_cycle, nmonth_total, month_cycle, month_of_year, source_month_code =
        time_block!(timings, :prepare_output_schedule) do
            annual_output_local = options.write_netcdf ? build_annual_output_schedule(forcing.time_values) : nothing
            unique_month_keys_local = options.write_netcdf ? unique((year(t), month(t)) for t in forcing.time_values) : Tuple{Int, Int}[]
            month_lookup = Dict{Tuple{Int, Int}, Int}()
            if options.write_netcdf
                for (idx, key) in enumerate(unique_month_keys_local)
                    month_lookup[key] = idx
                end
            end
            step_month_local = options.write_netcdf ? [month_lookup[(year(t), month(t))] for t in forcing.time_values] : Int[]
            nmonth_per_cycle_local = options.write_netcdf ? length(unique_month_keys_local) : 0
            nmonth_total_local = options.write_netcdf ? options.cycles * nmonth_per_cycle_local : 0
            month_cycle_local = Int32[]
            month_of_year_local = Int32[]
            source_month_code_local = Int32[]
            if options.write_netcdf
                for cyc in 1:options.cycles, key in unique_month_keys_local
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

    initial_thickness_vec = write_final_fields ? Vector{Float64}(undef, ncol) : Float64[]
    if write_final_fields
        time_block!(timings, :prepare_initial_output_fields) do
            @threads :static for idx in 1:ncol
                summary = summarize_column_state(domain, idx)
                initial_thickness_vec[idx] = summary.thickness
            end
        end
    end

    if options.backend == :threads && !options.write_outputs && !options.write_netcdf
        step_fields = SM.SnowpackStepFields(forcing)
        workspaces = time_block!(timings, :create_workspaces) do
            threaded_workspaces(domain)
        end
        sim = redirect_stdout(io) do
            run_case_cycles_threads_no_netcdf!(timings, options, domain, workspaces, step_fields)
        end
        run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
        _print_run_report(io, options, forcing.time_values, sim.history, sim.status, sim.simulation_wall_sec, run_wall_sec, timings)
        return RunResult(sim.history, sim.status, timings, sim.simulation_wall_sec, run_wall_sec, "", "", "", domain, options)
    end

    step_fields = SM.SnowpackStepFields(forcing)
    if options.backend == :gpu
        SM.cuda_available() || error("`backend=gpu` requested, but CUDA is not functional in the current environment.")
        domain = time_block!(timings, :gpu_transfer) do
            SM.gpu_domain(domain)
        end
        step_fields = time_block!(timings, :gpu_transfer) do
            SM.adapt(CUDA.CuArray, step_fields)
        end
        workspaces = time_block!(timings, :gpu_transfer) do
            ColumnarStepWorkspace(domain)
        end
        if !options.write_outputs && !options.write_netcdf
            sim = redirect_stdout(io) do
                run_case_cycles_gpu_no_netcdf!(
                    timings,
                    options,
                    domain,
                    workspaces,
                    step_fields,
                )
            end
            run_wall_sec = (time_ns() - run_wall_t0) * 1.0e-9
            _print_run_report(io, options, forcing.time_values, sim.history, sim.status, sim.simulation_wall_sec, run_wall_sec, timings)
            return RunResult(sim.history, sim.status, timings, sim.simulation_wall_sec, run_wall_sec, "", "", "", domain, options)
        end
    else
        workspaces = time_block!(timings, :create_workspaces) do
            threaded_workspaces(domain)
        end
    end

    gpu_netcdf_diagnostics = options.backend == :gpu && need_netcdf_step_diagnostics
    gpu_stage_sync = options.backend == :gpu ? CUDA.synchronize : nothing
    prev, final, last_delta_thickness_vec, last_delta_wet_mass_vec, last_delta_base_mass_vec, last_delta_ice_sheet_smb_vec =
        time_block!(timings, :allocate_cycle_buffers) do
            prev_local = options.backend == :gpu ? allocate_cycle_summary_buffers(domain, ncol) : allocate_cycle_summary_buffers(ncol)
            final_local = options.backend == :gpu ? allocate_cycle_summary_buffers(domain, ncol) : allocate_cycle_summary_buffers(ncol)
            last_delta_thickness_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            last_delta_wet_mass_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            last_delta_base_mass_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            last_delta_ice_sheet_smb_vec_local = options.backend == :gpu ? similar(domain.mass, Float64, ncol) : fill(NaN, ncol)
            return (
                prev_local,
                final_local,
                last_delta_thickness_vec_local,
                last_delta_wet_mass_vec_local,
                last_delta_base_mass_vec_local,
                last_delta_ice_sheet_smb_vec_local,
            )
        end
    cycle_metrics_workspace = options.backend == :gpu ? CycleMetricsWorkspace(domain) : nothing
    time_block!(timings, :summarize_columns_initial; synchronize=gpu_stage_sync) do
        if options.backend == :gpu
            SM.summarize_cycle_state!(
                prev.thickness,
                prev.wet_mass,
                prev.bulk_density,
                prev.base_mass,
                domain;
                backend=:kernelabstractions,
            )
        else
            summarize_cycle_columns!(prev, domain; backend=:threads)
        end
    end
    previous_base_mass_vec, previous_smb_ice_vec, previous_runoff_vec, previous_cycle_smb_ice_vec =
        time_block!(timings, :initialize_cycle_tracking) do
            previous_base_mass_vec_local = if gpu_netcdf_diagnostics
                copy(prev.base_mass)
            elseif need_netcdf_step_diagnostics
                copy(prev.base_mass)
            else
                Float64[]
            end
            previous_smb_ice_vec_local = if gpu_netcdf_diagnostics
                copy(domain.smb_ice)
            elseif need_netcdf_step_diagnostics
                copy(domain.smb_ice)
            else
                Float64[]
            end
            previous_runoff_vec_local = if gpu_netcdf_diagnostics
                copy(domain.runoff)
            elseif need_netcdf_step_diagnostics
                copy(domain.runoff)
            else
                Float64[]
            end
            previous_cycle_smb_ice_vec_local = need_last_cycle_smb_delta ? copy(domain.smb_ice) : Float64[]
            return (
                previous_base_mass_vec_local,
                previous_smb_ice_vec_local,
                previous_runoff_vec_local,
                previous_cycle_smb_ice_vec_local,
            )
        end

    initial_thickness = write_final_fields ? time_block!(timings, :prepare_initial_output_fields_grid) do
        scatter_to_grid(initial_thickness_vec, layout.js, layout.is, _grid_shape(layout))
    end : Matrix{Float64}(undef, 0, 0)
    nc_path = resolve_case_netcdf_path(options)
    writer = options.write_netcdf ? time_block!(timings, :init_netcdf) do
        init_case_netcdf(
            nc_path,
            options,
            forcing.time_values,
            domain.Ntot,
            layout,
            initial_thickness,
            month_cycle,
            month_of_year,
            source_month_code,
            annual_output.source_indices,
            annual_output.source_codes,
        )
    end : nothing
    write_step_fields = options.write_netcdf && !isnothing(writer) &&
        (haskey(writer.vars, :step_export_to_ice) || haskey(writer.vars, :step_ice_sheet_smb))
    step_summary, device_step_summary, annual_export_to_ice, annual_ice_sheet_smb, steps_written, history, status,
    monthly_sum_thickness, monthly_sum_wet_mass, monthly_sum_bulk_density, monthly_sum_base_mass,
    monthly_sum_ice_sheet_smb, monthly_sum_export, monthly_sum_net_ice_sheet_forcing, monthly_sum_runoff, monthly_count =
        time_block!(timings, :allocate_output_buffers) do
            step_summary_local = need_netcdf_step_diagnostics && !gpu_netcdf_diagnostics ? allocate_summary_buffers(ncol) : nothing
            device_step_summary_local = gpu_netcdf_diagnostics ? allocate_summary_buffers(domain, ncol) : nothing
            annual_export_to_ice_local = if gpu_netcdf_diagnostics && need_step_outputs
                CUDA.zeros(Float64, ncol)
            elseif need_step_outputs
                zeros(Float64, ncol)
            else
                Float64[]
            end
            annual_ice_sheet_smb_local = if gpu_netcdf_diagnostics && need_step_outputs
                CUDA.zeros(Float64, ncol)
            elseif need_step_outputs
                zeros(Float64, ncol)
            else
                Float64[]
            end
            history_local = NamedTuple[]
            status_local = :cycles
            monthly_sum_thickness_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_wet_mass_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_bulk_density_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_base_mass_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_ice_sheet_smb_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_export_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_net_ice_sheet_forcing_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_sum_runoff_local = if gpu_netcdf_diagnostics && need_monthly_outputs
                CUDA.zeros(Float64, nmonth_total, ncol)
            elseif need_monthly_outputs
                zeros(Float64, nmonth_total, ncol)
            else
                Matrix{Float64}(undef, 0, 0)
            end
            monthly_count_local = need_monthly_outputs ? zeros(Int32, nmonth_total) : Int32[]
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
    ntime = length(forcing.time_values)
    for cycle in 1:options.cycles
        for t in 1:ntime
            month_idx = need_monthly_outputs ? (cycle - 1) * nmonth_per_cycle + step_month[t] : 0
            if options.backend == :gpu
                time_counted_block!(timings, :model_step_wall, ncol; synchronize=gpu_stage_sync) do
                    SM.step!(domain, step_fields, t, workspaces)
                end
                if need_netcdf_step_diagnostics
                    time_counted_block!(timings, :step_diagnostics, ncol; synchronize=gpu_stage_sync) do
                        if gpu_netcdf_diagnostics
                            SM.summarize_domain_state!(
                                device_step_summary.thickness,
                                device_step_summary.wet_mass,
                                device_step_summary.bulk_density,
                                device_step_summary.base_mass,
                                device_step_summary.smb_ice,
                                device_step_summary.liquid_water,
                                device_step_summary.runoff,
                                domain;
                                backend=:kernelabstractions,
                            )
                            current_base = device_step_summary.base_mass
                            current_smb_ice = device_step_summary.smb_ice
                            current_runoff = device_step_summary.runoff
                            if need_monthly_outputs
                                monthly_sum_thickness[month_idx, :] .+= device_step_summary.thickness
                                monthly_sum_wet_mass[month_idx, :] .+= device_step_summary.wet_mass
                                monthly_sum_bulk_density[month_idx, :] .+= device_step_summary.bulk_density
                                monthly_sum_base_mass[month_idx, :] .+= current_base
                                monthly_sum_ice_sheet_smb[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_export[month_idx, :] .+= current_base .- previous_base_mass_vec
                                monthly_sum_net_ice_sheet_forcing[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_runoff[month_idx, :] .+= current_runoff .- previous_runoff_vec
                            end
                            if need_step_outputs
                                annual_export_to_ice .+= current_base .- previous_base_mass_vec
                                annual_ice_sheet_smb .+= current_smb_ice .- previous_smb_ice_vec
                            end
                            previous_base_mass_vec .= current_base
                            previous_smb_ice_vec .= current_smb_ice
                            previous_runoff_vec .= current_runoff
                        else
                            summarize_columns!(step_summary, domain; backend=:kernelabstractions, device_summary=device_step_summary)
                            current_base = step_summary.base_mass
                            current_smb_ice = step_summary.smb_ice
                            current_runoff = step_summary.runoff
                            if need_monthly_outputs
                                monthly_sum_thickness[month_idx, :] .+= step_summary.thickness
                                monthly_sum_wet_mass[month_idx, :] .+= step_summary.wet_mass
                                monthly_sum_bulk_density[month_idx, :] .+= step_summary.bulk_density
                                monthly_sum_base_mass[month_idx, :] .+= current_base
                                monthly_sum_ice_sheet_smb[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_export[month_idx, :] .+= current_base .- previous_base_mass_vec
                                monthly_sum_net_ice_sheet_forcing[month_idx, :] .+= current_smb_ice .- previous_smb_ice_vec
                                monthly_sum_runoff[month_idx, :] .+= current_runoff .- previous_runoff_vec
                            end
                            if need_step_outputs
                                annual_export_to_ice .+= current_base .- previous_base_mass_vec
                                annual_ice_sheet_smb .+= current_smb_ice .- previous_smb_ice_vec
                            end
                            previous_base_mass_vec .= current_base
                            previous_smb_ice_vec .= current_smb_ice
                            previous_runoff_vec .= current_runoff
                        end
                    end
                end
            elseif need_netcdf_step_diagnostics
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)

                diag_t0 = time_ns()
                summarize_columns!(step_summary, domain; backend=:threads)
                @threads :static for idx in 1:ncol
                    current_base = step_summary.base_mass[idx]
                    current_smb_ice = step_summary.smb_ice[idx]
                    current_runoff = step_summary.runoff[idx]
                    if need_monthly_outputs
                        monthly_sum_thickness[month_idx, idx] += step_summary.thickness[idx]
                        monthly_sum_wet_mass[month_idx, idx] += step_summary.wet_mass[idx]
                        monthly_sum_bulk_density[month_idx, idx] += step_summary.bulk_density[idx]
                        monthly_sum_base_mass[month_idx, idx] += current_base
                        monthly_sum_ice_sheet_smb[month_idx, idx] += current_smb_ice - previous_smb_ice_vec[idx]
                        monthly_sum_export[month_idx, idx] += current_base - previous_base_mass_vec[idx]
                        monthly_sum_net_ice_sheet_forcing[month_idx, idx] += current_smb_ice - previous_smb_ice_vec[idx]
                        monthly_sum_runoff[month_idx, idx] += current_runoff - previous_runoff_vec[idx]
                    end
                    if need_step_outputs
                        annual_export_to_ice[idx] += current_base - previous_base_mass_vec[idx]
                        annual_ice_sheet_smb[idx] += current_smb_ice - previous_smb_ice_vec[idx]
                    end
                    previous_base_mass_vec[idx] = current_base
                    previous_smb_ice_vec[idx] = current_smb_ice
                    previous_runoff_vec[idx] = current_runoff
                end
                add_timing!(timings, :step_diagnostics, (time_ns() - diag_t0) * 1.0e-9, ncol)
            else
                t0 = time_ns()
                SM.step!(domain, step_fields, t, workspaces)
                add_timing!(timings, :model_step_wall, (time_ns() - t0) * 1.0e-9, ncol)
            end
            if need_monthly_outputs
                monthly_count[month_idx] += 1
            end
            if need_step_outputs && annual_output.write_output[t]
                steps_written += 1
                if write_step_fields
                    step_export_to_ice_vec, step_ice_sheet_smb_vec, step_export_to_ice, step_ice_sheet_smb =
                        time_block!(timings, :step_output_prepare) do
                            step_export_to_ice_vec_local = gpu_netcdf_diagnostics ? Array(annual_export_to_ice) : annual_export_to_ice
                            step_ice_sheet_smb_vec_local = gpu_netcdf_diagnostics ? Array(annual_ice_sheet_smb) : annual_ice_sheet_smb
                            step_export_to_ice_local = scatter_to_grid(step_export_to_ice_vec_local, layout.js, layout.is, _grid_shape(layout))
                            step_ice_sheet_smb_local = scatter_to_grid(step_ice_sheet_smb_vec_local, layout.js, layout.is, _grid_shape(layout))
                            return (
                                step_export_to_ice_vec_local,
                                step_ice_sheet_smb_vec_local,
                                step_export_to_ice_local,
                                step_ice_sheet_smb_local,
                            )
                        end
                    time_block!(timings, :step_output_write) do
                        maybe_write_step_output!(writer, steps_written, :step_export_to_ice, step_export_to_ice)
                        maybe_write_step_output!(writer, steps_written, :step_ice_sheet_smb, step_ice_sheet_smb)
                    end
                end
                fill!(annual_export_to_ice, 0.0)
                fill!(annual_ice_sheet_smb, 0.0)
            end
        end

        if options.backend == :gpu
            time_block!(timings, :summarize_columns_cycle; synchronize=gpu_stage_sync) do
                SM.summarize_cycle_state!(
                    final.thickness,
                    final.wet_mass,
                    final.bulk_density,
                    final.base_mass,
                    domain;
                    backend=:kernelabstractions,
                )
            end
        else
            time_block!(timings, :summarize_columns_cycle) do
                summarize_cycle_columns!(final, domain; backend=:threads)
            end
        end
        if should_record_cycle_metrics(cycle, options.cycles, options.history_stride)
            record = if options.backend == :gpu
                time_block!(timings, :cycle_metrics; synchronize=gpu_stage_sync) do
                    record_local = make_cycle_record_and_deltas!(
                        cycle_metrics_workspace,
                        cycle,
                        last_delta_thickness_vec,
                        last_delta_wet_mass_vec,
                        last_delta_base_mass_vec,
                        final.thickness,
                        final.wet_mass,
                        final.bulk_density,
                        final.base_mass,
                        prev.thickness,
                        prev.wet_mass,
                        prev.base_mass,
                    )
                    if need_last_cycle_smb_delta
                        last_delta_ice_sheet_smb_vec .= domain.smb_ice .- previous_cycle_smb_ice_vec
                        previous_cycle_smb_ice_vec .= domain.smb_ice
                    end
                    record_local
                end
            else
                time_block!(timings, :cycle_metrics) do
                    record_local = make_cycle_record_and_deltas!(
                        cycle,
                        last_delta_thickness_vec,
                        last_delta_wet_mass_vec,
                        last_delta_base_mass_vec,
                        final.thickness,
                        final.wet_mass,
                        final.bulk_density,
                        final.base_mass,
                        prev.thickness,
                        prev.wet_mass,
                        prev.base_mass,
                    )
                    if need_last_cycle_smb_delta
                        current_smb_ice_vec = copy(domain.smb_ice)
                        last_delta_ice_sheet_smb_vec .= current_smb_ice_vec .- previous_cycle_smb_ice_vec
                        previous_cycle_smb_ice_vec .= current_smb_ice_vec
                    end
                    record_local
                end
            end
            push!(history, record)
            time_block!(timings, :cycle_logging) do
                println(io, cycle_log_line(record))
            end
        elseif options.backend == :gpu
            time_block!(timings, :cycle_state_deltas; synchronize=gpu_stage_sync) do
                last_delta_thickness_vec .= final.thickness .- prev.thickness
                last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
                last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
                if need_last_cycle_smb_delta
                    last_delta_ice_sheet_smb_vec .= domain.smb_ice .- previous_cycle_smb_ice_vec
                    previous_cycle_smb_ice_vec .= domain.smb_ice
                end
            end
        else
            time_block!(timings, :cycle_state_deltas) do
                last_delta_thickness_vec .= final.thickness .- prev.thickness
                last_delta_wet_mass_vec .= final.wet_mass .- prev.wet_mass
                last_delta_base_mass_vec .= final.base_mass .- prev.base_mass
                if need_last_cycle_smb_delta
                    current_smb_ice_vec = copy(domain.smb_ice)
                    last_delta_ice_sheet_smb_vec .= current_smb_ice_vec .- previous_cycle_smb_ice_vec
                    previous_cycle_smb_ice_vec .= current_smb_ice_vec
                end
            end
        end
        prev, final = final, prev
    end
    simulation_wall_sec = (time_ns() - simulation_wall_t0) * 1.0e-9
    final_state = if status == :cycles
        history[end].cycle == options.cycles ? final : prev
    else
        final
    end

    final_thickness = Matrix{Float64}(undef, 0, 0)
    final_wet_mass = Matrix{Float64}(undef, 0, 0)
    final_bulk_density = Matrix{Float64}(undef, 0, 0)
    final_base_mass = Matrix{Float64}(undef, 0, 0)
    final_ice_sheet_smb = Matrix{Float64}(undef, 0, 0)
    final_runoff = Matrix{Float64}(undef, 0, 0)
    last_delta_thickness = Matrix{Float64}(undef, 0, 0)
    last_delta_wet_mass = Matrix{Float64}(undef, 0, 0)
    last_delta_base_mass = Matrix{Float64}(undef, 0, 0)
    last_delta_ice_sheet_smb = Matrix{Float64}(undef, 0, 0)
    if write_final_fields
        final_state_host, final_smb_ice_vec, final_runoff_vec, last_delta_thickness_host, last_delta_wet_mass_host,
        last_delta_base_mass_host, last_delta_ice_sheet_smb_host = time_block!(timings, :finalize_state_transfer) do
            final_state_host_local = options.backend == :gpu ? cpu_cycle_summary(final_state) : final_state
            final_smb_ice_vec_local = options.backend == :gpu ? Array(domain.smb_ice) : copy(domain.smb_ice)
            final_runoff_vec_local = options.backend == :gpu ? Array(domain.runoff) : copy(domain.runoff)
            last_delta_thickness_host_local = options.backend == :gpu ? Array(last_delta_thickness_vec) : last_delta_thickness_vec
            last_delta_wet_mass_host_local = options.backend == :gpu ? Array(last_delta_wet_mass_vec) : last_delta_wet_mass_vec
            last_delta_base_mass_host_local = options.backend == :gpu ? Array(last_delta_base_mass_vec) : last_delta_base_mass_vec
            last_delta_ice_sheet_smb_host_local = if need_last_cycle_smb_delta
                options.backend == :gpu ? Array(last_delta_ice_sheet_smb_vec) : last_delta_ice_sheet_smb_vec
            else
                fill(NaN, ncol)
            end
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
                    scatter_to_grid(final_state_host.thickness, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_state_host.wet_mass, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_state_host.bulk_density, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_state_host.base_mass, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_smb_ice_vec, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(final_runoff_vec, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_thickness_host, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_wet_mass_host, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_base_mass_host, layout.js, layout.is, _grid_shape(layout)),
                    scatter_to_grid(last_delta_ice_sheet_smb_host, layout.js, layout.is, _grid_shape(layout)),
                )
            end
    end

    layer_grids = _empty_layer_grids()
    monthly_mean_thickness_grid = _empty_monthly_grid()
    monthly_mean_wet_mass_grid = _empty_monthly_grid()
    monthly_mean_bulk_density_grid = _empty_monthly_grid()
    monthly_mean_base_mass_grid = _empty_monthly_grid()
    monthly_mean_ice_sheet_smb_grid = _empty_monthly_grid()
    monthly_export_to_ice_grid = _empty_monthly_grid()
    monthly_net_ice_sheet_forcing_grid = _empty_monthly_grid()
    monthly_runoff_grid = _empty_monthly_grid()
    if options.write_netcdf
        if need_layer_outputs
            final_domain = options.backend == :gpu ? time_block!(timings, :gpu_transfer) do
                SM.cpu_domain(domain)
            end : domain
            layer_grids = time_block!(timings, :collect_final_layer_grids) do
                collect_final_layer_grids(final_domain, layout.js, layout.is, _grid_shape(layout), final_domain.Ntot)
            end
        end
        if need_monthly_outputs
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
                    monthly_mean_base_mass[m, :] .= monthly_sum_base_mass_host[m, :]
                    monthly_mean_ice_sheet_smb[m, :] .= monthly_sum_ice_sheet_smb_host[m, :]
                    monthly_export_to_ice[m, :] .= monthly_sum_export_host[m, :]
                    monthly_net_ice_sheet_forcing[m, :] .= monthly_sum_net_ice_sheet_forcing_host[m, :]
                    monthly_runoff[m, :] .= monthly_sum_runoff_host[m, :]
                end
                return (
                    monthly_vectors_to_grids(monthly_mean_thickness, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_wet_mass, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_bulk_density, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_base_mass, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_mean_ice_sheet_smb, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_export_to_ice, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_net_ice_sheet_forcing, layout.js, layout.is, _grid_shape(layout)),
                    monthly_vectors_to_grids(monthly_runoff, layout.js, layout.is, _grid_shape(layout)),
                )
            end
        end
    end

    summary_path = ""
    history_csv_path = ""
    if options.write_outputs
        mkpath(options.output_dir)
        summary_path = joinpath(options.output_dir, "$(options.name)_summary.txt")
        history_csv_path = joinpath(options.output_dir, "$(options.name)_history.csv")
        time_block!(timings, :write_summary_text) do
            write_case_summary(summary_path, options, forcing.time_values, ncol, history, status, timings)
        end
        time_block!(timings, :write_history_csv) do
            write_case_history_csv(history_csv_path, history)
        end
    end
    if options.write_netcdf
        time_block!(timings, :write_netcdf) do
            finalize_case_netcdf!(
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
    _print_run_report(
        io,
        options,
        forcing.time_values,
        history,
        status,
        simulation_wall_sec,
        run_wall_sec,
        timings;
        nc_path=nc_path,
        summary_path=summary_path,
        history_csv_path=history_csv_path,
    )
    return RunResult(history, status, timings, simulation_wall_sec, run_wall_sec, options.write_netcdf ? nc_path : "", summary_path, history_csv_path, domain, options)
end
