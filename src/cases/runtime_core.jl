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

High-level execution settings for [`SnowpackCase`](@ref) and [`run_case`](@ref),
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

function _case_slug(name::AbstractString)
    slug = replace(lowercase(strip(String(name))), r"[^a-z0-9]+" => "_")
    slug = strip(slug, '_')
    return isempty(slug) ? "snowpack_case" : slug
end

_default_case_output_dir(name::AbstractString) = joinpath(pwd(), "case_output", _case_slug(name))

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
