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

mutable struct TimingStats
    totals::Dict{Symbol, Float64}
    counts::Dict{Symbol, Int}
end

TimingStats() = TimingStats(Dict{Symbol, Float64}(), Dict{Symbol, Int}())

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

struct GridLayout
    x::Vector{Float64}
    y::Vector{Float64}
    js::Vector{Int}
    is::Vector{Int}
    mask::Matrix{Float64}
end

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
    dataset::NCDataset
    vars::Dict{Symbol, Any}
    max_steps::Int
    cycles::Int
end

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
    share_total = isnothing(total_wall_sec) ? sum(values(stats.totals)) : total_wall_sec
    for key in keys(stats.totals)
        total_sec = stats.totals[key]
        count = stats.counts[key]
        push!(rows, (
            key=key,
            total_sec=total_sec,
            count=count,
            mean_sec=count > 0 ? total_sec / count : NaN,
            share_pct=share_total > 0.0 ? 100.0 * total_sec / share_total : 0.0,
        ))
    end
    sort!(rows; by=row -> row.total_sec, rev=true)
    return rows, sum(values(stats.totals))
end

function print_timing_summary(io::IO, stats::TimingStats; total_wall_sec::Union{Nothing, Float64}=nothing)
    rows, total = timing_rows(stats; total_wall_sec=total_wall_sec)
    println(io, "Timing summary")
    println(io, @sprintf("  %-24s %12s %9s %12s %10s", "stage", "total [s]", "share", "mean [ms]", "count"))
    for row in rows
        println(io, @sprintf("  %-24s %12.3f %8.1f%% %12.3f %10d", String(row.key), row.total_sec, row.share_pct, row.mean_sec * 1.0e3, row.count))
    end
    if isnothing(total_wall_sec)
        println(io, @sprintf("  %-24s %12.3f", "total_accounted", total))
        return
    end
    unaccounted = max(total_wall_sec - total, 0.0)
    println(io, @sprintf("  %-24s %12.3f %8.1f%% %12s %10s", "unaccounted", unaccounted, total_wall_sec > 0.0 ? 100.0 * unaccounted / total_wall_sec : 0.0, "", ""))
    println(io, @sprintf("  %-24s %12.3f", "total_accounted", total))
    println(io, @sprintf("  %-24s %12.3f", "run_wall_total", total_wall_sec))
end

@inline function normalize_case_backend(backend)
    value = lowercase(strip(String(backend)))
    value == "cpu" && return :threads
    value in ("threads", "gpu") || error("Unsupported backend '$backend'. Use `threads`, `cpu`, or `gpu`.")
    return Symbol(value)
end

@inline normalize_history_stride(stride::Integer) = Int(stride) >= 0 ? Int(stride) : error("`history_stride` must be >= 0.")
@inline should_record_cycle_metrics(cycle::Int, cycles::Int, stride::Int) = cycle == cycles || (stride > 0 && mod(cycle, stride) == 0)
@inline completed_cycle_count(history::Vector{NamedTuple}, ::Symbol, cycles::Int) = isempty(history) ? 0 : min(history[end].cycle, cycles)
@inline cycle_metrics_schedule_label(stride::Int) = stride == 0 ? "final cycle only" : stride == 1 ? "every cycle" : "every $(stride) cycles + final"
@inline _looks_like_directory_path(path::AbstractString) = !isempty(path) && (endswith(path, '/') || endswith(path, '\\'))

function _case_slug(name::AbstractString)
    slug = strip(replace(lowercase(strip(String(name))), r"[^a-z0-9]+" => "_"), '_')
    return isempty(slug) ? "snowpack_case" : slug
end

_default_case_output_dir(name::AbstractString) = joinpath(pwd(), "case_output", _case_slug(name))

function resolve_case_netcdf_path(options::RunConfig)
    default_name = "$(options.name)_final_state.nc"
    isempty(options.netcdf_path) && return joinpath(options.output_dir, default_name)
    return isdir(options.netcdf_path) || _looks_like_directory_path(options.netcdf_path) ?
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
    allowed_groups = String.(propertynames(CASE_OUTPUT_GROUPS))
    for token in tokens
        stripped = strip(token)
        isempty(stripped) && continue
        key = Symbol(lowercase(stripped))
        if key == :all
            append!(selected, CASE_NETCDF_VARIABLES)
        elseif key == :none
            continue
        elseif hasproperty(CASE_OUTPUT_GROUPS, key)
            append!(selected, getproperty(CASE_OUTPUT_GROUPS, key))
        elseif key in CASE_NETCDF_VARIABLES
            push!(selected, key)
        else
            error("Unsupported NetCDF variable selector '$token'. Use `all`, `none`, a group ($(join(sort!(allowed_groups), ", "))), or an explicit variable name.")
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
    resolved_name = String(name)
    resolved_output_dir = isempty(output_dir) ? _default_case_output_dir(resolved_name) : String(output_dir)
    return RunConfig(
        resolved_name,
        String(input_label),
        resolved_output_dir,
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
    size(data, 2) == ntime || error("Matrix forcing fields must have $ntime columns, got $(size(data, 2)).")
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
    has_native_inputs && has_user_inputs && error("Pass either model-native forcing fields or user-facing forcing fields, not both.")

    dt_days_v = Float64.(collect(dt_days))
    isempty(dt_days_v) && error("`dt_days` must not be empty.")
    all(>(0.0), dt_days_v) || error("All `dt_days` entries must be positive.")
    ntime = length(dt_days_v)
    column_count = isnothing(ncol) ? 1 : Int(ncol)
    if isnothing(ncol)
        for field in (air_temperature, snowfall_rate, rainfall_rate, air_temperature_c, snowfall_mm_day, rainfall_mm_day, shortwave_down, wind_speed, q_lw_down, q_sh, q_lh)
            isnothing(field) || ((column_count = _forcing_column_count(field, ntime)); break)
        end
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

    dims = size(air_temperature)
    shortwave_down_m = _forcing_numeric_matrix(shortwave_down, column_count, ntime, "shortwave_down")
    wind_speed_m = isnothing(wind_speed) ? fill(5.0, dims) : _forcing_numeric_matrix(wind_speed, column_count, ntime, "wind_speed")
    q_lw_down_m = isnothing(q_lw_down) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_lw_down, column_count, ntime, "q_lw_down")
    has_q_lw_down_m = isnothing(q_lw_down) ? fill(false, dims) : isnothing(has_q_lw_down) ? fill(true, dims) : _forcing_bool_matrix(has_q_lw_down, column_count, ntime, "has_q_lw_down")
    q_sh_m = isnothing(q_sh) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_sh, column_count, ntime, "q_sh")
    has_q_sh_m = isnothing(q_sh) ? fill(false, dims) : isnothing(has_q_sh) ? fill(true, dims) : _forcing_bool_matrix(has_q_sh, column_count, ntime, "has_q_sh")
    q_lh_m = isnothing(q_lh) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_lh, column_count, ntime, "q_lh")
    has_q_lh_m = isnothing(q_lh) ? fill(false, dims) : isnothing(has_q_lh) ? fill(true, dims) : _forcing_bool_matrix(has_q_lh, column_count, ntime, "has_q_lh")

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

    time_values_v = isnothing(time_values) ? _synthesized_time_values(dt_days_v) : DateTime.(collect(time_values))
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

function GridLayout(x::AbstractVector, y::AbstractVector, js::AbstractVector{<:Integer}, is::AbstractVector{<:Integer}, mask)
    x_v, y_v = Float64.(x), Float64.(y)
    js_v, is_v = Int.(js), Int.(is)
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
    mass_m, mass_w_m = Matrix{Float64}(mass), Matrix{Float64}(mass_w)
    density_m, temperature_m = Matrix{Float64}(density), Matrix{Float64}(temperature)
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

@inline function case_selected(options::RunConfig, group::Symbol)
    group_vars = getproperty(CASE_OUTPUT_GROUPS, group)
    return any(var -> var in group_vars, options.netcdf_variables)
end

@inline _grid_shape(layout::GridLayout) = size(layout.mask)

const SUMMARY_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass, :smb_ice, :liquid_water, :runoff)
const CYCLE_BUFFER_NAMES = (:thickness, :wet_mass, :bulk_density, :base_mass)

_named_buffers(names::NTuple{N, Symbol}, build::F) where {N, F <: Function} = NamedTuple{names}(ntuple(_ -> build(), N))

allocate_summary_buffers(n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int) = _named_buffers(SUMMARY_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))
allocate_cycle_summary_buffers(n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> Vector{Float64}(undef, n))
allocate_cycle_summary_buffers(domain::SM.AbstractSnowpackDomain, n::Int) = _named_buffers(CYCLE_BUFFER_NAMES, () -> similar(domain.mass, Float64, n))

@inline _host_vector(data::Vector{Float64}; copy_array::Bool=false) = copy_array ? copy(data) : data
@inline _host_vector(data; copy_array::Bool=false) = Float64.(Array(data))

function _finite_mean(data)
    total = 0.0
    count = 0
    for value in data
        isfinite(value) || continue
        total += value
        count += 1
    end
    return count == 0 ? NaN : total / count
end

function _delta_stats(data)
    signed_total = 0.0
    abs_total = 0.0
    max_abs = 0.0
    count = 0
    for value in data
        isfinite(value) || continue
        abs_value = abs(value)
        signed_total += value
        abs_total += abs_value
        max_abs = max(max_abs, abs_value)
        count += 1
    end
    return (
        mean_signed=count == 0 ? NaN : signed_total / count,
        mean_abs=count == 0 ? NaN : abs_total / count,
        max_abs=count == 0 ? NaN : max_abs,
    )
end

function make_cycle_record_and_deltas!(
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
    delta_thickness .= thickness .- prev_thickness
    delta_wet_mass .= wet_mass .- prev_wet_mass
    delta_base_mass .= base_mass .- prev_base_mass
    dth = _delta_stats(_host_vector(delta_thickness))
    dwet = _delta_stats(_host_vector(delta_wet_mass))
    dbase = _delta_stats(_host_vector(delta_base_mass))
    return (
        cycle=cycle,
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

function summarize_columns!(summary, domain::SM.SnowpackDomain; backend::Symbol=:threads, device_summary=nothing)
    if backend == :kernelabstractions
        isnothing(device_summary) && (device_summary = allocate_summary_buffers(domain, length(summary.thickness)))
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
        return _copy_summary_fields!(summary, device_summary, SUMMARY_BUFFER_NAMES)
    end
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
    return summary
end

function summarize_cycle_columns!(summary, domain::SM.SnowpackDomain; backend::Symbol=:threads, device_summary=nothing)
    if backend == :kernelabstractions
        isnothing(device_summary) && (device_summary = allocate_cycle_summary_buffers(domain, length(summary.thickness)))
        SM.summarize_cycle_state!(
            device_summary.thickness,
            device_summary.wet_mass,
            device_summary.bulk_density,
            device_summary.base_mass,
            domain;
            backend=backend,
        )
        return _copy_summary_fields!(summary, device_summary, CYCLE_BUFFER_NAMES)
    end
    SM.summarize_cycle_state!(
        summary.thickness,
        summary.wet_mass,
        summary.bulk_density,
        summary.base_mass,
        domain;
        backend=backend,
    )
    return summary
end

function build_annual_output_schedule(time_values::Vector{DateTime})
    write_output = falses(length(time_values))
    output_slot = zeros(Int, length(time_values))
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
    return (write_output=write_output, output_slot=output_slot, source_indices=source_indices, source_codes=source_codes, years=years)
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

const HISTORY_CSV_SPECS = (
    (key=:cycle, label="cycle", integer=true),
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

function write_case_history_csv(out_path::AbstractString, history::Vector{NamedTuple})
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        println(io, join((spec.label for spec in HISTORY_CSV_SPECS), ","))
        for rec in history
            values = String[]
            for spec in HISTORY_CSV_SPECS
                value = getfield(rec, spec.key)
                push!(values, spec.integer ? string(value) : @sprintf("%.10f", value))
            end
            println(io, join(values, ","))
        end
    end
end

const SUMMARY_REPORT_SPECS = (
    (title="Final domain means", fields=(
        ("Thickness (m)", :mean_thickness),
        ("Wet mass (mmWE)", :mean_wet_mass),
        ("Bulk density (kg m-3)", :mean_bulk_density),
        ("Firn-to-ice mass (mmWE)", :mean_base_mass),
    )),
    (title="Last cycle deltas", fields=(
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
        println(io, "Cycles completed   : ", completed_cycle_count(history, status, options.cycles))
        for section in SUMMARY_REPORT_SPECS
            println(io)
            println(io, section.title)
            for (label, key) in section.fields
                println(io, @sprintf("%-28s : %.6f", label, getfield(last_record, key)))
            end
        end
        println(io)
        println(io, "Interpretation     : Requested cycles completed.")
        println(io)
        print_timing_summary(io, timings)
    end
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
        j, i = js[idx], is[idx]
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
        n_active=n_active,
        layer_density=layer_density,
        layer_thickness=layer_thickness,
        layer_snow_mass=layer_snow_mass,
        layer_liquid_mass=layer_liquid_mass,
        layer_temperature_c=layer_temperature_c,
    )
end

empty_final_grids() = NamedTuple{FINAL_GRID_KEYS}(ntuple(_ -> Matrix{Float64}(undef, 0, 0), length(FINAL_GRID_KEYS)))
empty_monthly_grids() = NamedTuple{MONTHLY_GRID_KEYS}(ntuple(_ -> Array{Float64}(undef, 0, 0, 0), length(MONTHLY_GRID_KEYS)))

function _empty_layer_grids()
    return (
        n_active=Matrix{Int32}(undef, 0, 0),
        layer_density=Array{Float64}(undef, 0, 0, 0),
        layer_thickness=Array{Float64}(undef, 0, 0, 0),
        layer_snow_mass=Array{Float64}(undef, 0, 0, 0),
        layer_liquid_mass=Array{Float64}(undef, 0, 0, 0),
        layer_temperature_c=Array{Float64}(undef, 0, 0, 0),
    )
end
