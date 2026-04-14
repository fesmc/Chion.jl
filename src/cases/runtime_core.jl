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
    timings::StepTimingStats
    simulation_wall_sec::Float64
    run_wall_sec::Float64
    netcdf_path::String
    summary_path::String
    history_csv_path::String
    domain
    run::RunConfig
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
