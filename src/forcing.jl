@inline function _synthesized_time_values(dt_days::Vector{Float64})
    base = DateTime(2000, 1, 1, 12)
    out = Vector{DateTime}(undef, length(dt_days))
    elapsed_ms = 0
    for idx in eachindex(dt_days)
        out[idx] = base + Dates.Millisecond(elapsed_ms)
        elapsed_ms += round(Int, dt_days[idx] * 86_400_000)
    end
    return out
end

@inline function _ensure_matching_field_sizes(reference::Tuple{Int,Int}, name::AbstractString, field)
    size(field) == reference || error("`$name` must have shape $(reference), got $(size(field)).")
end

@inline function _forcing_column_count(field, ntime::Int)
    field isa Number && return 1
    data = collect(field)
    ndims(data) == 1 && return 1
    ndims(data) == 2 || error("Forcing fields must be scalars, vectors, or matrices.")
    size(data, 2) == ntime || error("Matrix forcing fields must have $ntime columns, got $(size(data, 2)).")
    return size(data, 1)
end

@inline function _forcing_numeric_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
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

@inline function _forcing_bool_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
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

@inline function _forcing_column_metadata_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Number
        return fill(Float64(field), ncol, ntime)
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ncol || error("`$name` vector input must have length $ncol.")
        return repeat(reshape(Float64.(data), ncol, 1), 1, ntime)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Float64}(data)
    end
    error("`$name` must be a scalar, a vector of length $ncol, or a matrix of size ($ncol, $ntime).")
end

@inline function _calendar_day_of_year(t::DateTime)
    seconds_today = hour(t) * 3600 + minute(t) * 60 + second(t) + millisecond(t) / 1000
    return Float64(dayofyear(t)) + seconds_today / 86_400.0
end

@inline function _solar_longitude_deg_from_calendar_day(day_of_year)

    ##check this again
    days_since_j2000_like_year_start = day_of_year - 1.0
    mean_longitude = 280.46646 + 0.98564736 * days_since_j2000_like_year_start
    mean_anomaly = 357.52911 + 0.98560028 * days_since_j2000_like_year_start
    true_longitude = mean_longitude +
                     1.914602 * sind(mean_anomaly) +
                     0.019993 * sind(2.0 * mean_anomaly)
    return mod(true_longitude, 360.0)
end

"""Time-varying atmospheric boundary conditions for a snowpack model."""
struct SnowpackForcing
    time_values::Vector{DateTime}
    dt_days::Vector{Float64}
    day_of_year::Vector{Float64}
    solar_longitude_deg::Vector{Float64}
    air_temperature
    snowfall_rate
    rainfall_rate
    shortwave_down
    latitude_deg
    wind_speed
    q_lw_down
    has_q_lw_down
    q_sh
    has_q_sh
    q_lh
    has_q_lh
    relative_humidity
    has_relative_humidity
    air_pressure
    prescribed_albedo
    has_prescribed_albedo
end

@adapt_structure SnowpackForcing

@inline function forcing_step_kind(forcing::SnowpackForcing, time_index::Int)
    ## for monthly pdd
    dt = _step_dt(forcing.dt_days, time_index)
    return 27.0 <= dt <= 32.0 ? :monthly : :scheduled
end

function SnowpackForcing(;
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
    relative_humidity=nothing,
    has_relative_humidity=nothing,
    air_pressure=nothing,
    prescribed_albedo=nothing,
    has_prescribed_albedo=nothing,
    latitude_deg=nothing,
    time_values=nothing,
)
    has_native = !isnothing(air_temperature) || !isnothing(snowfall_rate) || !isnothing(rainfall_rate)
    has_user = !isnothing(air_temperature_c) || !isnothing(snowfall_mm_day) || !isnothing(rainfall_mm_day)
    has_native && has_user && error("Pass either model-native forcing fields or user-facing fields, not both.")

    dt_days_v = Float64.(collect(dt_days))
    isempty(dt_days_v) && error("`dt_days` must not be empty.")
    all(>(0.0), dt_days_v) || error("All `dt_days` entries must be positive.")
    ntime = length(dt_days_v)

    column_count = isnothing(ncol) ? 1 : Int(ncol)
    if isnothing(ncol)
        if !isnothing(latitude_deg) && !(latitude_deg isa Number)
            latitude_data = collect(latitude_deg)
            ndims(latitude_data) == 1 && (column_count = length(latitude_data))
            ndims(latitude_data) == 2 && (column_count = size(latitude_data, 1))
        end
        if isnothing(latitude_deg) || latitude_deg isa Number
            for field in (air_temperature, snowfall_rate, rainfall_rate, air_temperature_c,
                snowfall_mm_day, rainfall_mm_day, shortwave_down, wind_speed, q_lw_down, q_sh, q_lh,
                relative_humidity, air_pressure, prescribed_albedo)
                isnothing(field) || ((column_count = _forcing_column_count(field, ntime)); break)
            end
        end
    end
    column_count > 0 || error("`ncol` must be positive.")

    if has_user
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
    latitude_deg_m = isnothing(latitude_deg) ? fill(NaN, dims) : _forcing_column_metadata_matrix(latitude_deg, column_count, ntime, "latitude_deg")
    q_lw_down_m = isnothing(q_lw_down) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_lw_down, column_count, ntime, "q_lw_down")
    has_q_lw_down_m = isnothing(q_lw_down) ? fill(false, dims) : isnothing(has_q_lw_down) ? fill(true, dims) : _forcing_bool_matrix(has_q_lw_down, column_count, ntime, "has_q_lw_down")
    q_sh_m = isnothing(q_sh) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_sh, column_count, ntime, "q_sh")
    has_q_sh_m = isnothing(q_sh) ? fill(false, dims) : isnothing(has_q_sh) ? fill(true, dims) : _forcing_bool_matrix(has_q_sh, column_count, ntime, "has_q_sh")
    q_lh_m = isnothing(q_lh) ? zeros(Float64, dims) : _forcing_numeric_matrix(q_lh, column_count, ntime, "q_lh")
    has_q_lh_m = isnothing(q_lh) ? fill(false, dims) : isnothing(has_q_lh) ? fill(true, dims) : _forcing_bool_matrix(has_q_lh, column_count, ntime, "has_q_lh")
    relative_humidity_m = if !isnothing(relative_humidity)
        _forcing_numeric_matrix(relative_humidity, column_count, ntime, "relative_humidity")
    else
        zeros(Float64, dims)
    end
    has_relative_humidity_m = if isnothing(relative_humidity)
        fill(false, dims)
    elseif isnothing(has_relative_humidity)
        isfinite.(relative_humidity_m)
    else
        _forcing_bool_matrix(has_relative_humidity, column_count, ntime, "has_relative_humidity")
    end
    relative_humidity_m[.!has_relative_humidity_m] .= 0.0
    air_pressure_m = isnothing(air_pressure) ? fill(101_325.0, dims) : _forcing_numeric_matrix(air_pressure, column_count, ntime, "air_pressure")
    prescribed_albedo_m = isnothing(prescribed_albedo) ? zeros(Float64, dims) : _forcing_numeric_matrix(prescribed_albedo, column_count, ntime, "prescribed_albedo")
    has_prescribed_albedo_m = isnothing(prescribed_albedo) ? fill(false, dims) : isnothing(has_prescribed_albedo) ? fill(true, dims) : _forcing_bool_matrix(has_prescribed_albedo, column_count, ntime, "has_prescribed_albedo")

    for (name, field) in (
        ("snowfall_rate", snowfall_rate),
        ("rainfall_rate", rainfall_rate),
        ("shortwave_down", shortwave_down_m),
        ("latitude_deg", latitude_deg_m),
        ("wind_speed", wind_speed_m),
        ("q_lw_down", q_lw_down_m),
        ("has_q_lw_down", has_q_lw_down_m),
        ("q_sh", q_sh_m),
        ("has_q_sh", has_q_sh_m),
        ("q_lh", q_lh_m),
        ("has_q_lh", has_q_lh_m),
        ("relative_humidity", relative_humidity_m),
        ("has_relative_humidity", has_relative_humidity_m),
        ("air_pressure", air_pressure_m),
        ("prescribed_albedo", prescribed_albedo_m),
        ("has_prescribed_albedo", has_prescribed_albedo_m),
    )
        _ensure_matching_field_sizes(dims, name, field)
    end

    time_values_v = isnothing(time_values) ? _synthesized_time_values(dt_days_v) : DateTime.(collect(time_values))
    length(time_values_v) == dims[2] || error("`time_values` must have one entry per forcing timestep.")
    day_of_year_v = _calendar_day_of_year.(time_values_v)
    solar_longitude_deg_v = _solar_longitude_deg_from_calendar_day.(day_of_year_v)

    return SnowpackForcing(
        time_values_v,
        dt_days_v,
        day_of_year_v,
        solar_longitude_deg_v,
        air_temperature,
        snowfall_rate,
        rainfall_rate,
        shortwave_down_m,
        latitude_deg_m,
        wind_speed_m,
        q_lw_down_m,
        has_q_lw_down_m,
        q_sh_m,
        has_q_sh_m,
        q_lh_m,
        has_q_lh_m,
        relative_humidity_m,
        has_relative_humidity_m,
        air_pressure_m,
        prescribed_albedo_m,
        has_prescribed_albedo_m,
    )
end
