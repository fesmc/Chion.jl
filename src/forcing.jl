"""Read-only matrix whose entries all have the same value."""
struct ConstantForcingMatrix{T} <: AbstractMatrix{T}
    value::T
    dims::Tuple{Int, Int}
end

Base.size(field::ConstantForcingMatrix) = field.dims
@inline function Base.getindex(field::ConstantForcingMatrix, column::Int, time::Int)
    @boundscheck checkbounds(field, column, time)
    return field.value
end
Base.IndexStyle(::Type{<:ConstantForcingMatrix}) = IndexCartesian()
Adapt.@adapt_structure ConstantForcingMatrix

"""Logical matrix that repeats one value per column over all timesteps."""
struct ColumnForcingMatrix{T, V <: AbstractVector{T}} <: AbstractMatrix{T}
    values::V
    ntime::Int
end

Base.size(field::ColumnForcingMatrix) = (length(field.values), field.ntime)
@inline function Base.getindex(field::ColumnForcingMatrix, column::Int, time::Int)
    @boundscheck checkbounds(field, column, time)
    return @inbounds field.values[column]
end
@inline function Base.setindex!(field::ColumnForcingMatrix, value, column::Int, time::Int)
    @boundscheck checkbounds(field, column, time)
    @inbounds field.values[column] = value
    return value
end
Base.IndexStyle(::Type{<:ColumnForcingMatrix}) = IndexCartesian()
Adapt.@adapt_structure ColumnForcingMatrix

"""Logical matrix that repeats one value per timestep over all columns."""
struct TimeForcingMatrix{T, V <: AbstractVector{T}} <: AbstractMatrix{T}
    values::V
    ncol::Int
end

Base.size(field::TimeForcingMatrix) = (field.ncol, length(field.values))
@inline function Base.getindex(field::TimeForcingMatrix, column::Int, time::Int)
    @boundscheck checkbounds(field, column, time)
    return @inbounds field.values[time]
end
@inline function Base.setindex!(field::TimeForcingMatrix, value, column::Int, time::Int)
    @boundscheck checkbounds(field, column, time)
    @inbounds field.values[time] = value
    return value
end
Base.IndexStyle(::Type{<:TimeForcingMatrix}) = IndexCartesian()
Adapt.@adapt_structure TimeForcingMatrix

@inline _constant_forcing_matrix(value, dims::Tuple{Int, Int}) =
    ConstantForcingMatrix(value, dims)

@inline _map_forcing_field(f, field::ConstantForcingMatrix) =
    ConstantForcingMatrix(f(field.value), field.dims)
@inline _map_forcing_field(f, field::ColumnForcingMatrix) =
    ColumnForcingMatrix(f.(field.values), field.ntime)
@inline _map_forcing_field(f, field::TimeForcingMatrix) =
    TimeForcingMatrix(f.(field.values), field.ncol)
@inline _map_forcing_field(f, field) = f.(field)

@inline function _synthesized_time_values(dt_days::Vector{Float64})
    base = DateTime(2001, 1, 1, 12)
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
        return _constant_forcing_matrix(Float64(field), (ncol, ntime))
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return TimeForcingMatrix(Float64.(data), ncol)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Float64}(data)
    end
    error("`$name` must be a scalar, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

@inline function _forcing_bool_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Bool
        return _constant_forcing_matrix(field, (ncol, ntime))
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ntime || error("`$name` must have length $ntime.")
        return TimeForcingMatrix(Bool.(data), ncol)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Bool}(Bool.(data))
    end
    error("`$name` must be a Bool, a vector of length $ntime, or a matrix of size ($ncol, $ntime).")
end

@inline function _forcing_column_metadata_matrix(field, ncol::Int, ntime::Int, name::AbstractString)
    if field isa Number
        return _constant_forcing_matrix(Float64(field), (ncol, ntime))
    end
    data = collect(field)
    if ndims(data) == 1
        length(data) == ncol || error("`$name` vector input must have length $ncol.")
        return ColumnForcingMatrix(Float64.(data), ntime)
    elseif ndims(data) == 2
        size(data) == (ncol, ntime) || error("`$name` must have size ($ncol, $ntime).")
        return Matrix{Float64}(data)
    end
    error("`$name` must be a scalar, a vector of length $ncol, or a matrix of size ($ncol, $ntime).")
end

@inline function _normalize_air_pressure_temperature_mode(mode)
    mode_sym = mode isa Symbol ? mode : Symbol(lowercase(strip(String(mode))))
    mode_sym in (:annual_mean, :instantaneous) ||
        error("Unsupported air-pressure temperature mode '$mode'. Use :annual_mean or :instantaneous.")
    return mode_sym
end

@inline function _barometric_air_pressure(
    surface_height,
    air_temperature;
    sea_level_pressure::Real=DEFAULT_SEA_LEVEL_AIR_PRESSURE,
    gravity::Real=DEFAULT_GRAVITY,
    molar_mass_air::Real=DEFAULT_MOLAR_MASS_DRY_AIR,
    gas_constant::Real=DEFAULT_UNIVERSAL_GAS_CONSTANT,
)
    z = Float64(surface_height)
    T = Float64(air_temperature)
    if !isfinite(z) || !isfinite(T) || T <= 0.0
        return Float64(sea_level_pressure)
    end
    return Float64(sea_level_pressure) * exp(
        -Float64(gravity) * Float64(molar_mass_air) * z /
        (Float64(gas_constant) * T),
    )
end

function air_pressure_from_surface_height(
    surface_height,
    air_temperature;
    dt_days=nothing,
    time_values=nothing,
    temperature_mode=:annual_mean,
    sea_level_pressure::Real=DEFAULT_SEA_LEVEL_AIR_PRESSURE,
    gravity::Real=DEFAULT_GRAVITY,
    molar_mass_air::Real=DEFAULT_MOLAR_MASS_DRY_AIR,
    gas_constant::Real=DEFAULT_UNIVERSAL_GAS_CONSTANT,
)
    ncol, ntime = size(air_temperature)
    height = _forcing_column_metadata_matrix(surface_height, ncol, ntime, "surface_height")
    pressure = Matrix{Float64}(undef, ncol, ntime)
    mode = _normalize_air_pressure_temperature_mode(temperature_mode)

    if mode == :instantaneous
        @inbounds for t in 1:ntime, col in 1:ncol
            pressure[col, t] = _barometric_air_pressure(
                height[col, t],
                air_temperature[col, t];
                sea_level_pressure=sea_level_pressure,
                gravity=gravity,
                molar_mass_air=molar_mass_air,
                gas_constant=gas_constant,
            )
        end
        return pressure
    end

    isnothing(time_values) && error("`time_values` is required when temperature_mode=:annual_mean.")
    length(time_values) == ntime || error("`time_values` must have one entry per forcing timestep.")
    weights = isnothing(dt_days) ? ones(Float64, ntime) : Float64.(collect(dt_days))
    length(weights) == ntime || error("`dt_days` must have one entry per forcing timestep.")

    years = Dates.year.(DateTime.(collect(time_values)))
    start_idx = 1
    while start_idx <= ntime
        y = years[start_idx]
        stop_idx = start_idx
        while stop_idx < ntime && years[stop_idx + 1] == y
            stop_idx += 1
        end
        @inbounds for col in 1:ncol
            t_sum = 0.0
            valid_weight_sum = 0.0
            for t in start_idx:stop_idx
                T = air_temperature[col, t]
                w = weights[t]
                if isfinite(T) && isfinite(w) && w > 0.0
                    t_sum += T * w
                    valid_weight_sum += w
                end
            end
            Tmean = valid_weight_sum > 0.0 ? t_sum / valid_weight_sum : NaN
            p = _barometric_air_pressure(
                height[col, start_idx],
                Tmean;
                sea_level_pressure=sea_level_pressure,
                gravity=gravity,
                molar_mass_air=molar_mass_air,
                gas_constant=gas_constant,
            )
            for t in start_idx:stop_idx
                pressure[col, t] = p
            end
        end
        start_idx = stop_idx + 1
    end
    return pressure
end

@inline function _calendar_day_of_year(t::DateTime)
    seconds_today = hour(t) * 3600 + minute(t) * 60 + second(t) + millisecond(t) / 1000
    return Float64(dayofyear(t)) + seconds_today / 86_400.0
end

@inline function _solar_longitude_deg_from_calendar_day(day_of_year)
    days_since_j2000_like_year_start = day_of_year - 1.0
    mean_longitude = 280.46646 + 0.98564736 * days_since_j2000_like_year_start
    mean_anomaly = 357.52911 + 0.98560028 * days_since_j2000_like_year_start
    true_longitude = mean_longitude +
                     1.914602 * sind(mean_anomaly) +
                     0.019993 * sind(2.0 * mean_anomaly)
    return mod(true_longitude, 360.0)
end

struct ForcingCalendar{TV,DV,YV,SV}
    time_values::TV
    dt_days::DV
    day_of_year::YV
    solar_longitude_deg::SV
end

struct OptionalForcingField{V,M}
    values::V
    available::M
end

"""Time-varying atmospheric boundary conditions for a snowpack model."""
struct SnowpackForcing{C<:ForcingCalendar,F<:NamedTuple}
    calendar::C
    fields::F
end

struct PDDForcing{D,AT,SF,RF}
    dt_days::D
    air_temperature::AT
    snowfall_rate::SF
    rainfall_rate::RF
end

PDDForcing(forcing::SnowpackForcing) = PDDForcing(
    forcing.dt_days,
    forcing.air_temperature,
    forcing.snowfall_rate,
    forcing.rainfall_rate,
)

struct ITMForcing{D,AT,SF,RF,SW,QSW,HQSW,LAT,SH,HI,PDD}
    dt_days::D
    air_temperature::AT
    snowfall_rate::SF
    rainfall_rate::RF
    shortwave_down::SW
    q_sw_net::QSW
    has_q_sw_net::HQSW
    latitude_deg::LAT
    surface_height::SH
    ice_thickness::HI
    annual_pdd::PDD
end

ITMForcing(forcing::SnowpackForcing) = ITMForcing(
    forcing.dt_days,
    forcing.air_temperature,
    forcing.snowfall_rate,
    forcing.rainfall_rate,
    forcing.shortwave_down,
    forcing.q_sw_net,
    forcing.has_q_sw_net,
    forcing.latitude_deg,
    forcing.surface_height,
    forcing.ice_thickness,
    forcing.annual_pdd,
)

Adapt.@adapt_structure ForcingCalendar
Adapt.@adapt_structure OptionalForcingField
@adapt_structure SnowpackForcing
Adapt.@adapt_structure PDDForcing
Adapt.@adapt_structure ITMForcing

const _FORCING_CALENDAR_FIELD_NAMES = (:time_values, :dt_days, :day_of_year, :solar_longitude_deg)
const _FORCING_OPTIONAL_FIELD_NAMES = (
    :q_sw_net,
    :q_lw_down,
    :q_sh,
    :q_lh,
    :relative_humidity,
    :prescribed_albedo,
)
const _FORCING_AVAILABILITY_TO_FIELD = (
    has_q_sw_net=:q_sw_net,
    has_q_lw_down=:q_lw_down,
    has_q_sh=:q_sh,
    has_q_lh=:q_lh,
    has_relative_humidity=:relative_humidity,
    has_prescribed_albedo=:prescribed_albedo,
)

@generated function _forcing_field_property(fields::NamedTuple, ::Val{Name}) where {Name}
    availability_names = keys(_FORCING_AVAILABILITY_TO_FIELD)
    if Name in availability_names
        field_name = getfield(_FORCING_AVAILABILITY_TO_FIELD, Name)
        return :(getfield(getfield(fields, $(QuoteNode(field_name))), :available))
    elseif Name in _FORCING_OPTIONAL_FIELD_NAMES
        return :(getfield(getfield(fields, $(QuoteNode(Name))), :values))
    end
    return :(getfield(fields, $(QuoteNode(Name))))
end

@inline function Base.getproperty(forcing::SnowpackForcing, name::Symbol)
    name === :calendar && return getfield(forcing, :calendar)
    name === :fields && return getfield(forcing, :fields)
    if name in _FORCING_CALENDAR_FIELD_NAMES
        return getfield(getfield(forcing, :calendar), name)
    end
    return _forcing_field_property(getfield(forcing, :fields), Val(name))
end

Base.propertynames(::SnowpackForcing, private::Bool=false) = private ?
    (:calendar, :fields, _FORCING_CALENDAR_FIELD_NAMES..., _FORCING_MATRIX_FIELD_NAMES...) :
    (_FORCING_CALENDAR_FIELD_NAMES..., _FORCING_MATRIX_FIELD_NAMES...)

const _FORCING_MATRIX_FIELD_NAMES = (
    :air_temperature,
    :snowfall_rate,
    :rainfall_rate,
    :shortwave_down,
    :q_sw_net,
    :has_q_sw_net,
    :latitude_deg,
    :wind_speed,
    :q_lw_down,
    :has_q_lw_down,
    :q_sh,
    :has_q_sh,
    :q_lh,
    :has_q_lh,
    :relative_humidity,
    :has_relative_humidity,
    :surface_height,
    :ice_thickness,
    :annual_pdd,
    :air_pressure,
    :prescribed_albedo,
    :has_prescribed_albedo,
)

const _FORCING_COPY_FIELD_NAMES = (
    :dt_days,
    :day_of_year,
    :solar_longitude_deg,
    _FORCING_MATRIX_FIELD_NAMES...,
)

const _BESSI_DEVICE_FORCING_FIELD_NAMES = (
    :dt_days,
    :day_of_year,
    :solar_longitude_deg,
    _FORCING_MATRIX_FIELD_NAMES...,
)

@generated function get_fields(forcing::SnowpackForcing)
    entries = [
        :($(name) = getproperty(forcing, $(QuoteNode(name))))
        for name in _BESSI_DEVICE_FORCING_FIELD_NAMES
    ]
    return :((; $(entries...)))
end

@inline _device_forcing_fields(forcing::SnowpackForcing) = get_fields(forcing)
@inline _device_forcing_fields(forcing::NamedTuple) = forcing

struct SnowpackStepForcing{NF <: AbstractFloat}
    air_temperature::NF
    dt_days::NF
    snowfall_rate::NF
    rainfall_rate::NF
    shortwave_down::NF
    wind_speed::NF
    q_sw_net::NF
    q_lw_down::NF
    q_sh::NF
    q_lh::NF
    has_q_sw_net::Bool
    has_q_lw_down::Bool
    has_q_sh::Bool
    has_q_lh::Bool
    relative_humidity::NF
    has_relative_humidity::Bool
    air_pressure::NF
    prescribed_albedo::NF
    has_prescribed_albedo::Bool
    latitude_deg::NF
    day_of_year::NF
    solar_longitude_deg::NF
end

function SnowpackStepForcing(
    air_temperature,
    dt_days,
    snowfall_rate,
    rainfall_rate,
    shortwave_down,
    wind_speed;
    q_sw_net=zero(air_temperature),
    q_lw_down=zero(air_temperature),
    q_sh=zero(air_temperature),
    q_lh=zero(air_temperature),
    has_q_sw_net::Bool=false,
    has_q_lw_down::Bool=false,
    has_q_sh::Bool=false,
    has_q_lh::Bool=false,
    relative_humidity=zero(air_temperature),
    has_relative_humidity::Bool=false,
    air_pressure=oftype(air_temperature, 101_325.0),
    prescribed_albedo=zero(air_temperature),
    has_prescribed_albedo::Bool=false,
    latitude_deg=zero(air_temperature),
    day_of_year=zero(air_temperature),
    solar_longitude_deg=_solar_longitude_deg_from_calendar_day(day_of_year),
)
    return SnowpackStepForcing(
        air_temperature,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        q_sw_net,
        q_lw_down,
        q_sh,
        q_lh,
        has_q_sw_net,
        has_q_lw_down,
        has_q_sh,
        has_q_lh,
        relative_humidity,
        has_relative_humidity,
        air_pressure,
        prescribed_albedo,
        has_prescribed_albedo,
        latitude_deg,
        day_of_year,
        solar_longitude_deg,
    )
end

@adapt_structure SnowpackStepForcing

@inline _step_time_count(fields) = size(fields.air_temperature, 2)
@inline _step_dt(dt_days::Number, ::Int) = dt_days
@inline _step_dt(dt_days::AbstractVector, time_index::Int) = @inbounds dt_days[time_index]

Base.@propagate_inbounds function _step_forcing_at(forcing, idx::Int, time_index::Int)
    @inbounds begin
        air_temperature = forcing.air_temperature[idx, time_index]
        return SnowpackStepForcing(
            air_temperature,
            _step_dt(forcing.dt_days, time_index),
            forcing.snowfall_rate[idx, time_index],
            forcing.rainfall_rate[idx, time_index],
            forcing.shortwave_down[idx, time_index],
            forcing.wind_speed[idx, time_index];
            q_sw_net=forcing.q_sw_net[idx, time_index],
            has_q_sw_net=forcing.has_q_sw_net[idx, time_index],
            q_lw_down=forcing.q_lw_down[idx, time_index],
            has_q_lw_down=forcing.has_q_lw_down[idx, time_index],
            q_sh=forcing.q_sh[idx, time_index],
            has_q_sh=forcing.has_q_sh[idx, time_index],
            q_lh=forcing.q_lh[idx, time_index],
            has_q_lh=forcing.has_q_lh[idx, time_index],
            relative_humidity=forcing.relative_humidity[idx, time_index],
            has_relative_humidity=forcing.has_relative_humidity[idx, time_index],
            air_pressure=forcing.air_pressure[idx, time_index],
            prescribed_albedo=forcing.prescribed_albedo[idx, time_index],
            has_prescribed_albedo=forcing.has_prescribed_albedo[idx, time_index],
            latitude_deg=forcing.latitude_deg[idx, time_index],
            day_of_year=forcing.day_of_year[time_index],
            solar_longitude_deg=forcing.solar_longitude_deg[time_index],
        )
    end
end

function _optional_forcing_field(value, has_value, default, dims, column_count, ntime, name)
    values = isnothing(value) ?
             _constant_forcing_matrix(default, dims) :
             _forcing_numeric_matrix(value, column_count, ntime, name)
    available = isnothing(value) ?
                _constant_forcing_matrix(false, dims) :
                isnothing(has_value) ?
                _constant_forcing_matrix(true, dims) :
                _forcing_bool_matrix(has_value, column_count, ntime, "has_$(name)")
    return values, available
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
    q_sw_net=nothing,
    has_q_sw_net=nothing,
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
    surface_height=nothing,
    ice_thickness=nothing,
    annual_pdd=nothing,
    air_pressure_temperature_mode=:instantaneous,
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
                relative_humidity, air_pressure, surface_height, ice_thickness, annual_pdd, prescribed_albedo)
                isnothing(field) || ((column_count = _forcing_column_count(field, ntime)); break)
            end
        end
    end
    column_count > 0 || error("`ncol` must be positive.")

    if has_user
        isnothing(air_temperature_c) && error("`air_temperature_c` is required.")
        isnothing(snowfall_mm_day) && error("`snowfall_mm_day` is required.")
        isnothing(rainfall_mm_day) && error("`rainfall_mm_day` is required.")
        air_temperature = _map_forcing_field(
            temperature -> temperature + 273.15,
            _forcing_numeric_matrix(air_temperature_c, column_count, ntime, "air_temperature_c"),
        )
        snowfall_rate = _map_forcing_field(
            rate -> rate / 86_400.0,
            _forcing_numeric_matrix(snowfall_mm_day, column_count, ntime, "snowfall_mm_day"),
        )
        rainfall_rate = _map_forcing_field(
            rate -> rate / 86_400.0,
            _forcing_numeric_matrix(rainfall_mm_day, column_count, ntime, "rainfall_mm_day"),
        )
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
    q_sw_net_m, has_q_sw_net_m = _optional_forcing_field(
        q_sw_net, has_q_sw_net, 0.0, dims, column_count, ntime, "q_sw_net",
    )
    wind_speed_m = isnothing(wind_speed) ? _constant_forcing_matrix(5.0, dims) : _forcing_numeric_matrix(wind_speed, column_count, ntime, "wind_speed")
    latitude_deg_m = isnothing(latitude_deg) ? fill(NaN, dims) : _forcing_column_metadata_matrix(latitude_deg, column_count, ntime, "latitude_deg")
    q_lw_down_m, has_q_lw_down_m = _optional_forcing_field(
        q_lw_down, has_q_lw_down, 0.0, dims, column_count, ntime, "q_lw_down",
    )
    q_sh_m, has_q_sh_m = _optional_forcing_field(
        q_sh, has_q_sh, 0.0, dims, column_count, ntime, "q_sh",
    )
    q_lh_m, has_q_lh_m = _optional_forcing_field(
        q_lh, has_q_lh, 0.0, dims, column_count, ntime, "q_lh",
    )
    relative_humidity_m = if !isnothing(relative_humidity)
        _forcing_numeric_matrix(relative_humidity, column_count, ntime, "relative_humidity")
    else
        _constant_forcing_matrix(0.0, dims)
    end
    has_relative_humidity_m = if isnothing(relative_humidity)
        _constant_forcing_matrix(false, dims)
    elseif isnothing(has_relative_humidity)
        _map_forcing_field(isfinite, relative_humidity_m)
    else
        _forcing_bool_matrix(has_relative_humidity, column_count, ntime, "has_relative_humidity")
    end
    if !isnothing(relative_humidity)
        if isnothing(has_relative_humidity)
            relative_humidity_m = _map_forcing_field(
                value -> isfinite(value) ? value : 0.0,
                relative_humidity_m,
            )
        else
            relative_humidity_m = Matrix(relative_humidity_m)
            relative_humidity_m[.!has_relative_humidity_m] .= 0.0
        end
    end
    time_values_v = isnothing(time_values) ? _synthesized_time_values(dt_days_v) : DateTime.(collect(time_values))
    length(time_values_v) == dims[2] || error("`time_values` must have one entry per forcing timestep.")
    surface_height_m = isnothing(surface_height) ? fill(NaN, dims) : _forcing_column_metadata_matrix(surface_height, column_count, ntime, "surface_height")
    ice_thickness_m = isnothing(ice_thickness) ? fill(NaN, dims) : _forcing_column_metadata_matrix(ice_thickness, column_count, ntime, "ice_thickness")
    annual_pdd_m = isnothing(annual_pdd) ? fill(NaN, dims) : _forcing_column_metadata_matrix(annual_pdd, column_count, ntime, "annual_pdd")
    air_pressure_m = if !isnothing(air_pressure)
        _forcing_numeric_matrix(air_pressure, column_count, ntime, "air_pressure")
    elseif !isnothing(surface_height)
        air_pressure_from_surface_height(
            surface_height_m,
            air_temperature;
            dt_days=dt_days_v,
            time_values=time_values_v,
            temperature_mode=air_pressure_temperature_mode,
        )
    else
        fill(DEFAULT_SEA_LEVEL_AIR_PRESSURE, dims)
    end
    prescribed_albedo_m, has_prescribed_albedo_m = _optional_forcing_field(
        prescribed_albedo, has_prescribed_albedo, 0.0, dims, column_count, ntime, "prescribed_albedo",
    )

    for (name, field) in (
        ("snowfall_rate", snowfall_rate),
        ("rainfall_rate", rainfall_rate),
        ("shortwave_down", shortwave_down_m),
        ("q_sw_net", q_sw_net_m),
        ("has_q_sw_net", has_q_sw_net_m),
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
        ("surface_height", surface_height_m),
        ("ice_thickness", ice_thickness_m),
        ("annual_pdd", annual_pdd_m),
        ("air_pressure", air_pressure_m),
        ("prescribed_albedo", prescribed_albedo_m),
        ("has_prescribed_albedo", has_prescribed_albedo_m),
    )
        _ensure_matching_field_sizes(dims, name, field)
    end

    day_of_year_v = _calendar_day_of_year.(time_values_v)
    solar_longitude_deg_v = _solar_longitude_deg_from_calendar_day.(day_of_year_v)

    calendar = ForcingCalendar(
        time_values_v,
        dt_days_v,
        day_of_year_v,
        solar_longitude_deg_v,
    )
    fields = (
        air_temperature=air_temperature,
        snowfall_rate=snowfall_rate,
        rainfall_rate=rainfall_rate,
        shortwave_down=shortwave_down_m,
        q_sw_net=OptionalForcingField(q_sw_net_m, has_q_sw_net_m),
        latitude_deg=latitude_deg_m,
        wind_speed=wind_speed_m,
        q_lw_down=OptionalForcingField(q_lw_down_m, has_q_lw_down_m),
        q_sh=OptionalForcingField(q_sh_m, has_q_sh_m),
        q_lh=OptionalForcingField(q_lh_m, has_q_lh_m),
        relative_humidity=OptionalForcingField(relative_humidity_m, has_relative_humidity_m),
        surface_height=surface_height_m,
        ice_thickness=ice_thickness_m,
        annual_pdd=annual_pdd_m,
        air_pressure=air_pressure_m,
        prescribed_albedo=OptionalForcingField(prescribed_albedo_m, has_prescribed_albedo_m),
    )
    return SnowpackForcing(calendar, fields)
end

function update_air_pressure!(
    forcing::SnowpackForcing;
    temperature_mode=:instantaneous,
    sea_level_pressure::Real=DEFAULT_SEA_LEVEL_AIR_PRESSURE,
    gravity::Real=DEFAULT_GRAVITY,
    molar_mass_air::Real=DEFAULT_MOLAR_MASS_DRY_AIR,
    gas_constant::Real=DEFAULT_UNIVERSAL_GAS_CONSTANT,
)
    forcing.air_pressure .= air_pressure_from_surface_height(
        forcing.surface_height,
        forcing.air_temperature;
        dt_days=forcing.dt_days,
        time_values=forcing.time_values,
        temperature_mode=temperature_mode,
        sea_level_pressure=sea_level_pressure,
        gravity=gravity,
        molar_mass_air=molar_mass_air,
        gas_constant=gas_constant,
    )
    return forcing
end
