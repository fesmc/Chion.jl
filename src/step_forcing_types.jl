"""
Step-forcing container types and field-to-forcing conversion helpers.
"""

struct SnowpackStepForcing{NF <: AbstractFloat}
    air_temperature::NF
    precipitation_rate::NF
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
    diurnal_shortwave::Bool
    latitude::NF
    day_of_year::NF
end

struct SnowpackStepFields{
        DT,
        AT,
        ST,
        RT,
        SWT,
        WST,
        QLWT,
        HQLWT,
        QSHT,
        HQSHT,
        QLHT,
        HQLHT,
    }
    dt_days::DT
    air_temperature::AT
    snowfall_rate::ST
    rainfall_rate::RT
    shortwave_down::SWT
    wind_speed::WST
    q_lw_down::QLWT
    has_q_lw_down::HQLWT
    q_sh::QSHT
    has_q_sh::HQSHT
    q_lh::QLHT
    has_q_lh::HQLHT
end

@inline function _assert_step_field_shape(name::AbstractString, field, field_shape)
    size(field) == field_shape || error("`$name` must match `air_temperature`.")
    return nothing
end

function SnowpackStepFields(;
    dt_days,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    shortwave_down,
    wind_speed,
    q_lw_down,
    has_q_lw_down,
    q_sh,
    has_q_sh,
    q_lh,
    has_q_lh,
)
    field_shape = size(air_temperature)
    _assert_step_field_shape("snowfall_rate", snowfall_rate, field_shape)
    _assert_step_field_shape("rainfall_rate", rainfall_rate, field_shape)
    _assert_step_field_shape("shortwave_down", shortwave_down, field_shape)
    _assert_step_field_shape("wind_speed", wind_speed, field_shape)
    _assert_step_field_shape("q_lw_down", q_lw_down, field_shape)
    _assert_step_field_shape("has_q_lw_down", has_q_lw_down, field_shape)
    _assert_step_field_shape("q_sh", q_sh, field_shape)
    _assert_step_field_shape("has_q_sh", has_q_sh, field_shape)
    _assert_step_field_shape("q_lh", q_lh, field_shape)
    _assert_step_field_shape("has_q_lh", has_q_lh, field_shape)

    ntime = size(air_temperature, 2)
    if dt_days isa AbstractVector
        length(dt_days) == ntime || error("`dt_days` must have one entry per forcing time step.")
    else
        ntime == 1 || error("Scalar `dt_days` requires a single forcing time step.")
    end

    return SnowpackStepFields(
        dt_days,
        air_temperature,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        q_lw_down,
        has_q_lw_down,
        q_sh,
        has_q_sh,
        q_lh,
        has_q_lh,
    )
end

function SnowpackStepFields(forcing)
    return SnowpackStepFields(
        dt_days=forcing.dt_days,
        air_temperature=forcing.air_temperature,
        snowfall_rate=forcing.snowfall_rate,
        rainfall_rate=forcing.rainfall_rate,
        shortwave_down=forcing.shortwave_down,
        wind_speed=forcing.wind_speed,
        q_lw_down=forcing.q_lw_down,
        has_q_lw_down=forcing.has_q_lw_down,
        q_sh=forcing.q_sh,
        has_q_sh=forcing.has_q_sh,
        q_lh=forcing.q_lh,
        has_q_lh=forcing.has_q_lh,
    )
end

@inline _step_time_count(fields::SnowpackStepFields) = size(fields.air_temperature, 2)
@inline _step_dt(dt_days::Number, ::Int) = dt_days
@inline _step_dt(dt_days::AbstractVector, time_index::Int) = @inbounds dt_days[time_index]
@inline _step_dt(dt_days::CUDA.CuArray, time_index::Int) = CUDA.@allowscalar dt_days[time_index]

@inline function _step_forcing_from_fields(
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    dt_days,
    shortwave_down,
    wind_speed,
    q_lw_down,
    has_q_lw_down::Bool,
    q_sh,
    has_q_sh::Bool,
    q_lh,
    has_q_lh::Bool,
)
    return SnowpackStepForcing(
        air_temperature,
        snowfall_rate + rainfall_rate,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        zero(air_temperature),
        q_lw_down,
        q_sh,
        q_lh,
        false,
        has_q_lw_down,
        has_q_sh,
        has_q_lh,
        false,
        zero(air_temperature),
        zero(air_temperature),
    )
end

@inline function _step_forcing_from_fields(
    fields::SnowpackStepFields,
    idx::Int,
    time_index::Int,
)
    return _step_forcing_from_fields(
        @inbounds(fields.air_temperature[idx, time_index]),
        @inbounds(fields.snowfall_rate[idx, time_index]),
        @inbounds(fields.rainfall_rate[idx, time_index]),
        _step_dt(fields.dt_days, time_index),
        @inbounds(fields.shortwave_down[idx, time_index]),
        @inbounds(fields.wind_speed[idx, time_index]),
        @inbounds(fields.q_lw_down[idx, time_index]),
        @inbounds(fields.has_q_lw_down[idx, time_index]),
        @inbounds(fields.q_sh[idx, time_index]),
        @inbounds(fields.has_q_sh[idx, time_index]),
        @inbounds(fields.q_lh[idx, time_index]),
        @inbounds(fields.has_q_lh[idx, time_index]),
    )
end

function SnowpackStepForcing(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate,
    dt_days;
    snowfall_rate=zero(air_temperature),
    rainfall_rate=zero(air_temperature),
    shortwave_down=oftype(air_temperature, 400.0),
    wind_speed=oftype(air_temperature, 10.0),
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
    diurnal_shortwave::Bool=false,
    latitude=zero(air_temperature),
    day_of_year=zero(air_temperature),
)
    NF = number_type(c)
    return SnowpackStepForcing(
        convert(NF, air_temperature),
        convert(NF, precipitation_rate),
        convert(NF, dt_days),
        convert(NF, snowfall_rate),
        convert(NF, rainfall_rate),
        convert(NF, shortwave_down),
        convert(NF, wind_speed),
        convert(NF, isnothing(q_sw_net) ? zero(NF) : q_sw_net),
        convert(NF, isnothing(q_lw_down) ? zero(NF) : q_lw_down),
        convert(NF, isnothing(q_sh) ? zero(NF) : q_sh),
        convert(NF, isnothing(q_lh) ? zero(NF) : q_lh),
        !isnothing(q_sw_net),
        !isnothing(q_lw_down),
        !isnothing(q_sh),
        !isnothing(q_lh),
        diurnal_shortwave,
        convert(NF, latitude),
        convert(NF, day_of_year),
    )
end

@adapt_structure SnowpackStepForcing
@adapt_structure SnowpackStepFields

@inline _default_snow_fraction(c::SnowpackPhysicalConstants, air_temperature) =
    air_temperature > c.T0 ? zero(air_temperature) : one(air_temperature)

@inline _diagnosed_shortwave_down(shortwave_down) =
    isnothing(shortwave_down) ? 400.0 : max(shortwave_down, zero(shortwave_down))

function _resolve_step_partition(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate;
    snow_fraction=nothing,
    f_s=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    shortwave_down=nothing,
    p_snow=nothing,
    p_rain=nothing,
    s_boa=nothing,
)
    resolved_snow_fraction = _resolve_keyword_alias(snow_fraction, f_s, "snow_fraction", "f_s")
    resolved_snowfall_rate = _resolve_keyword_alias(snowfall_rate, p_snow, "snowfall_rate", "p_snow")
    resolved_rainfall_rate = _resolve_keyword_alias(rainfall_rate, p_rain, "rainfall_rate", "p_rain")
    resolved_shortwave_down = _resolve_keyword_alias(shortwave_down, s_boa, "shortwave_down", "s_boa")

    if !isnothing(resolved_snowfall_rate) || !isnothing(resolved_rainfall_rate)
        return (
            snowfall_rate=isnothing(resolved_snowfall_rate) ? zero(precipitation_rate) : resolved_snowfall_rate,
            rainfall_rate=isnothing(resolved_rainfall_rate) ? zero(precipitation_rate) : resolved_rainfall_rate,
            shortwave_down=resolved_shortwave_down,
        )
    end

    snowfall_fraction = isnothing(resolved_snow_fraction) ?
        _default_snow_fraction(c, air_temperature) :
        resolved_snow_fraction
    rainfall = precipitation_rate * (one(precipitation_rate) - snowfall_fraction)
    return (
        snowfall_rate=precipitation_rate - rainfall,
        rainfall_rate=rainfall,
        shortwave_down=resolved_shortwave_down,
    )
end

function _resolved_step_forcing(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate,
    dt_days;
    snow_fraction=nothing,
    f_s=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    shortwave_down=nothing,
    p_snow=nothing,
    p_rain=nothing,
    s_boa=nothing,
    wind_speed=oftype(air_temperature, 10.0),
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
    diurnal_shortwave::Bool=false,
    latitude=nothing,
    day_of_year=nothing,
)
    if diurnal_shortwave && (isnothing(latitude) || isnothing(day_of_year))
        error("`diurnal_shortwave=true` requires both `latitude` and `day_of_year`.")
    end

    resolved = _resolve_step_partition(
        c,
        air_temperature,
        precipitation_rate;
        snow_fraction=snow_fraction,
        f_s=f_s,
        snowfall_rate=snowfall_rate,
        rainfall_rate=rainfall_rate,
        shortwave_down=shortwave_down,
        p_snow=p_snow,
        p_rain=p_rain,
        s_boa=s_boa,
    )
    return SnowpackStepForcing(
        c,
        air_temperature,
        precipitation_rate,
        dt_days;
        snowfall_rate=resolved.snowfall_rate,
        rainfall_rate=resolved.rainfall_rate,
        shortwave_down=_diagnosed_shortwave_down(resolved.shortwave_down),
        wind_speed=wind_speed,
        q_sw_net=q_sw_net,
        q_lw_down=q_lw_down,
        q_sh=q_sh,
        q_lh=q_lh,
        diurnal_shortwave=diurnal_shortwave,
        latitude=isnothing(latitude) ? zero(air_temperature) : latitude,
        day_of_year=isnothing(day_of_year) ? zero(air_temperature) : day_of_year,
    )
end
