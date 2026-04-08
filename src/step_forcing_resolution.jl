"""
Keyword resolution and normalization for scalar stepping inputs.
"""

@inline _default_snow_fraction(c::SnowpackPhysicalConstants, air_temperature) =
    air_temperature > c.T0 ? zero(air_temperature) : one(air_temperature)

function _resolve_step_forcing(
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
    snowfall = precipitation_rate - rainfall
    return (
        snowfall_rate=snowfall,
        rainfall_rate=rainfall,
        shortwave_down=resolved_shortwave_down,
    )
end

@inline _diagnosed_shortwave_down(shortwave_down) =
    isnothing(shortwave_down) ? 400.0 : max(shortwave_down, zero(shortwave_down))

function _validate_diurnal_configuration(diurnal_shortwave::Bool, latitude, day_of_year)
    if diurnal_shortwave && (isnothing(latitude) || isnothing(day_of_year))
        error("`diurnal_shortwave=true` requires both `latitude` and `day_of_year`.")
    end
    return nothing
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
    _validate_diurnal_configuration(diurnal_shortwave, latitude, day_of_year)
    resolved = _resolve_step_forcing(
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
