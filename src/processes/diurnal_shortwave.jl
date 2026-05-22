"""
Energy-conserving adaptive diurnal shortwave substep helpers.
"""

@inline _diurnal_obliquity_deg(x) = oftype(x, 23.439291)
@inline _solar_declination_deg(solar_longitude_deg) =
    asind(sind(_diurnal_obliquity_deg(solar_longitude_deg)) * sind(solar_longitude_deg))

@inline function _sunset_hour_angle(latitude_deg, declination_deg)
    cos_h0 = -tand(latitude_deg) * tand(declination_deg)
    if cos_h0 >= one(cos_h0)
        return zero(cos_h0)
    elseif cos_h0 <= -one(cos_h0)
        return oftype(cos_h0, π)
    end
    return acos(cos_h0)
end

@inline function _diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg)
    declination_deg = _solar_declination_deg(solar_longitude_deg)
    h0 = _sunset_hour_angle(latitude_deg, declination_deg)
    lat_rad = deg2rad(latitude_deg)
    dec_rad = deg2rad(declination_deg)
    sin_lat_sin_dec = sin(lat_rad) * sin(dec_rad)
    cos_lat_cos_dec = cos(lat_rad) * cos(dec_rad)
    daylight_integral = oftype(latitude_deg, 2) *
                        (h0 * sin_lat_sin_dec + cos_lat_cos_dec * sin(h0))
    return (
        declination_deg=declination_deg,
        sunset_hour_angle=h0,
        daylight_integral=daylight_integral,
        sin_lat_sin_dec=sin_lat_sin_dec,
        cos_lat_cos_dec=cos_lat_cos_dec,
    )
end

@inline function _diurnal_shortwave_interval_average(
    shortwave_daily_mean,
    latitude_deg,
    solar_longitude_deg,
    hour_angle_start,
    hour_angle_end,
)
    interval_width = hour_angle_end - hour_angle_start
    if shortwave_daily_mean <= zero(shortwave_daily_mean) ||
       interval_width <= zero(interval_width) ||
       !isfinite(latitude_deg) ||
       !isfinite(solar_longitude_deg)
        return zero(shortwave_daily_mean)
    end

    terms = _diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg)
    if terms.daylight_integral <= eps(typeof(float(terms.daylight_integral))) ||
       terms.sunset_hour_angle <= zero(terms.sunset_hour_angle)
        return zero(shortwave_daily_mean)
    end

    daylight_start = max(hour_angle_start, -terms.sunset_hour_angle)
    daylight_end = min(hour_angle_end, terms.sunset_hour_angle)
    daylight_end <= daylight_start && return zero(shortwave_daily_mean)

    daylight_integral = (daylight_end - daylight_start) * terms.sin_lat_sin_dec +
                        terms.cos_lat_cos_dec * (sin(daylight_end) - sin(daylight_start))
    scale = shortwave_daily_mean * oftype(shortwave_daily_mean, 2π) / terms.daylight_integral
    return max(scale * daylight_integral / interval_width, zero(shortwave_daily_mean))
end

@inline function _diurnal_shortwave_peak_flux(shortwave_daily_mean, latitude_deg, solar_longitude_deg)
    if shortwave_daily_mean <= zero(shortwave_daily_mean) ||
       !isfinite(latitude_deg) ||
       !isfinite(solar_longitude_deg)
        return zero(shortwave_daily_mean)
    end
    terms = _diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg)
    if terms.daylight_integral <= eps(typeof(float(terms.daylight_integral))) ||
       terms.sunset_hour_angle <= zero(terms.sunset_hour_angle)
        return zero(shortwave_daily_mean)
    end
    scale = shortwave_daily_mean * oftype(shortwave_daily_mean, 2π) / terms.daylight_integral
    return max(
        scale * (terms.sin_lat_sin_dec + terms.cos_lat_cos_dec),
        zero(shortwave_daily_mean),
    )
end

@inline function _diurnal_temperature_interval_average(
    air_temperature_daily_mean,
    amplitude,
    hour_angle_start,
    hour_angle_end,
)
    interval_width = hour_angle_end - hour_angle_start
    if amplitude <= zero(amplitude) || interval_width <= zero(interval_width)
        return air_temperature_daily_mean
    end
    return air_temperature_daily_mean +
           amplitude * (sin(hour_angle_end) - sin(hour_angle_start)) / interval_width
end

@inline function _diurnal_shortwave_substep_count(
    dt_days,
    shortwave_daily_mean,
    air_temperature,
    min_air_temperature,
    latitude_deg,
    solar_longitude_deg,
    threshold,
    max_substeps::Int,
)
    max_substeps <= 1 && return 1
    if dt_days < oftype(dt_days, 0.75) ||
       dt_days > oftype(dt_days, 1.25) ||
       shortwave_daily_mean <= zero(shortwave_daily_mean) ||
       air_temperature <= min_air_temperature ||
       !isfinite(air_temperature) ||
       !isfinite(latitude_deg) ||
       !isfinite(solar_longitude_deg)
        return 1
    end

    peak_flux = _diurnal_shortwave_peak_flux(shortwave_daily_mean, latitude_deg, solar_longitude_deg)
    max(peak_flux - shortwave_daily_mean, zero(shortwave_daily_mean)) <= threshold && return 1

    terms = _diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg)
    terms.sunset_hour_angle <= zero(terms.sunset_hour_angle) && return 1
    return max_substeps
end

@inline function _diurnal_substep_forcing(
    forcing::SnowpackStepForcing,
    fraction,
    air_temperature,
    shortwave_down,
    q_sw_net,
)
    return SnowpackStepForcing(
        air_temperature,
        forcing.precipitation_rate,
        forcing.dt_days * fraction,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        shortwave_down,
        forcing.wind_speed,
        q_sw_net,
        forcing.q_lw_down,
        forcing.q_sh,
        forcing.q_lh,
        forcing.has_q_sw_net,
        forcing.has_q_lw_down,
        forcing.has_q_sh,
        forcing.has_q_lh,
        forcing.relative_humidity,
        forcing.has_relative_humidity,
        forcing.air_pressure,
        forcing.prescribed_albedo,
        forcing.has_prescribed_albedo,
        false,
        forcing.latitude_deg,
        forcing.day_of_year,
        forcing.solar_longitude_deg,
        forcing.diurnal_shortwave_threshold,
        forcing.diurnal_shortwave_max_substeps,
        forcing.diurnal_shortwave_min_air_temperature,
        false,
        forcing.diurnal_temperature_amplitude,
    )
end
