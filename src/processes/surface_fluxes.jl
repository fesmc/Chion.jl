"""
Bare-ice and diurnal surface-flux helpers used by `step.jl`.
"""

# Daily-mean top-of-atmosphere shortwave below which the shortwave cloud proxy
# is undefined (polar night) and the night cloud fraction is used instead.
const LONGWAVE_CLOUD_PROXY_MIN_TOA = 50.0
const SOLAR_CONSTANT = 1361.0

"""Daily-mean top-of-atmosphere shortwave (W m⁻²) from latitude and season."""
@inline function _daily_toa_shortwave(latitude_deg, solar_longitude_deg, day_of_year)
    terms = _diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg)
    eccentricity = one(latitude_deg) + oftype(latitude_deg, 0.033) * cos(oftype(latitude_deg, 2π) * day_of_year / oftype(latitude_deg, 365))
    return oftype(latitude_deg, SOLAR_CONSTANT) * eccentricity * max(terms.daylight_integral, zero(latitude_deg)) /
           oftype(latitude_deg, 2π)
end

"""
    _cloud_proxy_emissivity(c, forcing)

Effective atmospheric emissivity relative to the near-surface air temperature,
ε = ε₀ + ε_T (T_a − T₀) + ε_n n. The cloudiness n = 1 − τ/τ_clear(z) uses the
daily shortwave transmissivity τ = SW↓/SW_TOA, which the forcing already
provides; in the polar night (or for sub-daily forcing, where SW↓ is not a
daily mean) n falls back to a constant night cloudiness.
"""
@inline function _cloud_proxy_emissivity(c::SnowpackPhysicalConstants, forcing::SnowpackStepForcing)
    T = forcing.air_temperature
    toa = _daily_toa_shortwave(forcing.latitude_deg, forcing.solar_longitude_deg, forcing.day_of_year)
    # Forcing without surface height (NaN) uses the sea-level clear-sky transmissivity.
    height = ifelse(isfinite(forcing.surface_height), max(forcing.surface_height, zero(T)), zero(T))
    clear_transmissivity = c.lw_clear_sky_transmissivity +
                           c.lw_clear_sky_transmissivity_per_km * height / oftype(T, 1000)
    daily = (toa > oftype(T, LONGWAVE_CLOUD_PROXY_MIN_TOA)) & (forcing.dt_days >= oftype(T, 0.75)) &
            isfinite(forcing.latitude_deg)
    cloudiness = ifelse(
        daily,
        clamp(one(T) - forcing.shortwave_down / _safe_positive(toa * clear_transmissivity), zero(T), one(T)),
        c.lw_night_cloud_fraction,
    )
    emissivity = c.lw_emissivity_base + c.lw_emissivity_temperature_slope * (T - c.T0) +
                 c.lw_emissivity_cloud_slope * cloudiness
    return clamp(emissivity, oftype(T, 0.4), oftype(T, 1.3))
end

"""
    _with_parameterized_longwave(c, forcing)

For the `:cloud_proxy` longwave scheme, resolve the daily longwave down once from
the daily forcing and pass it on as a fixed flux, so the diurnal substeps reuse
the daily cloud proxy. Prescribed longwave and the graybody scheme pass through.
"""
@inline function _with_parameterized_longwave(c::SnowpackPhysicalConstants, forcing::SnowpackStepForcing)
    (c.longwave_scheme == LONGWAVE_CLOUD_PROXY) & !forcing.has_q_lw_down || return forcing
    q_lw_down = _cloud_proxy_emissivity(c, forcing) * c.σ * forcing.air_temperature^4
    return _with_step_forcing(forcing, (; q_lw_down, has_q_lw_down=true))
end

@inline function _surface_vapor_fluxes(vapor_mass, latent_heat_flux)
    return (
        vapor_mass=vapor_mass,
        sublimation_mass=max(-vapor_mass, zero(vapor_mass)),
        latent_heat_flux=latent_heat_flux,
    )
end

"""
    _resolved_nonshortwave_surface_flux_components(c, air_temperature, rainfall_rate, dt_seconds, surface_temperature, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Resolve longwave, sensible, latent, and rain heat flux components for a known
surface temperature. Returned fluxes are instantaneous energy fluxes in model
surface-flux units.
"""
@inline function _resolved_nonshortwave_surface_flux_components(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    surface_temperature,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
    use_relative_humidity::Bool,
    relative_humidity,
    air_pressure,
    wind_speed,
    z0m,
    emissivity,
)
    longwave_down = use_q_lw_down ? q_lw_down_value : c.σ * c.ϵ_air * air_temperature^4
    longwave_flux = _uses_semix_seb(c) ?
        emissivity * (longwave_down - c.σ * surface_temperature^4) :
        longwave_down - c.σ * c.ϵ_snow * surface_temperature^4
    semix_sensible_constant, semix_sensible_linear, semix_latent_constant, semix_latent_linear =
        _semix_turbulent_flux_linearized(
            surface_temperature, c, air_temperature, relative_humidity,
            air_pressure, wind_speed, z0m,
        )
    sensible_heat_flux = use_q_sh ? q_sh_value :
                         _uses_semix_turbulence(c) ? semix_sensible_constant - semix_sensible_linear * surface_temperature :
                         c.D_sh * (air_temperature - surface_temperature)
    latent_heat_flux = use_q_lh ? q_lh_value :
                       use_relative_humidity ?
                       (_uses_semix_turbulence(c) ? semix_latent_constant - semix_latent_linear * surface_temperature :
                        _bessi_latent_vapor_flux(surface_temperature, c, air_temperature, relative_humidity, air_pressure)) :
                       zero(dt_seconds)
    rain_heat_flux = rainfall_rate * c.cw * (air_temperature - c.T0)
    return longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux
end

"""
    _resolved_bare_ice_surface_flux_components(c, air_temperature, rainfall_rate, dt_seconds, shortwave_down, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Resolve all bare-ice surface-flux components, including absorbed shortwave
energy.
"""
@inline function _resolved_bare_ice_surface_flux_components(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    shortwave_down,
    use_q_sw_net::Bool,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
    use_relative_humidity::Bool,
    relative_humidity,
    air_pressure,
    wind_speed,
    surface_albedo,
)
    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        max(shortwave_down, zero(dt_seconds)) *
        (one(dt_seconds) - clamp(surface_albedo, zero(surface_albedo), one(surface_albedo)))
    longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux =
        _resolved_nonshortwave_surface_flux_components(
            c,
            air_temperature,
            rainfall_rate,
            dt_seconds,
            c.T0,
            use_q_lw_down,
            q_lw_down_value,
            use_q_sh,
            q_sh_value,
            use_q_lh,
            q_lh_value,
            use_relative_humidity,
            relative_humidity,
            air_pressure,
            wind_speed,
            c.semix_z0m_ice,
            c.eps_ice,
        )
    # Bare ice remains a solid surface at the melting point: its direct
    # vapour exchange is sublimation/deposition and therefore carries Lᵥ+Lₘ.
    # (Snow with liquid water uses Lᵥ in the snow-surface routine.)
    latent_heat_flux = !use_q_lh && use_relative_humidity && !_uses_semix_turbulence(c) ?
        _bessi_latent_vapor_flux(
            c.T0, c, air_temperature, relative_humidity, air_pressure, c.Lv + c.Lm,
        ) : latent_heat_flux
    return (
        absorbed_shortwave,
        longwave_flux,
        sensible_heat_flux,
        latent_heat_flux,
        rain_heat_flux,
    )
end

"""
    _bare_ice_surface_mass_fluxes_resolved(c, air_temperature, rainfall_rate, dt_seconds, shortwave_down, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value, use_relative_humidity, relative_humidity, air_pressure, surface_albedo)

Return bare-ice melt mass, vapor-mass correction, and net SMB mass change.
The vapor term follows the BESSI sign convention: positive values add mass by
deposition, negative values remove mass by sublimation.
"""
function _bare_ice_surface_mass_fluxes_resolved(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    shortwave_down,
    use_q_sw_net::Bool,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
    use_relative_humidity::Bool,
    relative_humidity,
    air_pressure,
    wind_speed,
    surface_albedo,
)
    absorbed_shortwave, longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux =
        _resolved_bare_ice_surface_flux_components(
            c,
            air_temperature,
            rainfall_rate,
            dt_seconds,
            shortwave_down,
            use_q_sw_net,
            q_sw_net_value,
            use_q_lw_down,
            q_lw_down_value,
            use_q_sh,
            q_sh_value,
            use_q_lh,
            q_lh_value,
            use_relative_humidity,
            relative_humidity,
            air_pressure,
            wind_speed,
            surface_albedo,
        )
    net_surface_flux = absorbed_shortwave + longwave_flux + sensible_heat_flux + latent_heat_flux + rain_heat_flux
    melt_mass = max(net_surface_flux, zero(net_surface_flux)) * dt_seconds / c.Lm
    vapor_mass = latent_heat_flux * dt_seconds / (c.Lv + c.Lm)
    return (
        melt_mass=melt_mass,
        _surface_vapor_fluxes(vapor_mass, latent_heat_flux)...,
        net_mass_change=vapor_mass - melt_mass,
        absorbed_shortwave=absorbed_shortwave,
        longwave_flux=longwave_flux,
        sensible_heat_flux=sensible_heat_flux,
        rain_heat_flux=rain_heat_flux,
    )
end

"""
    _bare_ice_ablation_mass(c, forcing, dt_seconds)

Dispatch bare-ice ablation. Diurnal shortwave corrections are snow-column-only,
so bare ice remains on the daily-mean surface-energy path.
"""
function _bare_ice_ablation_mass(
    c::SnowpackPhysicalConstants,
    forcing::SnowpackStepForcing,
    dt_seconds,
)
    return _bare_ice_surface_mass_fluxes_resolved(
        c,
        forcing.air_temperature,
        forcing.rainfall_rate,
        dt_seconds,
        forcing.shortwave_down,
        forcing.has_q_sw_net,
        forcing.q_sw_net,
        forcing.has_q_lw_down,
        forcing.q_lw_down,
        forcing.has_q_sh,
        forcing.q_sh,
        forcing.has_q_lh,
        forcing.q_lh,
        forcing.has_relative_humidity,
        forcing.relative_humidity,
        forcing.air_pressure,
        forcing.wind_speed,
        _bare_ice_albedo(c, forcing),
    )
end

"""Albedo of an exposed bare-ice surface."""
@inline function _bare_ice_albedo(c::SnowpackPhysicalConstants, forcing::SnowpackStepForcing)
    return _uses_semix_albedo(c) && forcing.has_prescribed_ice_albedo ?
        clamp(forcing.prescribed_ice_albedo, zero(forcing.prescribed_ice_albedo), one(forcing.prescribed_ice_albedo)) :
    _uses_prescribed_albedo(c) && forcing.has_prescribed_albedo ?
        _prescribed_surface_albedo(forcing) :
        c.alpha_ice
end

@inline function _resolved_turbulent_latent_heat_flux(
    c::SnowpackPhysicalConstants,
    surface_temperature,
    air_temperature,
    use_q_lh::Bool,
    q_lh_value,
    use_relative_humidity::Bool,
    relative_humidity,
    air_pressure,
    wind_speed,
)
    if use_q_lh
        return q_lh_value
    elseif !use_relative_humidity
        return zero(surface_temperature)
    elseif _uses_semix_turbulence(c)
        _, _, latent_constant, latent_linear = _semix_turbulent_flux_linearized(
            surface_temperature, c, air_temperature, relative_humidity,
            air_pressure, wind_speed, c.semix_z0m_snow,
        )
        return latent_constant - latent_linear * surface_temperature
    end
    return _bessi_latent_vapor_flux(surface_temperature, c, air_temperature, relative_humidity, air_pressure)
end

"""
    _apply_snow_surface_vapor_mass_flux!(..., forcing, dt_seconds)

Apply the BESSI turbulent latent-heat mass bookkeeping to a snow-covered
surface. Positive latent flux deposits mass; negative latent flux removes mass.
Subfreezing surfaces exchange solid mass using `Lv + Lm`, while melting-point
surfaces exchange liquid water using `Lv`.
"""
function _apply_snow_surface_vapor_mass_flux!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    runoff,
    Tsrf,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    forcing::SnowpackStepForcing,
    dt_seconds,
    mass_split,
    mass_min,
)
    if !_surface_has_snow(N_storage, mass, idx)
        return _surface_vapor_fluxes(zero(dt_seconds), zero(dt_seconds))
    end

    surface_temperature = _get_scalar(Tsrf, idx)
    latent_heat_flux = _resolved_turbulent_latent_heat_flux(
        c,
        surface_temperature,
        forcing.air_temperature,
        forcing.has_q_lh,
        forcing.q_lh,
        forcing.has_relative_humidity,
        forcing.relative_humidity,
        forcing.air_pressure,
        forcing.wind_speed,
    )
    if latent_heat_flux == zero(latent_heat_flux)
        return _surface_vapor_fluxes(zero(latent_heat_flux), latent_heat_flux)
    end

    # A prescribed latent-heat flux must also control its associated mass
    # exchange. This preserves the MAR LHF/SU consistency in all-prescribed
    # runs. For parameterized fluxes, SEMIX converts its resolved energy flux
    # with the phase-appropriate latent heat, while BESSI retains its direct
    # vapour-gradient bookkeeping.
    vapor_mass = forcing.has_q_lh ?
                 latent_heat_flux * dt_seconds / _surface_vapor_latent_heat(surface_temperature, c) :
                 _uses_semix_turbulence(c) ?
                 latent_heat_flux * dt_seconds / _surface_vapor_latent_heat(surface_temperature, c) :
                 _bessi_vapor_mass_flux(
                     surface_temperature,
                     c,
                     forcing.air_temperature,
                     forcing.relative_humidity,
                     forcing.air_pressure,
                 ) * dt_seconds

    if surface_temperature < c.T0
        previous_surface_mass = _get_layer(mass, 1, idx)
        updated_surface_mass = max(previous_surface_mass + vapor_mass, zero(vapor_mass))
        _set_layer!(mass, 1, idx, updated_surface_mass)
        # A strongly sublimating interval can exhaust a thin surface layer.
        # In that case the exchange diagnosed to the atmosphere must be the
        # mass actually available, rather than the unconstrained flux demand.
        vapor_mass = updated_surface_mass - previous_surface_mass
        while _n_active(N_storage, idx) > 0 && _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
            _remove_depleted_surface_and_route_water!(N_storage, mass, mass_w, density, temperature, runoff, idx, c)
        end
        while _n_active(N_storage, idx) > 1 && _get_layer(mass, 1, idx) < mass_min
            _merge_surface_layer!(N_storage, mass, mass_w, density, temperature, idx, typemax(Int), mass_split, mass_min, c)
        end
    else
        previous_surface_water = _get_layer(mass_w, 1, idx)
        updated_surface_water = max(previous_surface_water + vapor_mass, zero(vapor_mass))
        _set_layer!(mass_w, 1, idx, updated_surface_water)
        vapor_mass = updated_surface_water - previous_surface_water
    end

    if _n_active(N_storage, idx) == 0
        _set_scalar!(Tsrf, idx, c.T0)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end
    return _surface_vapor_fluxes(vapor_mass, latent_heat_flux)
end
