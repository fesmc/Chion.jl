"""
Bare-ice and diurnal surface-flux helpers used by `step.jl`.
"""

# Daily-mean top-of-atmosphere shortwave below which the shortwave cloud proxy
# is undefined (polar night) and the night cloud fraction is used instead.
const LONGWAVE_CLOUD_PROXY_MIN_TOA = 50.0
const SOLAR_CONSTANT = 1361.0

"""Daily-mean top-of-atmosphere shortwave (W m⁻²) from latitude and season."""
@inline _daily_toa_shortwave(latitude_deg, solar_longitude_deg, day_of_year) =
    _daily_toa_shortwave(_diurnal_shortwave_integral_terms(latitude_deg, solar_longitude_deg), day_of_year)

@inline function _daily_toa_shortwave(terms::NamedTuple, day_of_year)
    x = terms.daylight_integral
    eccentricity = one(x) + oftype(x, 0.033) * cos(oftype(x, 2π) * day_of_year / oftype(x, 365))
    return oftype(x, SOLAR_CONSTANT) * eccentricity * max(x, zero(x)) / oftype(x, 2π)
end

"""
    clear_sky_transmissivity(params, height, toa)

Daily shortwave transmissivity of a clear sky (Konzelmann et al., 1994) at the
surface height `height` (m) and daily-mean top-of-atmosphere shortwave `toa`
(W m⁻²), τ_cs = min(1, τ₀ + τ_z z + τ_S ln(SW_TOA / 300 W m⁻²)) with z in km,
the reference of the cloud proxy's cloud fraction. A height that is not finite
(no surface height) or negative counts as sea level; `toa` is floored at
$(LONGWAVE_CLOUD_PROXY_MIN_TOA) W m⁻², as in the proxy.
"""
@inline function clear_sky_transmissivity(params::Params, height, toa)
    c = params
    height_km = ifelse(isfinite(height), max(height, zero(height)), zero(height)) / oftype(height, 1000)
    safe_toa = max(toa, oftype(toa, LONGWAVE_CLOUD_PROXY_MIN_TOA))
    transmissivity = c.longwave_scheme.clear_sky_transmissivity_sea_level + c.longwave_scheme.clear_sky_transmissivity_per_km * height_km +
                     c.longwave_scheme.clear_sky_transmissivity_insolation * log(safe_toa / oftype(safe_toa, 300))
    return min(one(transmissivity), transmissivity)
end

"""Daily cloud fraction diagnosed from the shortwave transmissivity of the atmosphere.

The clear-sky transmissivity τ_cs (`clear_sky_transmissivity`) increases with
elevation and insolation, and clouds reduce it by the factor
1 − a nᵖ exp(−b z) (Konzelmann et al., 1994), which is inverted for n ∈ [0, 1].
In the polar night, or for sub-daily forcing where SW↓ is not a daily mean,
n is the constant night cloud fraction.
"""
@inline function _cloud_proxy_cloud_fraction(c::Params, forcing::SnowpackStepForcing, terms)
    T = forcing.air_temperature
    toa = _daily_toa_shortwave(terms, forcing.day_of_year)
    # Forcing without surface height (NaN) uses the sea-level coefficients.
    height_km = ifelse(isfinite(forcing.surface_height), max(forcing.surface_height, zero(T)), zero(T)) / oftype(T, 1000)
    daily = (toa > oftype(T, LONGWAVE_CLOUD_PROXY_MIN_TOA)) & (forcing.dt_days >= oftype(T, 0.75)) &
            isfinite(forcing.latitude_deg)
    safe_toa = max(toa, oftype(T, LONGWAVE_CLOUD_PROXY_MIN_TOA))
    clear_transmissivity = clear_sky_transmissivity(c, forcing.surface_height, toa)
    deficit = clamp(one(T) - forcing.shortwave_down / _safe_positive(safe_toa * clear_transmissivity), zero(T), one(T))
    cloud_effect = c.longwave_scheme.cloud_transmission_coefficient * exp(-c.longwave_scheme.cloud_transmission_decay_per_km * height_km)
    cloud_fraction = _safe_pow(clamp(deficit / _safe_positive(cloud_effect), zero(T), one(T)),
        one(T) / c.longwave_scheme.cloud_transmission_exponent)
    return ifelse(daily, cloud_fraction, c.longwave_scheme.night_cloud_fraction)
end

"""
    _cloud_proxy_emissivity(c, forcing)

Effective atmospheric emissivity relative to the near-surface air temperature,
weighted between clear-sky and overcast emissivities (Konzelmann et al., 1994),
ε = ε_cs(T)(1 − n^q) + ε_oc(T) n^q, with both emissivities linear in T − T₀ and
the cloud fraction n from `_cloud_proxy_cloud_fraction`.
"""
@inline _cloud_proxy_emissivity(c::Params, forcing::SnowpackStepForcing) =
    _cloud_proxy_emissivity(c, forcing, _diurnal_shortwave_integral_terms(forcing.latitude_deg, forcing.solar_longitude_deg))

@inline function _cloud_proxy_emissivity(c::Params, forcing::SnowpackStepForcing, terms)
    T = forcing.air_temperature
    weight = _safe_pow(_cloud_proxy_cloud_fraction(c, forcing, terms), c.longwave_scheme.cloud_emissivity_exponent)
    clear = c.longwave_scheme.clear_sky_emissivity + c.longwave_scheme.clear_sky_emissivity_temperature_slope * (T - c.T0)
    overcast = c.longwave_scheme.overcast_emissivity + c.longwave_scheme.overcast_emissivity_temperature_slope * (T - c.T0)
    return clamp(clear * (one(T) - weight) + overcast * weight, oftype(T, 0.4), c.longwave_scheme.emissivity_max)
end

"""
    _with_parameterized_longwave(c, forcing, terms)

For the `CloudProxyLongwave` scheme, resolve the daily longwave down once from
the daily forcing and pass it on as a fixed flux, so the diurnal substeps reuse
the daily cloud proxy. Prescribed longwave and the graybody scheme pass through.
"""
@inline _with_parameterized_longwave(c::Params, forcing::SnowpackStepForcing, terms) =
    _with_parameterized_longwave(c.longwave_scheme, c, forcing, terms)

@inline _with_parameterized_longwave(::GraybodyLongwave, c, forcing::SnowpackStepForcing, terms) = forcing

@inline function _with_parameterized_longwave(::CloudProxyLongwave, c, forcing::SnowpackStepForcing, terms)
    forcing.has_q_lw_down && return forcing
    q_lw_down = _cloud_proxy_emissivity(c, forcing, terms) * c.σ * forcing.air_temperature^4
    return _with_step_forcing(forcing, (; q_lw_down, has_q_lw_down=true))
end

"""Longwave down (W m⁻²): the forcing's, or the graybody sky's ϵ_air σ T_a⁴."""
@inline _longwave_down(c::Params, air_temperature, use_q_lw_down::Bool, q_lw_down_value) =
    use_q_lw_down ? q_lw_down_value : c.σ * c.ϵ_air * air_temperature^4

"""
    _net_longwave(scheme, c, emissivity, longwave_down, surface_temperature)

Net longwave flux into the surface (W m⁻²) at the surface temperature `T_s` (K).
SEMIX: ε (L↓ − σ T_s⁴); BESSI: L↓ − ε_snow σ T_s⁴.
"""
@inline _net_longwave(::SemixSurfaceEnergyBalance, c, emissivity, longwave_down, surface_temperature) =
    emissivity * (longwave_down - c.σ * surface_temperature^4)
@inline _net_longwave(::BESSISurfaceEnergyBalance, c, emissivity, longwave_down, surface_temperature) =
    longwave_down - c.σ * c.ϵ_snow * surface_temperature^4

"""
    _roughness_length(scheme, surface_is_ice)

Momentum roughness length z₀ₘ (m) of the snow or bare-ice surface for the turbulent fluxes. Only
`SemixTurbulentFluxes` uses it; the BESSI fluxes have none (zero, unused).
"""
@inline _roughness_length(s::SemixTurbulentFluxes, surface_is_ice::Bool) = surface_is_ice ? s.z0m_ice : s.z0m_snow
@inline _roughness_length(s::BESSITurbulentFluxes, surface_is_ice::Bool) = zero(s.D_sh)

"""
    _turbulent_heat_fluxes(scheme, c, surface_temperature, air_temperature, relative_humidity, air_pressure, wind_speed, z0m)

Parameterized sensible and latent heat fluxes into the surface (W m⁻²) at the surface temperature
(K). SEMIX: the bulk formulation with stability (`_semix_turbulent_flux_linearized`); BESSI:
`D_sh (T_a − T_s)` and the vapour-gradient latent heat (`_bessi_latent_vapor_flux`).
"""
@inline function _turbulent_heat_fluxes(::SemixTurbulentFluxes, c, surface_temperature, air_temperature,
        relative_humidity, air_pressure, wind_speed, z0m)
    sensible_constant, sensible_linear, latent_constant, latent_linear = _semix_turbulent_flux_linearized(
        surface_temperature, c, air_temperature, relative_humidity, air_pressure, wind_speed, z0m,
    )
    return (
        sensible_constant - sensible_linear * surface_temperature,
        latent_constant - latent_linear * surface_temperature,
    )
end

@inline function _turbulent_heat_fluxes(::BESSITurbulentFluxes, c, surface_temperature, air_temperature,
        relative_humidity, air_pressure, wind_speed, z0m)
    return (
        c.turbulent_flux_scheme.D_sh * (air_temperature - surface_temperature),
        _bessi_latent_vapor_flux(surface_temperature, c, air_temperature, relative_humidity, air_pressure),
    )
end

"""
    _bare_ice_latent_heat_flux(scheme, c, latent_heat_flux, use_q_lh, use_relative_humidity, air_temperature, relative_humidity, air_pressure)

Latent heat flux over bare ice (W m⁻²). Bare ice remains a solid surface at the melting point: its
direct vapour exchange is sublimation/deposition and therefore carries Lᵥ + Lₘ. BESSI evaluates its
vapour-gradient flux with that latent heat; SEMIX keeps its bulk flux.
"""
@inline _bare_ice_latent_heat_flux(::SemixTurbulentFluxes, c, latent_heat_flux, use_q_lh::Bool,
    use_relative_humidity::Bool, air_temperature, relative_humidity, air_pressure) = latent_heat_flux

@inline function _bare_ice_latent_heat_flux(::BESSITurbulentFluxes, c, latent_heat_flux, use_q_lh::Bool,
        use_relative_humidity::Bool, air_temperature, relative_humidity, air_pressure)
    return !use_q_lh && use_relative_humidity ?
        _bessi_latent_vapor_flux(c.T0, c, air_temperature, relative_humidity, air_pressure, c.Lv + c.Lm) :
        latent_heat_flux
end

"""
    _parameterized_vapor_mass(scheme, c, latent_heat_flux, dt_seconds, surface_temperature, air_temperature, relative_humidity, air_pressure)

Vapour mass exchanged with the atmosphere over `dt_seconds` (kg m⁻², positive = deposition) for
parameterized fluxes. SEMIX converts its resolved energy flux with the phase-appropriate latent
heat; BESSI keeps its direct vapour-gradient bookkeeping.
"""
@inline _parameterized_vapor_mass(::SemixTurbulentFluxes, c, latent_heat_flux, dt_seconds, surface_temperature,
    air_temperature, relative_humidity, air_pressure) =
    latent_heat_flux * dt_seconds / _surface_vapor_latent_heat(surface_temperature, c)

@inline _parameterized_vapor_mass(::BESSITurbulentFluxes, c, latent_heat_flux, dt_seconds, surface_temperature,
    air_temperature, relative_humidity, air_pressure) =
    _bessi_vapor_mass_flux(surface_temperature, c, air_temperature, relative_humidity, air_pressure) * dt_seconds

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
    c::Params,
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
    longwave_down = _longwave_down(c, air_temperature, use_q_lw_down, q_lw_down_value)
    longwave_flux = _net_longwave(c.surface_energy_balance_scheme, c, emissivity, longwave_down, surface_temperature)
    parameterized_sensible, parameterized_latent = _turbulent_heat_fluxes(
        c.turbulent_flux_scheme, c, surface_temperature, air_temperature, relative_humidity,
        air_pressure, wind_speed, z0m,
    )
    sensible_heat_flux = use_q_sh ? q_sh_value : parameterized_sensible
    latent_heat_flux = use_q_lh ? q_lh_value : use_relative_humidity ? parameterized_latent : zero(dt_seconds)
    rain_heat_flux = rainfall_rate * c.cw * (air_temperature - c.T0)
    return longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux
end

"""
    _resolved_bare_ice_surface_flux_components(c, air_temperature, rainfall_rate, dt_seconds, shortwave_down, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Resolve all bare-ice surface-flux components, including absorbed shortwave
energy.
"""
@inline function _resolved_bare_ice_surface_flux_components(
    c::Params,
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
            _roughness_length(c.turbulent_flux_scheme, true),
            c.eps_ice,
        )
    # (Snow with liquid water uses Lᵥ in the snow-surface routine.)
    latent_heat_flux = _bare_ice_latent_heat_flux(
        c.turbulent_flux_scheme, c, latent_heat_flux, use_q_lh, use_relative_humidity,
        air_temperature, relative_humidity, air_pressure,
    )
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
    c::Params,
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
    c::Params,
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
@inline _bare_ice_albedo(c::Params, forcing::SnowpackStepForcing) =
    _bare_ice_albedo(c.albedo_scheme, c, forcing)

@inline _bare_ice_albedo(::SemixAlbedo, c::Params, forcing::SnowpackStepForcing) =
    forcing.has_prescribed_ice_albedo ?
        clamp(forcing.prescribed_ice_albedo, zero(forcing.prescribed_ice_albedo), one(forcing.prescribed_ice_albedo)) :
        c.alpha_ice
@inline _bare_ice_albedo(::PrescribedAlbedo, c::Params, forcing::SnowpackStepForcing) =
    forcing.has_prescribed_albedo ? _prescribed_surface_albedo(forcing) : c.alpha_ice
@inline _bare_ice_albedo(::AbstractAlbedoScheme, c::Params, forcing::SnowpackStepForcing) =
    c.alpha_ice

@inline function _resolved_turbulent_latent_heat_flux(
    c::Params,
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
    end
    _, latent_heat_flux = _turbulent_heat_fluxes(
        c.turbulent_flux_scheme, c, surface_temperature, air_temperature, relative_humidity,
        air_pressure, wind_speed, _roughness_length(c.turbulent_flux_scheme, false),
    )
    return latent_heat_flux
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
    c::Params,
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
                 _parameterized_vapor_mass(
                     c.turbulent_flux_scheme, c, latent_heat_flux, dt_seconds, surface_temperature,
                     forcing.air_temperature, forcing.relative_humidity, forcing.air_pressure,
                 )

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

    # Sublimation that removes the last layer leaves the (subfreezing) surface temperature
    # just solved. Resetting it to T₀ would make the next step's vapour flux, linearised about
    # T₀, sublimate any fresh snowfall at once: on dry, cold surfaces (East Antarctic plateau)
    # the column then never builds snow again.
    if _n_active(N_storage, idx) == 0
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end
    return _surface_vapor_fluxes(vapor_mass, latent_heat_flux)
end
