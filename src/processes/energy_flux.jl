#=
Energy-flux temperature solver for array-backed snowpack states.
=#

@inline _safe_positive(x) = x > EPS_TINY ? x : oftype(x, EPS_TINY)

@inline function _copy_column!(dst::AbstractMatrix, src::AbstractMatrix, idx::Int, n::Int)
    @inbounds for i in 1:n
        dst[i, idx] = src[i, idx]
    end
    return nothing
end
@inline function _copy_column!(dst::AbstractVector, src::AbstractVector, ::Int, n::Int)
    @inbounds @simd for i in 1:n
        dst[i] = src[i]
    end
    return nothing
end

"""
    _clamp_to_melt!(temperature_profile, idx, melting_temperature, n)

Clamp the first `n` layer temperatures of column `idx` to `melting_temperature`
from above. Mutates `temperature_profile` in-place and returns it.
"""
@inline function _clamp_to_melt!(
    temperature_profile,
    idx::Int,
    melting_temperature,
    n::Int,
)
    @inbounds @simd for layer_index in 1:n
        if _get_layer(temperature_profile, layer_index, idx) > melting_temperature
            _set_layer!(temperature_profile, layer_index, idx, melting_temperature)
        end
    end
    return temperature_profile
end

@inline function _snow_thermal_conductivity(ρ, temperature, ice_density)
    # Calonne et al. (2019), Eq. (5): a smooth blend of the snow and firn
    # regressions, scaled by the temperature dependences of ice and air.
    reference_ice_conductivity = oftype(ρ, 2.107)
    reference_air_conductivity = oftype(ρ, 0.024)
    transition = one(ρ) / (one(ρ) + exp(-oftype(ρ, 0.04) * (ρ - oftype(ρ, 450.0))))
    ice_conductivity = oftype(ρ, 9.828) * exp(-oftype(ρ, 5.7e-3) * temperature)
    air_conductivity = oftype(ρ, 2.334e-3) * temperature^oftype(ρ, 1.5) /
                       _safe_positive(oftype(ρ, 164.54) + temperature)
    snow_conductivity = reference_air_conductivity -
                         oftype(ρ, 1.23e-4) * ρ +
                         oftype(ρ, 2.5e-6) * ρ^2
    firn_conductivity = reference_ice_conductivity +
                         oftype(ρ, 3.618e-3) * (ρ - ice_density)
    snow_scale = ice_conductivity * air_conductivity /
                 (reference_ice_conductivity * reference_air_conductivity)
    firn_scale = ice_conductivity / reference_ice_conductivity
    return (one(ρ) - transition) * snow_scale * snow_conductivity +
           transition * firn_scale * firn_conductivity
end

"""Bulk turbulent vapour-mass transfer coefficient (kg m⁻² s⁻¹ Pa⁻¹)."""
@inline _bessi_vapor_exchange_coefficient(c::SnowpackPhysicalConstants) =
    c.latent_heat_flux_ratio * c.D_sh / c.cp_air * oftype(c.D_sh, 0.622)


"""Latent heat for vapour exchange at a snow surface of the given temperature."""
@inline _surface_vapor_latent_heat(surface_temperature, c::SnowpackPhysicalConstants) =
    surface_temperature < c.T0 ? c.Lv + c.Lm : c.Lv

@inline function _relative_humidity_fraction(relative_humidity)
    rh = relative_humidity > one(relative_humidity) ? relative_humidity / oftype(relative_humidity, 100.0) : relative_humidity
    return clamp(rh, zero(rh), one(rh))
end

@inline function _bessi_water_saturation_vapor_pressure(temperature, T0)
    temperature_c = temperature - T0
    return oftype(temperature, 611.2) *
           exp(oftype(temperature, 17.27) * temperature_c / (temperature_c + oftype(temperature, 243.12)))
end

@inline function _bessi_air_vapor_pressure(air_temperature, relative_humidity, T0)
    return _relative_humidity_fraction(relative_humidity) *
           _bessi_water_saturation_vapor_pressure(air_temperature, T0)
end

@inline function _bessi_ice_saturation_vapor_pressure(surface_temperature, T0)
    surface_c = surface_temperature - T0
    return oftype(surface_temperature, 611.2) *
           exp(oftype(surface_temperature, 22.46) * surface_c / (surface_c + oftype(surface_temperature, 272.62)))
end

@inline function _bessi_ice_saturation_vapor_pressure_derivative(surface_temperature, T0, es)
    denominator = surface_temperature - T0 + oftype(surface_temperature, 272.62)
    return es * oftype(surface_temperature, 22.46) * oftype(surface_temperature, 272.62) / _safe_positive(denominator * denominator)
end

@inline function _bessi_vapor_mass_flux(surface_temperature, c::SnowpackPhysicalConstants, air_temperature, relative_humidity, air_pressure)
    exchange = _bessi_vapor_exchange_coefficient(c) / _safe_positive(air_pressure)
    ea = _bessi_air_vapor_pressure(air_temperature, relative_humidity, c.T0)
    es = _bessi_ice_saturation_vapor_pressure(surface_temperature, c.T0)
    return exchange * (ea - es)
end

@inline function _bessi_latent_vapor_flux(surface_temperature, c::SnowpackPhysicalConstants, air_temperature, relative_humidity, air_pressure, latent_heat=_surface_vapor_latent_heat(surface_temperature, c))
    vapor_mass_flux = _bessi_vapor_mass_flux(surface_temperature, c, air_temperature, relative_humidity, air_pressure)
    return latent_heat * vapor_mass_flux
end

@inline function _bessi_latent_vapor_flux_linearized(surface_temperature, c::SnowpackPhysicalConstants, air_temperature, relative_humidity, air_pressure)
    exchange = _bessi_vapor_exchange_coefficient(c) / _safe_positive(air_pressure)
    ea = _bessi_air_vapor_pressure(air_temperature, relative_humidity, c.T0)
    es = _bessi_ice_saturation_vapor_pressure(surface_temperature, c.T0)
    des_dT = _bessi_ice_saturation_vapor_pressure_derivative(surface_temperature, c.T0, es)
    latent_heat = _surface_vapor_latent_heat(surface_temperature, c)
    linear = latent_heat * exchange * des_dT
    constant = latent_heat * exchange * (ea - es + des_dT * surface_temperature)
    return constant, linear
end

# SEMIX uses a neutral aerodynamic transfer coefficient corrected with the
# bulk Richardson number. Positive Ri denotes warm air over a colder surface
# (stable stratification) and suppresses turbulent exchange; negative Ri
# denotes an unstable surface layer and enhances it.
@inline function _semix_aerodynamic_resistance(c::SnowpackPhysicalConstants, surface_temperature, air_temperature, air_pressure, wind_speed, z0m)
    z0h = z0m / c.semix_zm_to_zh
    neutral_ch = c.semix_karman^2 /
                 _safe_positive(log(c.semix_surface_height / z0m) * log(c.semix_surface_height / z0h))
    wind = max(wind_speed, oftype(wind_speed, 0.1))
    bulk_richardson = oftype(surface_temperature, 9.80665) * c.semix_surface_height *
                      (air_temperature - surface_temperature) /
                      _safe_positive(air_temperature * wind^2)
    stability_factor = ifelse(
        bulk_richardson < zero(bulk_richardson),
        sqrt(max(one(bulk_richardson) - oftype(bulk_richardson, 16) * bulk_richardson, one(bulk_richardson))),
        one(bulk_richardson) /
        (one(bulk_richardson) + c.semix_stable_coefficient * bulk_richardson),
    )
    return one(wind) / _safe_positive(neutral_ch * stability_factor * wind)
end

@inline _semix_air_density(air_temperature, air_pressure) =
    air_pressure / _safe_positive(oftype(air_temperature, 287.05) * air_temperature)

@inline function _semix_turbulent_flux_linearized(surface_temperature, c::SnowpackPhysicalConstants, air_temperature, relative_humidity, air_pressure, wind_speed, z0m)
    air_density = _semix_air_density(air_temperature, air_pressure)
    resistance = _semix_aerodynamic_resistance(c, surface_temperature, air_temperature, air_pressure, wind_speed, z0m)
    sensible_coefficient = c.semix_sensible_exchange_factor * air_density * c.cp_air / resistance
    sensible_constant = sensible_coefficient * air_temperature
    q_air = oftype(surface_temperature, 0.622) * _relative_humidity_fraction(relative_humidity) *
            _bessi_ice_saturation_vapor_pressure(air_temperature, c.T0) / _safe_positive(air_pressure)
    q_surface_pressure = _bessi_ice_saturation_vapor_pressure(surface_temperature, c.T0)
    q_surface = oftype(surface_temperature, 0.622) * q_surface_pressure / _safe_positive(air_pressure)
    dq_surface = oftype(surface_temperature, 0.622) *
                 _bessi_ice_saturation_vapor_pressure_derivative(surface_temperature, c.T0, q_surface_pressure) /
                 _safe_positive(air_pressure)
    # Permit both sublimation and deposition. At a melting/wet surface this
    # carries L_v; otherwise exchange with the solid surface carries L_s.
    latent_exchange = c.semix_latent_exchange_factor *
                      _surface_vapor_latent_heat(surface_temperature, c) * air_density / resistance
    latent_constant = latent_exchange * (q_air - q_surface + dq_surface * surface_temperature)
    latent_linear = latent_exchange * dq_surface
    return sensible_constant, sensible_coefficient, latent_constant, latent_linear
end

"""
    shortwave_absorbed(shortwave_down; surface_albedo)

Return the absorbed shortwave flux after applying `surface_albedo`.
"""
@inline function shortwave_absorbed(
    shortwave_down;
    surface_albedo,
)
    return (one(shortwave_down) - clamp(surface_albedo, zero(surface_albedo), one(surface_albedo))) * shortwave_down
end
@inline shortwave_absorbed(shortwave_down, surface_albedo) =
    (one(shortwave_down) - clamp(surface_albedo, zero(surface_albedo), one(surface_albedo))) * shortwave_down

@inline function interface_conductance(Kᵢ, Δzᵢ, Kⱼ, Δzⱼ)
    # Thermal resistances over the two half layers are in series:
    # G = 1 / (Δzᵢ / (2Kᵢ) + Δzⱼ / (2Kⱼ)).
    return oftype(Kᵢ, 2.0) * Kᵢ * Kⱼ /
           _safe_positive(Kⱼ * Δzᵢ + Kᵢ * Δzⱼ)
end

"""
    _thomas_forward!(lower_diagonal, main_diagonal, upper_diagonal, right_hand_side, idx, n)

Forward elimination phase of the Thomas algorithm for column `idx` over the
first `n` rows. Mutates `main_diagonal` and `right_hand_side` in-place.
"""
function _thomas_forward!(
    lower_diagonal,
    main_diagonal,
    upper_diagonal,
    right_hand_side,
    idx::Int,
    n::Int,
)
    @inbounds for row_index in 2:n
        elimination_factor = _get_layer(lower_diagonal, row_index - 1, idx) / _get_layer(main_diagonal, row_index - 1, idx)
        _set_layer!(
            main_diagonal,
            row_index,
            idx,
            _get_layer(main_diagonal, row_index, idx) - elimination_factor * _get_layer(upper_diagonal, row_index - 1, idx),
        )
        _set_layer!(
            right_hand_side,
            row_index,
            idx,
            _get_layer(right_hand_side, row_index, idx) - elimination_factor * _get_layer(right_hand_side, row_index - 1, idx),
        )
    end
    return nothing
end

"""
    _thomas_backward!(main_diagonal, upper_diagonal, right_hand_side, idx, n)

Back-substitution phase of the Thomas algorithm for column `idx` over the
first `n` rows. The solution overwrites `right_hand_side`.
"""
function _thomas_backward!(
    main_diagonal,
    upper_diagonal,
    right_hand_side,
    idx::Int,
    n::Int,
)
    _set_layer!(
        right_hand_side,
        n,
        idx,
        _get_layer(right_hand_side, n, idx) / _get_layer(main_diagonal, n, idx),
    )
    @inbounds for offset in Base.OneTo(max(n - 1, 0))
        row_index = n - offset
        _set_layer!(
            right_hand_side,
            row_index,
            idx,
            (
                _get_layer(right_hand_side, row_index, idx) -
                _get_layer(upper_diagonal, row_index, idx) * _get_layer(right_hand_side, row_index + 1, idx)
            ) / _get_layer(main_diagonal, row_index, idx),
        )
    end
    return right_hand_side
end

"""
    _solve_tridiagonal_thomas_prefix!(lower_diagonal, main_diagonal, upper_diagonal, right_hand_side, idx, n)

Solve an in-place tridiagonal system for column `idx` using the Thomas
algorithm on the first `n` rows. The solution overwrites `right_hand_side`.
"""
@inline function _solve_tridiagonal_thomas_prefix!(
    lower_diagonal,
    main_diagonal,
    upper_diagonal,
    right_hand_side,
    idx::Int,
    n::Int,
)
    _thomas_forward!(lower_diagonal, main_diagonal, upper_diagonal, right_hand_side, idx, n)
    return _thomas_backward!(main_diagonal, upper_diagonal, right_hand_side, idx, n)
end

"""
    _energy_flux_result(; ...)

Package the key diagnostics returned by the energy-flux solver into a named
tuple.
"""
@inline _energy_flux_result(;
    needs_melt,
    melt_energy_available,
    heating,
    surface_flux_constant,
    surface_flux_linear,
    longwave_flux_constant,
    longwave_flux_linear,
    latent_heat_linear_coefficient,
    latent_heat_constant_term,
    absorbed_shortwave,
    sensible_heat_flux_constant,
    sensible_heat_flux_linear,
) = (; needs_melt, melt_energy_available, heating, surface_flux_constant, surface_flux_linear,
    longwave_flux_constant, longwave_flux_linear, latent_heat_linear_coefficient,
    latent_heat_constant_term, absorbed_shortwave, sensible_heat_flux_constant, sensible_heat_flux_linear)

"""
    _diagnose_latent_heat_flux_coefficients(has_surface_snow, c, air_temperature, snowfall_rate, rainfall_rate)

Return effective linear and constant latent-heat terms induced by snowfall or
rainfall forcing.
"""
@inline function _diagnose_latent_heat_flux_coefficients(
    has_surface_snow::Bool,
    c::SnowpackPhysicalConstants,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
)
    if snowfall_rate > zero(snowfall_rate)
        return snowfall_rate * c.ci, snowfall_rate * c.ci * air_temperature
    elseif has_surface_snow && rainfall_rate > zero(rainfall_rate)
        return zero(rainfall_rate), rainfall_rate * c.cw * (air_temperature - c.T0)
    else
        return zero(air_temperature), zero(air_temperature)
    end
end

"""
    _surface_has_snow(N_storage, mass, idx)

Return `true` when column `idx` has a nonempty surface snow layer.
"""
@inline function _surface_has_snow(N_storage, mass, idx::Int)
    return _n_active(N_storage, idx) > 0 && _get_layer(mass, 1, idx) > EPS_EMPTY_LAYER
end

"""
    _thermal_row(mass, density, temperature, ice_temperature, n_snow, ice_top_thickness, ice_density, idx, row)

Return `(mass, density, temperature)` of one thermal row: rows `1:n_snow` are
the snow/firn layers and the remaining rows are the fixed-geometry ice
substrate, whose layer thicknesses double from `ice_top_thickness`.
"""
@inline function _thermal_row(mass, density, temperature, ice_temperature, n_snow::Int, ice_top_thickness, ice_density, idx::Int, row::Int)
    if row <= n_snow
        return _get_layer(mass, row, idx), _get_layer(density, row, idx), _get_layer(temperature, row, idx)
    end
    ice_row = row - n_snow
    return ice_density * ldexp(ice_top_thickness, ice_row - 1), ice_density, _get_layer(ice_temperature, ice_row, idx)
end

@inline function _set_thermal_row_temperature!(temperature, ice_temperature, n_snow::Int, idx::Int, row::Int, value)
    if row <= n_snow
        _set_layer!(temperature, row, idx, value)
    else
        _set_layer!(ice_temperature, row - n_snow, idx, value)
    end
    return nothing
end

"""
    _go_energy_flux_resolved!(..., scratch, air_temperature, shortwave_down, latent_heat_linear_coefficient_eff, latent_heat_constant_term_eff, dt_seconds, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value, ..., ice_temperature, n_ice, ice_top_thickness)

Advance the temperature profile of column `idx` over one time step by solving
the surface energy balance and vertical heat diffusion problem. Mutates
temperature and surface-temperature state in-place and returns a named tuple of
energy diagnostics, including melt availability. With `n_ice > 0` the optional
ice substrate is solved below the snow/firn layers in the same implicit system;
without surface snow its top layer forms the (bare-ice) surface.
"""
function _go_energy_flux_resolved!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    Tsrf,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    scratch,
    air_temperature,
    shortwave_down,
    latent_heat_linear_coefficient_eff,
    latent_heat_constant_term_eff,
    dt_seconds,
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
    ice_temperature=temperature,
    n_ice::Int=0,
    ice_top_thickness=one(dt_seconds),
)
    n_active = _n_active(N_storage, idx)
    n_snow = n_active > 0 && _get_layer(mass, 1, idx) > EPS_EMPTY_LAYER ? n_active : 0
    n_layers = n_snow + n_ice
    if n_layers <= 0
        return _energy_flux_result(;
            needs_melt=false, melt_energy_available=zero(dt_seconds), heating=zero(dt_seconds),
            surface_flux_constant=zero(dt_seconds), surface_flux_linear=zero(dt_seconds),
            longwave_flux_constant=zero(dt_seconds), longwave_flux_linear=zero(dt_seconds),
            latent_heat_linear_coefficient=latent_heat_linear_coefficient_eff,
            latent_heat_constant_term=latent_heat_constant_term_eff,
            absorbed_shortwave=zero(dt_seconds), sensible_heat_flux_constant=zero(dt_seconds), sensible_heat_flux_linear=zero(dt_seconds),
        )
    end
    surface_is_ice = n_snow == 0
    emissivity = surface_is_ice ? c.eps_ice : c.ϵ_snow
    z0m = surface_is_ice ? c.semix_z0m_ice : c.semix_z0m_snow

    lower = scratch.lower
    diag = scratch.diag
    upper = scratch.upper
    rhs = scratch.rhs
    interface_terms = scratch.interface_conductance
    solver_diag = scratch.previous_temperature
    first_mass, first_density, first_temperature =
        _thermal_row(mass, density, temperature, ice_temperature, n_snow, ice_top_thickness, c.rho_i, idx, 1)
    surface_mass = _safe_positive(first_mass)
    # `Tsrf` is the physical atmosphere--snow interface; the first numerical
    # temperature is at the centre of its finite-volume snow cell.
    previous_surface_temperature = _get_scalar(Tsrf, idx)
    surface_temperature_sq = previous_surface_temperature * previous_surface_temperature
    surface_temperature_cube = surface_temperature_sq * previous_surface_temperature
    surface_temperature_fourth = surface_temperature_sq * surface_temperature_sq

    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        shortwave_absorbed(shortwave_down, _get_scalar(albedo_dynamic, idx))
    longwave_down = use_q_lw_down ? q_lw_down_value : c.σ * c.ϵ_air * air_temperature^4
    longwave_flux_constant = if _uses_semix_seb(c)
        emissivity * (longwave_down + c.σ * oftype(air_temperature, 3.0) * surface_temperature_fourth)
    else
        longwave_down + c.σ * c.ϵ_snow * oftype(air_temperature, 3.0) * surface_temperature_fourth
    end
    longwave_flux_linear = c.σ * (_uses_semix_seb(c) ? emissivity : c.ϵ_snow) *
                           oftype(air_temperature, 4.0) * surface_temperature_cube
    semix_sensible_constant, semix_sensible_linear, semix_latent_constant, semix_latent_linear =
        _semix_turbulent_flux_linearized(
            previous_surface_temperature, c, air_temperature, relative_humidity,
            air_pressure, wind_speed, z0m,
        )
    sensible_heat_flux_constant = use_q_sh ? q_sh_value :
                                  _uses_semix_turbulence(c) ? semix_sensible_constant : air_temperature * c.D_sh
    sensible_heat_flux_linear = use_q_sh ? zero(dt_seconds) :
                                _uses_semix_turbulence(c) ? semix_sensible_linear : c.D_sh
    turbulent_latent_heat_constant, turbulent_latent_heat_linear = if use_q_lh
        q_lh_value, zero(dt_seconds)
    elseif use_relative_humidity
        _uses_semix_turbulence(c) ?
        (semix_latent_constant, semix_latent_linear) :
        _bessi_latent_vapor_flux_linearized(previous_surface_temperature, c, air_temperature, relative_humidity, air_pressure)
    else
        zero(dt_seconds), zero(dt_seconds)
    end
    latent_heat_flux_constant = latent_heat_constant_term_eff + turbulent_latent_heat_constant
    latent_heat_flux_linear = latent_heat_linear_coefficient_eff + turbulent_latent_heat_linear

    surface_flux_constant = sensible_heat_flux_constant + longwave_flux_constant + absorbed_shortwave + latent_heat_flux_constant
    surface_flux_linear = sensible_heat_flux_linear + longwave_flux_linear + latent_heat_flux_linear
    needs_melt = false
    heating = zero(dt_seconds)

    previous_layer_thickness = surface_mass / _safe_positive(first_density)
    previous_layer_conductivity = _snow_thermal_conductivity(first_density, first_temperature, c.rho_i)
    surface_layer_thickness = previous_layer_thickness
    surface_layer_conductivity = previous_layer_conductivity
    _set_layer!(rhs, 1, idx, first_temperature)

    @inbounds for layer_index in 2:n_layers
        layer_mass, layer_density, layer_temperature =
            _thermal_row(mass, density, temperature, ice_temperature, n_snow, ice_top_thickness, c.rho_i, idx, layer_index)
        layer_thickness = layer_mass / _safe_positive(layer_density)
        layer_conductivity = _snow_thermal_conductivity(layer_density, layer_temperature, c.rho_i)
        _set_layer!(
            interface_terms,
            layer_index - 1,
            idx,
            interface_conductance(
                previous_layer_conductivity,
                previous_layer_thickness,
                layer_conductivity,
                layer_thickness,
            ),
        )
        _set_layer!(rhs, layer_index, idx, layer_temperature)
        previous_layer_thickness = layer_thickness
        previous_layer_conductivity = layer_conductivity
    end

    # `interface_conductance` is the physical conductance between layer
    # centres (the two half-layer resistances in series).  The finite-volume
    # energy balance therefore contributes -Δt G / (cᵢ mᵢ) off diagonal;
    # the legacy factor of two belonged to the former half-conductance
    # approximation and would double vertical heat diffusion here.
    β_scale = -dt_seconds / c.ci

    # Robin surface boundary: Q_SEB(Ts) = Gs * (Ts - T1). Eliminating the
    # algebraic interface temperature Ts leaves no numerical surface heat
    # capacity; the first layer remains a regular finite-volume snow cell.
    surface_conductance = oftype(surface_mass, 2) * surface_layer_conductivity /
                          _safe_positive(surface_layer_thickness)
    surface_denominator = _safe_positive(surface_flux_linear + surface_conductance)
    surface_constant = surface_flux_constant / surface_denominator
    surface_coefficient = surface_conductance / surface_denominator
    β1 = β_scale / surface_mass
    boundary_term = β1 * surface_conductance
    _set_layer!(rhs, 1, idx, _get_layer(rhs, 1, idx) - boundary_term * surface_constant)
    # Evaluate 1 - boundary_term*(1-surface_coefficient) without subtracting
    # two enormous, nearly equal terms for very thin surface layers. In the
    # zero-flux-derivative limit the diagonal must remain exactly one, not zero.
    _set_layer!(diag, 1, idx, one(dt_seconds) -
        boundary_term * (surface_flux_linear / surface_denominator))
    if n_layers > 1
        _set_layer!(upper, 1, idx, β1 * _get_layer(interface_terms, 1, idx))
        _set_layer!(diag, 1, idx, _get_layer(diag, 1, idx) - _get_layer(upper, 1, idx))
        last_mass, _, _ = _thermal_row(mass, density, temperature, ice_temperature, n_snow, ice_top_thickness, c.rho_i, idx, n_layers)
        βn = β_scale / _safe_positive(last_mass)
        _set_layer!(lower, n_layers - 1, idx, βn * _get_layer(interface_terms, n_layers - 1, idx))
        _set_layer!(diag, n_layers, idx, one(dt_seconds) - _get_layer(lower, n_layers - 1, idx))
    end

    @inbounds for layer_index in 2:(n_layers - 1)
        layer_mass, _, _ = _thermal_row(mass, density, temperature, ice_temperature, n_snow, ice_top_thickness, c.rho_i, idx, layer_index)
        βi = β_scale / _safe_positive(layer_mass)
        _set_layer!(lower, layer_index - 1, idx, βi * _get_layer(interface_terms, layer_index - 1, idx))
        _set_layer!(upper, layer_index, idx, βi * _get_layer(interface_terms, layer_index, idx))
        _set_layer!(
            diag,
            layer_index,
            idx,
            one(dt_seconds) - _get_layer(lower, layer_index - 1, idx) - _get_layer(upper, layer_index, idx),
        )
    end

    _copy_column!(solver_diag, diag, idx, n_layers)
    resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, solver_diag, upper, rhs, idx, n_layers)

    resolved_surface_temperature = surface_constant + surface_coefficient * _get_layer(resolved_temperature, 1, idx)
    if resolved_surface_temperature > c.T0
        needs_melt = true

        @inbounds for layer_index in 1:n_layers
            _, _, layer_temperature = _thermal_row(mass, density, temperature, ice_temperature, n_snow, ice_top_thickness, c.rho_i, idx, layer_index)
            _set_layer!(rhs, layer_index, idx, layer_temperature)
        end
        _set_layer!(rhs, 1, idx, _get_layer(rhs, 1, idx) - boundary_term * c.T0)
        _copy_column!(solver_diag, diag, idx, n_layers)
        # Retain the first-to-second-layer conductance from the unconstrained
        # matrix. Only the eliminated Robin-temperature contribution is
        # replaced by the fixed melting boundary.
        _set_layer!(solver_diag, 1, idx, _get_layer(diag, 1, idx) - boundary_term * surface_coefficient)

        resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, solver_diag, upper, rhs, idx, n_layers)
        _clamp_to_melt!(resolved_temperature, idx, c.T0, n_layers)
        resolved_surface_temperature = c.T0
        heating = dt_seconds * (surface_flux_constant - surface_flux_linear * c.T0)
    else
        _clamp_to_melt!(resolved_temperature, idx, c.T0, n_layers)
        heating = dt_seconds * (surface_flux_constant - surface_flux_linear * resolved_surface_temperature)
    end

    @inbounds for layer_index in 1:n_layers
        _set_thermal_row_temperature!(temperature, ice_temperature, n_snow, idx, layer_index,
            _get_layer(resolved_temperature, layer_index, idx))
    end
    _set_scalar!(Tsrf, idx, resolved_surface_temperature)

    melt_energy_available = needs_melt ? max(
        (surface_flux_constant - surface_flux_linear * c.T0 -
         surface_conductance * (c.T0 - _get_layer(resolved_temperature, 1, idx))) * dt_seconds,
        zero(dt_seconds),
    ) : zero(dt_seconds)

    return _energy_flux_result(;
        needs_melt, melt_energy_available, heating,
        surface_flux_constant, surface_flux_linear,
        longwave_flux_constant, longwave_flux_linear,
        latent_heat_linear_coefficient=latent_heat_linear_coefficient_eff,
        latent_heat_constant_term=latent_heat_constant_term_eff,
        absorbed_shortwave, sensible_heat_flux_constant, sensible_heat_flux_linear,
    )
end

"""
    go_energy_flux!(state, idx, air_temperature, shortwave_down, latent_heat_linear_coefficient, latent_heat_constant_term, dt_seconds; scratch, q_sw_net=nothing, q_lw_down=nothing, q_sh=nothing, q_lh=nothing)

Public wrapper for the column energy-flux solve on `state`. Reuses `scratch`
for temporary arrays and returns the same diagnostic named tuple as
`_go_energy_flux_resolved!`.
"""
function go_energy_flux!(
    state,
    idx::Int,
    air_temperature,
    shortwave_down,
    latent_heat_linear_coefficient,
    latent_heat_constant_term,
    dt_seconds;
    scratch,
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
)
    rows = _n_active(state.N, idx) + size(state.ice_temperature, 1)
    size(scratch.rhs, 1) >= rows ||
        error("`scratch` has $(size(scratch.rhs, 1)) rows but the column needs $rows (snow + ice substrate); use `EnergyWorkspace(state)`.")
    return _go_energy_flux_resolved!(
        state.N,
        state.mass,
        state.mass_w,
        state.density,
        state.temperature,
        state.Tsrf,
        state.albedo,
        idx,
        state.c,
        scratch,
        air_temperature,
        shortwave_down,
        latent_heat_linear_coefficient,
        latent_heat_constant_term,
        dt_seconds,
        !isnothing(q_sw_net),
        isnothing(q_sw_net) ? zero(dt_seconds) : q_sw_net,
        !isnothing(q_lw_down),
        isnothing(q_lw_down) ? zero(dt_seconds) : q_lw_down,
        !isnothing(q_sh),
        isnothing(q_sh) ? zero(dt_seconds) : q_sh,
        !isnothing(q_lh),
        isnothing(q_lh) ? zero(dt_seconds) : q_lh,
        false,
        zero(dt_seconds),
        oftype(dt_seconds, 101_325.0),
        oftype(dt_seconds, 5),
        state.ice_temperature,
        size(state.ice_temperature, 1),
        state.parameters.ice_substrate_top_thickness_m,
    )
end
