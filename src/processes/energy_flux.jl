#=
Energy-flux temperature solver for array-backed snowpack states.
=#

@inline _safe_positive(x) = x > EPS_TINY ? x : oftype(x, EPS_TINY)

@inline _copy_column!(dst::AbstractMatrix, src::AbstractMatrix, idx::Int, n::Int) =
    copyto!(view(dst, 1:n, idx), view(src, 1:n, idx))
@inline function _copy_column!(dst::AbstractVector, src::AbstractVector, ::Int, n::Int)
    @inbounds @simd for i in 1:n
        dst[i] = src[i]
    end
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

@inline _snow_thermal_conductivity(ρ, Kᵢ) =
    Kᵢ * (ρ * oftype(ρ, 1.0e-3))^oftype(ρ, 1.88)

@inline _bessi_latent_exchange_coefficient(c::SnowpackPhysicalConstants) =
    c.latent_heat_flux_ratio * c.D_sh / c.cp_air * oftype(c.D_sh, 0.622) * (c.Lv + c.Lm)

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

@inline function _bessi_latent_vapor_flux(surface_temperature, c::SnowpackPhysicalConstants, air_temperature, relative_humidity, air_pressure)
    exchange = _bessi_latent_exchange_coefficient(c) / _safe_positive(air_pressure)
    ea = _bessi_air_vapor_pressure(air_temperature, relative_humidity, c.T0)
    es = _bessi_ice_saturation_vapor_pressure(surface_temperature, c.T0)
    return exchange * (ea - es)
end

@inline function _bessi_latent_vapor_flux_linearized(surface_temperature, c::SnowpackPhysicalConstants, air_temperature, relative_humidity, air_pressure)
    exchange = _bessi_latent_exchange_coefficient(c) / _safe_positive(air_pressure)
    ea = _bessi_air_vapor_pressure(air_temperature, relative_humidity, c.T0)
    es = _bessi_ice_saturation_vapor_pressure(surface_temperature, c.T0)
    des_dT = _bessi_ice_saturation_vapor_pressure_derivative(surface_temperature, c.T0, es)
    linear = exchange * des_dT
    constant = exchange * (ea - es + des_dT * surface_temperature)
    return constant, linear
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

@inline interface_conductance(Kᵢ, Δzᵢ, Kⱼ, Δzⱼ) =
    (Kᵢ * Δzᵢ + Kⱼ * Δzⱼ) / _safe_positive((Δzᵢ + Δzⱼ)^2)

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
    @assert n >= 1
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
    @inbounds for row_index in (n - 1):-1:1
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
@inline function _energy_flux_result(
    ;
    needs_melt,
    energy_to_melting,
    melt_energy_available,
    heating,
    surface_flux_constant,
    surface_flux_linear,
    latent_heat_linear_coefficient,
    latent_heat_constant_term,
)
    return (
        needs_melt=needs_melt,
        energy_to_melting=energy_to_melting,
        melt_energy_available=melt_energy_available,
        heating=heating,
        surface_flux_constant=surface_flux_constant,
        surface_flux_linear=surface_flux_linear,
        latent_heat_linear_coefficient=latent_heat_linear_coefficient,
        latent_heat_constant_term=latent_heat_constant_term,
    )
end
@inline function _energy_flux_result(
    needs_melt,
    energy_to_melting,
    melt_energy_available,
    heating,
    surface_flux_constant,
    surface_flux_linear,
    latent_heat_linear_coefficient,
    latent_heat_constant_term,
)
    return (
        needs_melt=needs_melt,
        energy_to_melting=energy_to_melting,
        melt_energy_available=melt_energy_available,
        heating=heating,
        surface_flux_constant=surface_flux_constant,
        surface_flux_linear=surface_flux_linear,
        latent_heat_linear_coefficient=latent_heat_linear_coefficient,
        latent_heat_constant_term=latent_heat_constant_term,
    )
end

"""
    _residual_melt_energy(surface_flux_constant, surface_flux_linear, surface_temperature, energy_to_melting, dt_seconds; needs_melt)

Return melt energy that remains after bringing the surface to the melting
point.
"""
@inline function _residual_melt_energy(
    surface_flux_constant,
    surface_flux_linear,
    surface_temperature,
    energy_to_melting,
    dt_seconds;
    needs_melt,
)
    needs_melt || return zero(dt_seconds)
    return max(
        (surface_flux_constant - surface_flux_linear * surface_temperature) * dt_seconds - energy_to_melting,
        zero(dt_seconds),
    )
end
@inline function _residual_melt_energy(
    surface_flux_constant,
    surface_flux_linear,
    surface_temperature,
    energy_to_melting,
    dt_seconds,
    needs_melt::Bool,
)
    needs_melt || return zero(dt_seconds)
    return max(
        (surface_flux_constant - surface_flux_linear * surface_temperature) * dt_seconds - energy_to_melting,
        zero(dt_seconds),
    )
end

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
    _go_energy_flux_resolved!(..., scratch, air_temperature, shortwave_down, latent_heat_linear_coefficient_eff, latent_heat_constant_term_eff, dt_seconds, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Advance the temperature profile of column `idx` over one time step by solving
the surface energy balance and vertical heat diffusion problem. Mutates
temperature and surface-temperature state in-place and returns a named tuple of
energy diagnostics, including melt availability.
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
)
    n_layers = _n_active(N_storage, idx)
    if n_layers <= 0 || _get_layer(mass, 1, idx) <= zero(eltype(mass))
        return _energy_flux_result(
            false,
            zero(dt_seconds),
            zero(dt_seconds),
            zero(dt_seconds),
            zero(dt_seconds),
            zero(dt_seconds),
            latent_heat_linear_coefficient_eff,
            latent_heat_constant_term_eff,
        )
    end

    lower = scratch.lower
    diag = scratch.diag
    upper = scratch.upper
    rhs = scratch.rhs
    interface_terms = scratch.interface_conductance
    solver_diag = scratch.previous_temperature
    surface_mass = _safe_positive(_get_layer(mass, 1, idx))
    previous_surface_temperature = _get_layer(temperature, 1, idx)
    surface_temperature_scale = dt_seconds / c.ci / surface_mass
    surface_temperature_sq = previous_surface_temperature * previous_surface_temperature
    surface_temperature_cube = surface_temperature_sq * previous_surface_temperature
    surface_temperature_fourth = surface_temperature_sq * surface_temperature_sq

    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        shortwave_absorbed(shortwave_down, _get_scalar(albedo_dynamic, idx))
    longwave_flux_constant = use_q_lw_down ?
        (q_lw_down_value + c.σ * c.ϵ_snow * oftype(air_temperature, 3.0) * surface_temperature_fourth) :
        (c.σ * (c.ϵ_air * air_temperature^4 + c.ϵ_snow * oftype(air_temperature, 3.0) * surface_temperature_fourth))
    longwave_flux_linear = c.σ * c.ϵ_snow * oftype(air_temperature, 4.0) * surface_temperature_cube
    sensible_heat_flux_constant = use_q_sh ? q_sh_value : air_temperature * c.D_sh
    sensible_heat_flux_linear = use_q_sh ? zero(dt_seconds) : c.D_sh
    turbulent_latent_heat_constant, turbulent_latent_heat_linear = if use_q_lh
        q_lh_value, zero(dt_seconds)
    elseif use_relative_humidity
        _bessi_latent_vapor_flux_linearized(previous_surface_temperature, c, air_temperature, relative_humidity, air_pressure)
    else
        zero(dt_seconds), zero(dt_seconds)
    end
    latent_heat_flux_constant = latent_heat_constant_term_eff + turbulent_latent_heat_constant
    latent_heat_flux_linear = latent_heat_linear_coefficient_eff + turbulent_latent_heat_linear

    surface_flux_constant = sensible_heat_flux_constant + longwave_flux_constant + absorbed_shortwave + latent_heat_flux_constant
    surface_flux_linear = sensible_heat_flux_linear + longwave_flux_linear + latent_heat_flux_linear
    surface_rhs_term = surface_temperature_scale * surface_flux_constant
    surface_diag_term = surface_temperature_scale * surface_flux_linear

    needs_melt = false
    energy_to_melting = zero(dt_seconds)
    heating = zero(dt_seconds)

    if n_layers == 1
        updated_surface_temperature = (previous_surface_temperature + surface_rhs_term) /
                                      _safe_positive(one(previous_surface_temperature) + surface_diag_term)
        if updated_surface_temperature > c.T0
            needs_melt = true
            energy_to_melting = (c.T0 - previous_surface_temperature) * c.ci * surface_mass
            updated_surface_temperature = c.T0
            heating = energy_to_melting
        else
            heating = dt_seconds * (surface_flux_constant - surface_flux_linear * updated_surface_temperature)
        end
        resolved_surface_temperature = min(updated_surface_temperature, c.T0)
        _set_layer!(temperature, 1, idx, resolved_surface_temperature)
        _set_scalar!(Tsrf, idx, resolved_surface_temperature)
        melt_energy_available = _residual_melt_energy(
            surface_flux_constant,
            surface_flux_linear,
            resolved_surface_temperature,
            energy_to_melting,
            dt_seconds,
            needs_melt,
        )
        return _energy_flux_result(
            needs_melt,
            energy_to_melting,
            melt_energy_available,
            heating,
            surface_flux_constant,
            surface_flux_linear,
            latent_heat_linear_coefficient_eff,
            latent_heat_constant_term_eff,
        )
    end

    previous_layer_density = _get_layer(density, 1, idx)
    previous_layer_thickness = surface_mass / _safe_positive(previous_layer_density)
    previous_layer_conductivity = _snow_thermal_conductivity(previous_layer_density, c.Ki)
    _set_layer!(rhs, 1, idx, previous_surface_temperature + surface_rhs_term)

    @inbounds for layer_index in 2:n_layers
        layer_density = _get_layer(density, layer_index, idx)
        layer_thickness = _get_layer(mass, layer_index, idx) / _safe_positive(layer_density)
        layer_conductivity = _snow_thermal_conductivity(layer_density, c.Ki)
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
        _set_layer!(rhs, layer_index, idx, _get_layer(temperature, layer_index, idx))
        previous_layer_thickness = layer_thickness
        previous_layer_conductivity = layer_conductivity
    end

    β_scale = -oftype(dt_seconds, 2.0) * dt_seconds / c.ci

    β1 = β_scale / _safe_positive(_get_layer(mass, 1, idx))
    _set_layer!(upper, 1, idx, β1 * _get_layer(interface_terms, 1, idx))
    _set_layer!(diag, 1, idx, one(dt_seconds) - _get_layer(upper, 1, idx) + surface_diag_term)

    βn = β_scale / _safe_positive(_get_layer(mass, n_layers, idx))
    _set_layer!(lower, n_layers - 1, idx, βn * _get_layer(interface_terms, n_layers - 1, idx))
    _set_layer!(diag, n_layers, idx, one(dt_seconds) - _get_layer(lower, n_layers - 1, idx))

    @inbounds for layer_index in 2:(n_layers - 1)
        βi = β_scale / _safe_positive(_get_layer(mass, layer_index, idx))
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

    if _get_layer(resolved_temperature, 1, idx) > c.T0
        needs_melt = true
        energy_to_melting = (c.T0 - previous_surface_temperature) * c.ci * surface_mass

        @inbounds for layer_index in 1:n_layers
            _set_layer!(rhs, layer_index, idx, _get_layer(temperature, layer_index, idx))
        end
        _set_layer!(rhs, 1, idx, c.T0)
        _copy_column!(solver_diag, diag, idx, n_layers)
        _set_layer!(solver_diag, 1, idx, _get_layer(solver_diag, 1, idx) - surface_diag_term)

        resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, solver_diag, upper, rhs, idx, n_layers)
        energy_to_melting += (c.T0 - _get_layer(resolved_temperature, 1, idx)) * c.ci * surface_mass
        _set_layer!(resolved_temperature, 1, idx, c.T0)
        _clamp_to_melt!(resolved_temperature, idx, c.T0, n_layers)
        heating = energy_to_melting
    else
        _clamp_to_melt!(resolved_temperature, idx, c.T0, n_layers)
        heating = dt_seconds * (surface_flux_constant - surface_flux_linear * _get_layer(resolved_temperature, 1, idx))
    end

    _copy_column!(temperature, resolved_temperature, idx, n_layers)
    resolved_surface_temperature = _get_layer(resolved_temperature, 1, idx)
    _set_scalar!(Tsrf, idx, resolved_surface_temperature)

    melt_energy_available = _residual_melt_energy(
        surface_flux_constant,
        surface_flux_linear,
        resolved_surface_temperature,
        energy_to_melting,
        dt_seconds,
        needs_melt,
    )

    return _energy_flux_result(
        needs_melt,
        energy_to_melting,
        melt_energy_available,
        heating,
        surface_flux_constant,
        surface_flux_linear,
        latent_heat_linear_coefficient_eff,
        latent_heat_constant_term_eff,
    )
end

"""
    go_energy_flux!(domain, idx, air_temperature, shortwave_down, latent_heat_linear_coefficient, latent_heat_constant_term, dt_seconds; scratch, q_sw_net=nothing, q_lw_down=nothing, q_sh=nothing, q_lh=nothing)

Public wrapper for the column energy-flux solve on `domain`. Reuses `scratch`
for temporary arrays and returns the same diagnostic named tuple as
`_go_energy_flux_resolved!`.
"""
function go_energy_flux!(
    domain::AbstractSnowpackDomain,
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
    return _go_energy_flux_resolved!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.Tsrf,
        domain.albedo_dynamic,
        idx,
        domain.c,
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
    )
end
