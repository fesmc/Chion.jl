#=
Energy-flux temperature solver for array-backed snowpack states.
=#

@inline _safe_positive(x) = x > EPS_TINY ? x : oftype(x, EPS_TINY)

@inline function _clamp_to_melt!(
    temperature_profile,
    idx::Int,
    melting_temperature,
    n::Int,
)
    @inbounds for layer_index in 1:n
        if _get_layer(temperature_profile, layer_index, idx) > melting_temperature
            _set_layer!(temperature_profile, layer_index, idx, melting_temperature)
        end
    end
    return temperature_profile
end

function _snow_thermal_conductivity(ρ, Kᵢ, diffusion_model::Int)
    if diffusion_model == 1
        return Kᵢ * (ρ / oftype(ρ, 1000.0))^oftype(ρ, 1.88)
    elseif diffusion_model == 2
        if ρ > oftype(ρ, 156.0)
            return oftype(ρ, 0.138) - oftype(ρ, 1.01e-3) * ρ + oftype(ρ, 3.233e-6) * ρ^2
        else
            return oftype(ρ, 0.023) + oftype(ρ, 0.234e-3) * ρ
        end
    else
        return oftype(ρ, 2.1e-2) + oftype(ρ, 4.2e-4) * ρ + oftype(ρ, 2.2e-9) * ρ^3
    end
end

@inline _snow_thermal_conductivity_model1(ρ, Kᵢ) =
    Kᵢ * (ρ * oftype(ρ, 1.0e-3))^oftype(ρ, 1.88)

@inline function _snow_thermal_conductivity_model2(ρ)
    if ρ > oftype(ρ, 156.0)
        return oftype(ρ, 0.138) - oftype(ρ, 1.01e-3) * ρ + oftype(ρ, 3.233e-6) * ρ^2
    end
    return oftype(ρ, 0.023) + oftype(ρ, 0.234e-3) * ρ
end

@inline _snow_thermal_conductivity_model3(ρ) =
    oftype(ρ, 2.1e-2) + oftype(ρ, 4.2e-4) * ρ + oftype(ρ, 2.2e-9) * ρ^3

@inline function shortwave_absorbed(
    shortwave_down;
    surface_albedo,
)
    return (one(shortwave_down) - clamp(surface_albedo, zero(surface_albedo), one(surface_albedo))) * shortwave_down
end

@inline interface_conductance(Kᵢ, Δzᵢ, Kⱼ, Δzⱼ) =
    (Kᵢ * Δzᵢ + Kⱼ * Δzⱼ) / _safe_positive((Δzᵢ + Δzⱼ)^2)

function _solve_tridiagonal_thomas_prefix!(
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

@inline function _surface_has_snow(N_storage, mass, idx::Int)
    return _n_active(N_storage, idx) > 0 && _get_layer(mass, 1, idx) > EPS_EMPTY_LAYER
end

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
    resolved_diffusion_model::Int,
    use_q_sw_net::Bool,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
)
    n_layers = _n_active(N_storage, idx)
    if n_layers <= 0 || _get_layer(mass, 1, idx) <= zero(eltype(mass))
        return _energy_flux_result(
            needs_melt=false,
            energy_to_melting=zero(dt_seconds),
            melt_energy_available=zero(dt_seconds),
            heating=zero(dt_seconds),
            surface_flux_constant=zero(dt_seconds),
            surface_flux_linear=zero(dt_seconds),
            latent_heat_linear_coefficient=latent_heat_linear_coefficient_eff,
            latent_heat_constant_term=latent_heat_constant_term_eff,
        )
    end

    lower = scratch.lower
    diag = scratch.diag
    upper = scratch.upper
    rhs = scratch.rhs
    interface_terms = scratch.interface_conductance
    previous_temperature = scratch.previous_temperature
    layer_thickness = scratch.layer_thickness
    thermal_conductivity = scratch.thermal_conductivity

    if resolved_diffusion_model == 1
        @inbounds for layer_index in 1:n_layers
            layer_density = _get_layer(density, layer_index, idx)
            layer_mass = _get_layer(mass, layer_index, idx)
            layer_temperature = _get_layer(temperature, layer_index, idx)
            _set_layer!(previous_temperature, layer_index, idx, layer_temperature)
            _set_layer!(rhs, layer_index, idx, layer_temperature)
            _set_layer!(layer_thickness, layer_index, idx, layer_mass / _safe_positive(layer_density))
            _set_layer!(thermal_conductivity, layer_index, idx, _snow_thermal_conductivity_model1(layer_density, c.Ki))
        end
    elseif resolved_diffusion_model == 2
        @inbounds for layer_index in 1:n_layers
            layer_density = _get_layer(density, layer_index, idx)
            layer_mass = _get_layer(mass, layer_index, idx)
            layer_temperature = _get_layer(temperature, layer_index, idx)
            _set_layer!(previous_temperature, layer_index, idx, layer_temperature)
            _set_layer!(rhs, layer_index, idx, layer_temperature)
            _set_layer!(layer_thickness, layer_index, idx, layer_mass / _safe_positive(layer_density))
            _set_layer!(thermal_conductivity, layer_index, idx, _snow_thermal_conductivity_model2(layer_density))
        end
    else
        @inbounds for layer_index in 1:n_layers
            layer_density = _get_layer(density, layer_index, idx)
            layer_mass = _get_layer(mass, layer_index, idx)
            layer_temperature = _get_layer(temperature, layer_index, idx)
            _set_layer!(previous_temperature, layer_index, idx, layer_temperature)
            _set_layer!(rhs, layer_index, idx, layer_temperature)
            _set_layer!(layer_thickness, layer_index, idx, layer_mass / _safe_positive(layer_density))
            _set_layer!(thermal_conductivity, layer_index, idx, _snow_thermal_conductivity_model3(layer_density))
        end
    end

    surface_mass = _safe_positive(_get_layer(mass, 1, idx))
    previous_surface_temperature = _get_layer(previous_temperature, 1, idx)
    surface_temperature_scale = dt_seconds / c.ci / surface_mass
    surface_temperature_sq = previous_surface_temperature * previous_surface_temperature
    surface_temperature_cube = surface_temperature_sq * previous_surface_temperature
    surface_temperature_fourth = surface_temperature_sq * surface_temperature_sq

    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        shortwave_absorbed(shortwave_down; surface_albedo=_get_scalar(albedo_dynamic, idx))
    longwave_flux_constant = use_q_lw_down ?
        (q_lw_down_value + c.σ * c.ϵ_snow * oftype(air_temperature, 3.0) * surface_temperature_fourth) :
        (c.σ * (c.ϵ_air * air_temperature^4 + c.ϵ_snow * oftype(air_temperature, 3.0) * surface_temperature_fourth))
    longwave_flux_linear = c.σ * c.ϵ_snow * oftype(air_temperature, 4.0) * surface_temperature_cube
    sensible_heat_flux_constant = use_q_sh ? q_sh_value : air_temperature * c.D_sh
    sensible_heat_flux_linear = use_q_sh ? zero(dt_seconds) : c.D_sh
    latent_heat_flux_constant = use_q_lh ? q_lh_value : latent_heat_constant_term_eff
    latent_heat_flux_linear = use_q_lh ? zero(dt_seconds) : latent_heat_linear_coefficient_eff

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
            dt_seconds;
            needs_melt=needs_melt,
        )
        return _energy_flux_result(
            needs_melt=needs_melt,
            energy_to_melting=energy_to_melting,
            melt_energy_available=melt_energy_available,
            heating=heating,
            surface_flux_constant=surface_flux_constant,
            surface_flux_linear=surface_flux_linear,
            latent_heat_linear_coefficient=latent_heat_linear_coefficient_eff,
            latent_heat_constant_term=latent_heat_constant_term_eff,
        )
    end

    @inbounds for layer_index in 1:(n_layers - 1)
        interface_term = interface_conductance(
            _get_layer(thermal_conductivity, layer_index, idx),
            _get_layer(layer_thickness, layer_index, idx),
            _get_layer(thermal_conductivity, layer_index + 1, idx),
            _get_layer(layer_thickness, layer_index + 1, idx),
        )
        _set_layer!(interface_terms, layer_index, idx, interface_term)
    end

    function assemble_system!(surface_diag, use_melt_rhs::Bool)
        β_scale = -oftype(dt_seconds, 2.0) * dt_seconds / c.ci

        β1 = β_scale / _safe_positive(_get_layer(mass, 1, idx))
        _set_layer!(upper, 1, idx, β1 * _get_layer(interface_terms, 1, idx))
        _set_layer!(diag, 1, idx, one(dt_seconds) - _get_layer(upper, 1, idx) + surface_diag)

        βn = β_scale / _safe_positive(_get_layer(mass, n_layers, idx))
        _set_layer!(lower, n_layers - 1, idx, βn * _get_layer(interface_terms, n_layers - 1, idx))
        _set_layer!(diag, n_layers, idx, one(dt_seconds) - _get_layer(lower, n_layers - 1, idx))

        for layer_index in 2:(n_layers - 1)
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

        if use_melt_rhs
            @inbounds for layer_index in 1:n_layers
                _set_layer!(rhs, layer_index, idx, _get_layer(previous_temperature, layer_index, idx))
            end
            _set_layer!(rhs, 1, idx, c.T0)
        else
            _set_layer!(rhs, 1, idx, previous_surface_temperature + surface_rhs_term)
        end
        return nothing
    end

    assemble_system!(surface_diag_term, false)
    resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, diag, upper, rhs, idx, n_layers)

    if _get_layer(resolved_temperature, 1, idx) > c.T0
        needs_melt = true
        energy_to_melting = (c.T0 - previous_surface_temperature) * c.ci * surface_mass

        assemble_system!(zero(surface_diag_term), true)
        resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, diag, upper, rhs, idx, n_layers)
        energy_to_melting += (c.T0 - _get_layer(resolved_temperature, 1, idx)) * c.ci * surface_mass
        _set_layer!(resolved_temperature, 1, idx, c.T0)
        _clamp_to_melt!(resolved_temperature, idx, c.T0, n_layers)
        heating = energy_to_melting
    else
        _clamp_to_melt!(resolved_temperature, idx, c.T0, n_layers)
        heating = dt_seconds * (surface_flux_constant - surface_flux_linear * _get_layer(resolved_temperature, 1, idx))
    end

    @inbounds for layer_index in 1:n_layers
        _set_layer!(temperature, layer_index, idx, _get_layer(resolved_temperature, layer_index, idx))
    end
    resolved_surface_temperature = _get_layer(resolved_temperature, 1, idx)
    _set_scalar!(Tsrf, idx, resolved_surface_temperature)

    melt_energy_available = _residual_melt_energy(
        surface_flux_constant,
        surface_flux_linear,
        resolved_surface_temperature,
        energy_to_melting,
        dt_seconds;
        needs_melt=needs_melt,
    )

    return _energy_flux_result(
        needs_melt=needs_melt,
        energy_to_melting=energy_to_melting,
        melt_energy_available=melt_energy_available,
        heating=heating,
        surface_flux_constant=surface_flux_constant,
        surface_flux_linear=surface_flux_linear,
        latent_heat_linear_coefficient=latent_heat_linear_coefficient_eff,
        latent_heat_constant_term=latent_heat_constant_term_eff,
    )
end

function go_energy_flux!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    air_temperature,
    shortwave_down,
    latent_heat_linear_coefficient,
    latent_heat_constant_term,
    dt_seconds;
    scratch=EnergyWorkspace(domain),
    diffusion_model::Int=2,
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
        diffusion_model,
        !isnothing(q_sw_net),
        isnothing(q_sw_net) ? zero(dt_seconds) : q_sw_net,
        !isnothing(q_lw_down),
        isnothing(q_lw_down) ? zero(dt_seconds) : q_lw_down,
        !isnothing(q_sh),
        isnothing(q_sh) ? zero(dt_seconds) : q_sh,
        !isnothing(q_lh),
        isnothing(q_lh) ? zero(dt_seconds) : q_lh,
    )
end
