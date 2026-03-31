#=
Energy-flux temperature solver for array-backed snowpack states.
=#

@inline _safe_positive(x) = x > EPS_TINY ? x : oftype(x, EPS_TINY)

@inline function _clamp_to_melt!(
    temperature_profile::AbstractVector,
    melting_temperature,
)
    @inbounds for layer_index in eachindex(temperature_profile)
        if temperature_profile[layer_index] > melting_temperature
            temperature_profile[layer_index] = melting_temperature
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

@inline function shortwave_absorbed(
    shortwave_down;
    surface_albedo,
)
    return (one(shortwave_down) - clamp(surface_albedo, zero(surface_albedo), one(surface_albedo))) * shortwave_down
end

@inline interface_conductance(Kᵢ, Δzᵢ, Kⱼ, Δzⱼ) =
    (Kᵢ * Δzᵢ + Kⱼ * Δzⱼ) / _safe_positive((Δzᵢ + Δzⱼ)^2)

function _solve_tridiagonal_thomas_prefix!(
    lower_diagonal::AbstractVector,
    main_diagonal::AbstractVector,
    upper_diagonal::AbstractVector,
    right_hand_side::AbstractVector,
    n::Int,
)
    @assert n >= 1
    @assert length(main_diagonal) >= n
    @assert length(right_hand_side) >= n
    @assert length(lower_diagonal) >= n - 1
    @assert length(upper_diagonal) >= n - 1

    @inbounds for row_index in 2:n
        elimination_factor = lower_diagonal[row_index - 1] / main_diagonal[row_index - 1]
        main_diagonal[row_index] -= elimination_factor * upper_diagonal[row_index - 1]
        right_hand_side[row_index] -= elimination_factor * right_hand_side[row_index - 1]
    end

    right_hand_side[n] /= main_diagonal[n]
    @inbounds for row_index in (n - 1):-1:1
        right_hand_side[row_index] = (
            right_hand_side[row_index] - upper_diagonal[row_index] * right_hand_side[row_index + 1]
        ) / main_diagonal[row_index]
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

    @inbounds for layer_index in 1:n_layers
        previous_temperature[layer_index] = _get_layer(temperature, layer_index, idx)
        layer_thickness[layer_index] = _get_layer(mass, layer_index, idx) / _safe_positive(_get_layer(density, layer_index, idx))
    end

    surface_mass = _safe_positive(_get_layer(mass, 1, idx))
    surface_temperature_scale = dt_seconds / c.ci / surface_mass

    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        shortwave_absorbed(shortwave_down; surface_albedo=_get_scalar(albedo_dynamic, idx))
    longwave_flux_constant = use_q_lw_down ?
        (q_lw_down_value + c.σ * c.ϵ_snow * oftype(air_temperature, 3.0) * previous_temperature[1]^4) :
        (c.σ * (c.ϵ_air * air_temperature^4 + c.ϵ_snow * oftype(air_temperature, 3.0) * previous_temperature[1]^4))
    longwave_flux_linear = c.σ * c.ϵ_snow * oftype(air_temperature, 4.0) * previous_temperature[1]^3
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
        updated_surface_temperature = (previous_temperature[1] + surface_rhs_term) / _safe_positive(one(previous_temperature[1]) + surface_diag_term)
        if updated_surface_temperature > c.T0
            needs_melt = true
            energy_to_melting = (c.T0 - previous_temperature[1]) * c.ci * surface_mass
            updated_surface_temperature = c.T0
            heating = energy_to_melting
        else
            heating = dt_seconds * (surface_flux_constant - surface_flux_linear * updated_surface_temperature)
        end
        _set_layer!(temperature, 1, idx, min(updated_surface_temperature, c.T0))
        _set_scalar!(Tsrf, idx, _get_layer(temperature, 1, idx))
        melt_energy_available = _residual_melt_energy(
            surface_flux_constant,
            surface_flux_linear,
            _get_layer(temperature, 1, idx),
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

    @inbounds for layer_index in 1:n_layers
        thermal_conductivity[layer_index] = _snow_thermal_conductivity(
            _get_layer(density, layer_index, idx),
            c.Ki,
            resolved_diffusion_model,
        )
    end
    @inbounds for layer_index in 1:(n_layers - 1)
        interface_terms[layer_index] = interface_conductance(
            thermal_conductivity[layer_index],
            layer_thickness[layer_index],
            thermal_conductivity[layer_index + 1],
            layer_thickness[layer_index + 1],
        )
    end

    function assemble_system!(surface_diag, use_melt_rhs::Bool)
        fill!(lower, zero(eltype(lower)))
        fill!(diag, zero(eltype(diag)))
        fill!(upper, zero(eltype(upper)))

        β1 = -oftype(dt_seconds, 2.0) * dt_seconds / (_safe_positive(_get_layer(density, 1, idx)) * c.ci * _safe_positive(layer_thickness[1]))
        upper[1] = β1 * interface_terms[1]
        diag[1] = one(dt_seconds) - upper[1] + surface_diag

        βn = -oftype(dt_seconds, 2.0) * dt_seconds / (_safe_positive(_get_layer(density, n_layers, idx)) * c.ci * _safe_positive(layer_thickness[n_layers]))
        lower[n_layers - 1] = βn * interface_terms[n_layers - 1]
        diag[n_layers] = one(dt_seconds) - lower[n_layers - 1]

        for layer_index in 2:(n_layers - 1)
            βi = -oftype(dt_seconds, 2.0) * dt_seconds / (_safe_positive(_get_layer(density, layer_index, idx)) * c.ci * _safe_positive(layer_thickness[layer_index]))
            lower[layer_index - 1] = βi * interface_terms[layer_index - 1]
            upper[layer_index] = βi * interface_terms[layer_index]
            diag[layer_index] = one(dt_seconds) - lower[layer_index - 1] - upper[layer_index]
        end

        copyto!(rhs, previous_temperature)
        if use_melt_rhs
            rhs[1] = c.T0
        else
            rhs[1] += surface_rhs_term
        end
        return nothing
    end

    assemble_system!(surface_diag_term, false)
    resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, diag, upper, rhs, n_layers)

    if resolved_temperature[1] > c.T0
        needs_melt = true
        energy_to_melting = (c.T0 - previous_temperature[1]) * c.ci * surface_mass

        assemble_system!(zero(surface_diag_term), true)
        resolved_temperature = _solve_tridiagonal_thomas_prefix!(lower, diag, upper, rhs, n_layers)
        energy_to_melting += (c.T0 - resolved_temperature[1]) * c.ci * surface_mass
        resolved_temperature[1] = c.T0
        _clamp_to_melt!(resolved_temperature, c.T0)
        heating = energy_to_melting
    else
        _clamp_to_melt!(resolved_temperature, c.T0)
        heating = dt_seconds * (surface_flux_constant - surface_flux_linear * resolved_temperature[1])
    end

    @inbounds for layer_index in 1:n_layers
        _set_layer!(temperature, layer_index, idx, resolved_temperature[layer_index])
    end
    _set_scalar!(Tsrf, idx, _get_layer(temperature, 1, idx))

    melt_energy_available = _residual_melt_energy(
        surface_flux_constant,
        surface_flux_linear,
        _get_layer(temperature, 1, idx),
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
