#=
Energy-flux temperature solver for `SnowpackColumn`.

Port of BESSI `go_energy_flux_new` to the Chion column nomenclature.
=#
using LinearAlgebra

@inline _safe_positive(x::Float64) = x > EPS_TINY ? x : EPS_TINY

@inline function _clamp_to_melt!(
    temperature_profile::AbstractVector{Float64},
    melting_temperature::Float64,
)
    @inbounds for layer_index in eachindex(temperature_profile)
        if temperature_profile[layer_index] > melting_temperature
            temperature_profile[layer_index] = melting_temperature
        end
    end
    return temperature_profile
end

function _snow_thermal_conductivity(ρ::Float64, Kᵢ::Float64, diffusion_model::Int)
    if diffusion_model == 1
        # Yen (1981)
        return Kᵢ * (ρ / 1000.0)^1.88
    elseif diffusion_model == 2
        # Sturm (1997), piecewise form
        if ρ > 156.0
            return 0.138 - 1.01e-3 * ρ + 3.233e-6 * ρ^2
        else
            return 0.023 + 0.234e-3 * ρ
        end
    else
        # Van Dusen (1929)
        return 2.1e-2 + 4.2e-4 * ρ + 2.2e-9 * ρ^3
    end
end

@inline function shortwave_absorbed(
    shortwave_down::Float64,
    surface_temperature::Float64,
    melting_temperature::Float64;
    snow_cover::Float64=1.0,
    alpha_dry::Float64=0.8,
    alpha_wet::Float64=0.6,
    alpha_ice::Float64=0.35,
)
    # Keep the legacy BESSI-style path for now: partial snow cover does not yet
    # blend in the ice albedo, even though the arguments are already available.
    snow_albedo = (surface_temperature >= melting_temperature) ? alpha_wet : alpha_dry
    return (1.0 - snow_albedo) * shortwave_down
end

@inline interface_conductance(Kᵢ, Δzᵢ, Kⱼ, Δzⱼ) =
    (Kᵢ * Δzᵢ + Kⱼ * Δzⱼ) / _safe_positive((Δzᵢ + Δzⱼ)^2)

@inline function _normalize_tridiagonal_solver(solver::Symbol)
    solver in (:linear_algebra, :thomas) ||
        error("Unsupported tridiagonal solver '$solver'. Use :linear_algebra or :thomas.")
    return solver
end

function _solve_tridiagonal_thomas!(
    lower_diagonal::AbstractVector{Float64},
    main_diagonal::AbstractVector{Float64},
    upper_diagonal::AbstractVector{Float64},
    right_hand_side::AbstractVector{Float64},
)
    n = length(main_diagonal)
    @assert length(lower_diagonal) == n - 1
    @assert length(upper_diagonal) == n - 1
    @assert length(right_hand_side) == n

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

function _solve_tridiagonal_system(
    lower_diagonal::AbstractVector{Float64},
    main_diagonal::AbstractVector{Float64},
    upper_diagonal::AbstractVector{Float64},
    right_hand_side::AbstractVector{Float64};
    solver::Symbol=:linear_algebra,
)
    normalized_solver = _normalize_tridiagonal_solver(solver)
    if normalized_solver == :linear_algebra
        system_matrix = Tridiagonal(lower_diagonal, main_diagonal, upper_diagonal)
        return system_matrix \ right_hand_side
    end

    return _solve_tridiagonal_thomas!(
        copy(lower_diagonal),
        copy(main_diagonal),
        copy(upper_diagonal),
        copy(right_hand_side),
    )
end

@inline function _energy_flux_result(
    ;
    needs_melt::Bool,
    energy_to_melting::Float64,
    heating::Float64,
    surface_flux_constant::Float64,
    surface_flux_linear::Float64,
    latent_heat_linear_coefficient::Float64,
    latent_heat_constant_term::Float64,
)
    return (
        needs_melt = needs_melt,
        china_syndrome = needs_melt,
        energy_to_melting = energy_to_melting,
        Q_heat = energy_to_melting,
        heating = heating,
        surface_flux_constant = surface_flux_constant,
        F_const = surface_flux_constant,
        surface_flux_linear = surface_flux_linear,
        F_lin = surface_flux_linear,
        latent_heat_linear_coefficient = latent_heat_linear_coefficient,
        H_lh = latent_heat_linear_coefficient,
        latent_heat_constant_term = latent_heat_constant_term,
        K_lh = latent_heat_constant_term,
    )
end

@inline function _diagnose_latent_heat_flux_coefficients(
    column::SnowpackColumn,
    air_temperature::Float64,
    snowfall_rate::Float64,
    rainfall_rate::Float64,
)
    # Match Fortran logic from accumulation branch:
    # - snowfall: H_lh = accum*c_i, K_lh = accum*c_i*T_air
    # - rainfall on existing snow: H_lh = 0, K_lh = rainman*c_w*(T_air-T0)
    if snowfall_rate > 0.0
        latent_heat_linear_coefficient = snowfall_rate * column.c.ci
        latent_heat_constant_term = snowfall_rate * column.c.ci * air_temperature
    elseif column.N > 0 && column.mass[1] > 0.0 && rainfall_rate > 0.0
        latent_heat_linear_coefficient = 0.0
        latent_heat_constant_term = rainfall_rate * column.c.cw * (air_temperature - column.c.T0)
    else
        latent_heat_linear_coefficient = 0.0
        latent_heat_constant_term = 0.0
    end
    return latent_heat_linear_coefficient, latent_heat_constant_term
end

@inline function bare_ice_ablation_mass(
    column::SnowpackColumn,
    air_temperature::Float64,
    rainfall_rate::Float64,
    dt_seconds::Float64;
    shortwave_down::Union{Nothing, Float64}=nothing,
    q_sw_net::Union{Nothing, Float64}=nothing,
    q_lw_down::Union{Nothing, Float64}=nothing,
    q_sh::Union{Nothing, Float64}=nothing,
    q_lh::Union{Nothing, Float64}=nothing,
)
    dt_seconds <= 0.0 && return 0.0

    c = column.c
    absorbed_shortwave = isnothing(q_sw_net) ?
        (isnothing(shortwave_down) ? 400.0 : max(shortwave_down, 0.0)) * (1.0 - c.alpha_ice) :
        q_sw_net
    longwave_flux = isnothing(q_lw_down) ?
        c.σ * (c.ϵ_air * air_temperature^4 - c.ϵ_snow * c.T0^4) :
        q_lw_down - c.σ * c.ϵ_snow * c.T0^4
    sensible_heat_flux = isnothing(q_sh) ? c.D_sh * (air_temperature - c.T0) : q_sh
    latent_heat_flux = isnothing(q_lh) ? 0.0 : q_lh
    rain_heat_flux = rainfall_rate * c.cw * (air_temperature - c.T0)
    net_surface_flux = absorbed_shortwave + longwave_flux + sensible_heat_flux + latent_heat_flux + rain_heat_flux
    return max(net_surface_flux, 0.0) * dt_seconds / c.Lm
end

"""
    go_energy_flux!(
        column::SnowpackColumn,
        air_temperature::Float64,
        shortwave_down::Float64,
        latent_heat_linear_coefficient::Union{Nothing, Float64},
        latent_heat_constant_term::Union{Nothing, Float64},
        dt_seconds::Float64;
        snowfall_rate::Float64=0.0,
        rainfall_rate::Float64=0.0,
        diffusion_model::Int = 2,
        tridiagonal_solver::Symbol = :linear_algebra,
        q_sw_net::Union{Nothing, Float64}=nothing,
        q_lw_down::Union{Nothing, Float64}=nothing,
        q_sh::Union{Nothing, Float64}=nothing,
        q_lh::Union{Nothing, Float64}=nothing,
    ) -> NamedTuple

Update snow temperatures with an implicit conductive solve and surface flux terms.

Inputs:
- `shortwave_down`: incoming solar radiation at the bottom of the atmosphere [W m^-2]
- `latent_heat_linear_coefficient`: optional latent-heat linear coefficient [W m^-2 K^-1]
- `latent_heat_constant_term`: optional latent-heat constant term [W m^-2]
- `snowfall_rate`: snowfall rate [kg m^-2 s^-1], used to diagnose latent-heat coefficients when they are `nothing`
- `rainfall_rate`: rainfall rate [kg m^-2 s^-1], used to diagnose latent-heat coefficients when they are `nothing`
- `q_sw_net`: optional observed net shortwave (down-up) [W m^-2]
- `q_lw_down`: optional observed downward longwave [W m^-2]
- `q_sh`: optional observed sensible heat flux, downward positive [W m^-2]
- `q_lh`: optional observed latent heat flux, downward positive [W m^-2]
- `tridiagonal_solver`: one of `:linear_algebra` or `:thomas`

Legacy keyword aliases `P_snow`, `P_rain`, and `diff_model` are still accepted.

Returns:
- `needs_melt` / `china_syndrome`: surface reached melt point and melt routine should run
- `energy_to_melting` / `Q_heat`: energy used to bring/cap surface at melting point [J m^-2]
- `heating`: diagnosed net heating term [J m^-2]
"""
function go_energy_flux!(
    column::SnowpackColumn,
    air_temperature::Float64,
    shortwave_down::Float64,
    latent_heat_linear_coefficient::Union{Nothing, Float64},
    latent_heat_constant_term::Union{Nothing, Float64},
    dt_seconds::Float64;
    snowfall_rate::Union{Nothing, Float64}=nothing,
    rainfall_rate::Union{Nothing, Float64}=nothing,
    P_snow::Union{Nothing, Float64}=nothing,
    P_rain::Union{Nothing, Float64}=nothing,
    diffusion_model::Union{Nothing, Int}=nothing,
    diff_model::Union{Nothing, Int}=nothing,
    tridiagonal_solver::Symbol=:linear_algebra,
    q_sw_net::Union{Nothing, Float64}=nothing,
    q_lw_down::Union{Nothing, Float64}=nothing,
    q_sh::Union{Nothing, Float64}=nothing,
    q_lh::Union{Nothing, Float64}=nothing,
)
    resolved_snowfall_rate = _resolve_keyword_alias(snowfall_rate, P_snow, "snowfall_rate", "P_snow")
    resolved_rainfall_rate = _resolve_keyword_alias(rainfall_rate, P_rain, "rainfall_rate", "P_rain")
    resolved_diffusion_model = _resolve_keyword_alias(
        diffusion_model,
        diff_model,
        "diffusion_model",
        "diff_model",
    )
    resolved_snowfall_rate = isnothing(resolved_snowfall_rate) ? 0.0 : resolved_snowfall_rate
    resolved_rainfall_rate = isnothing(resolved_rainfall_rate) ? 0.0 : resolved_rainfall_rate
    resolved_diffusion_model = isnothing(resolved_diffusion_model) ? 2 : resolved_diffusion_model
    resolved_tridiagonal_solver = _normalize_tridiagonal_solver(tridiagonal_solver)

    latent_heat_linear_coefficient_eff, latent_heat_constant_term_eff =
        if isnothing(latent_heat_linear_coefficient) || isnothing(latent_heat_constant_term)
            _diagnose_latent_heat_flux_coefficients(
                column,
                air_temperature,
                resolved_snowfall_rate,
                resolved_rainfall_rate,
            )
    else
        latent_heat_linear_coefficient, latent_heat_constant_term
    end

    update_snow_cover!(column)

    n_layers = column.N
    if n_layers <= 0 || column.mass[1] <= 0.0
        return _energy_flux_result(
            needs_melt=false,
            energy_to_melting=0.0,
            heating = 0.0,
            surface_flux_constant = 0.0,
            surface_flux_linear = 0.0,
            latent_heat_linear_coefficient = latent_heat_linear_coefficient_eff,
            latent_heat_constant_term = latent_heat_constant_term_eff,
        )
    end

    melting_temperature = column.c.T0
    ice_heat_capacity = column.c.ci
    sensible_heat_coefficient = column.c.D_sh
    stefan_boltzmann = column.c.σ
    air_emissivity = column.c.ϵ_air
    snow_emissivity = column.c.ϵ_snow
    ice_conductivity = column.c.Ki

    a_upper = zeros(Float64, n_layers)
    a_lower = zeros(Float64, n_layers)
    a_diag = zeros(Float64, n_layers)

    @views temperature_profile = column.temperature[1:n_layers]
    previous_temperature_profile = copy(temperature_profile)
    @views snow_density = column.density[1:n_layers]
    layer_thickness = zeros(Float64, n_layers)
    for layer_index in 1:n_layers
        layer_thickness[layer_index] = column.mass[layer_index] / _safe_positive(snow_density[layer_index])
    end

    surface_mass = _safe_positive(column.mass[1])
    surface_temperature_scale = dt_seconds / ice_heat_capacity / surface_mass

    absorbed_shortwave = isnothing(q_sw_net) ? shortwave_absorbed(
        shortwave_down,
        temperature_profile[1],
        melting_temperature;
        snow_cover=column.snow_cover,
        alpha_dry=column.c.alpha_dry,
        alpha_wet=column.c.alpha_wet,
        alpha_ice=column.c.alpha_ice,
    ) : q_sw_net

    longwave_flux_constant = isnothing(q_lw_down) ?
        (stefan_boltzmann * (
            air_emissivity * air_temperature^4 +
            snow_emissivity * 3.0 * temperature_profile[1]^4
        )) :
        (q_lw_down + stefan_boltzmann * snow_emissivity * 3.0 * temperature_profile[1]^4)
    longwave_flux_linear = stefan_boltzmann * snow_emissivity * 4.0 * temperature_profile[1]^3

    sensible_heat_flux_constant = isnothing(q_sh) ? (air_temperature * sensible_heat_coefficient) : q_sh
    sensible_heat_flux_linear = isnothing(q_sh) ? sensible_heat_coefficient : 0.0

    latent_heat_flux_constant = isnothing(q_lh) ? latent_heat_constant_term_eff : q_lh
    latent_heat_flux_linear = isnothing(q_lh) ? latent_heat_linear_coefficient_eff : 0.0

    # Net surface flux linearized as: F_const - F_lin * Ts
    surface_flux_constant = sensible_heat_flux_constant +
                            longwave_flux_constant +
                            absorbed_shortwave +
                            latent_heat_flux_constant
    surface_flux_linear = sensible_heat_flux_linear +
                          longwave_flux_linear +
                          latent_heat_flux_linear
    surface_rhs_term = surface_temperature_scale * surface_flux_constant
    surface_diag_term = surface_temperature_scale * surface_flux_linear

    needs_melt = false
    energy_to_melting = 0.0
    heating = 0.0

    if n_layers == 1
        updated_surface_temperature = (
            temperature_profile[1] + surface_rhs_term
        ) / _safe_positive(1.0 + surface_diag_term)
        if updated_surface_temperature > melting_temperature
            needs_melt = true
            energy_to_melting = (
                melting_temperature - previous_temperature_profile[1]
            ) * ice_heat_capacity * surface_mass
            updated_surface_temperature = melting_temperature
            heating = energy_to_melting
        else
            heating = dt_seconds * (
                surface_flux_constant - surface_flux_linear * updated_surface_temperature
            )
        end

        column.temperature[1] = min(updated_surface_temperature, melting_temperature)
        column.Tsrf = column.temperature[1]
        return _energy_flux_result(
            needs_melt=needs_melt,
            energy_to_melting=energy_to_melting,
            heating = heating,
            surface_flux_constant = surface_flux_constant,
            surface_flux_linear = surface_flux_linear,
            latent_heat_linear_coefficient = latent_heat_linear_coefficient_eff,
            latent_heat_constant_term = latent_heat_constant_term_eff,
        )
    end

    thermal_conductivity = zeros(Float64, n_layers)
    for layer_index in 1:n_layers
        thermal_conductivity[layer_index] = _snow_thermal_conductivity(
            snow_density[layer_index],
            ice_conductivity,
            resolved_diffusion_model,
        )
    end

    @inbounds begin
        # i = 1 boundary
        k₁₂ = interface_conductance(
            thermal_conductivity[1],
            layer_thickness[1],
            thermal_conductivity[2],
            layer_thickness[2],
        )
        β₁ = -2.0 * dt_seconds / (
            _safe_positive(snow_density[1]) *
            ice_heat_capacity *
            _safe_positive(layer_thickness[1])
        )
        a_upper[1] = β₁ * k₁₂
        a_diag[1] = 1.0 - a_upper[1]

        # i = n_layers boundary
        kₙ = interface_conductance(
            thermal_conductivity[n_layers],
            layer_thickness[n_layers],
            thermal_conductivity[n_layers - 1],
            layer_thickness[n_layers - 1],
        )
        βₙ = -2.0 * dt_seconds / (
            _safe_positive(snow_density[n_layers]) *
            ice_heat_capacity *
            _safe_positive(layer_thickness[n_layers])
        )
        a_lower[n_layers] = βₙ * kₙ
        a_diag[n_layers] = 1.0 - a_lower[n_layers]

        # interior
        for layer_index in 2:n_layers-1
            lower_interface_conductance = interface_conductance(
                thermal_conductivity[layer_index],
                layer_thickness[layer_index],
                thermal_conductivity[layer_index - 1],
                layer_thickness[layer_index - 1],
            )
            upper_interface_conductance = interface_conductance(
                thermal_conductivity[layer_index],
                layer_thickness[layer_index],
                thermal_conductivity[layer_index + 1],
                layer_thickness[layer_index + 1],
            )
            βᵢ = -2.0 * dt_seconds / (
                _safe_positive(snow_density[layer_index]) *
                ice_heat_capacity *
                _safe_positive(layer_thickness[layer_index])
            )
            a_lower[layer_index] = βᵢ * lower_interface_conductance
            a_upper[layer_index] = βᵢ * upper_interface_conductance
            a_diag[layer_index] = 1.0 - a_lower[layer_index] - a_upper[layer_index]
        end
    end

    a_diag₁_backup = a_diag[1]
    a_diag[1] += surface_diag_term

    right_hand_side = copy(temperature_profile)
    right_hand_side[1] += surface_rhs_term
    updated_temperature = _solve_tridiagonal_system(
        a_lower[2:end],
        a_diag,
        a_upper[1:(end - 1)],
        right_hand_side;
        solver=resolved_tridiagonal_solver,
    )

    if updated_temperature[1] > melting_temperature
        needs_melt = true
        energy_to_melting = (
            melting_temperature - previous_temperature_profile[1]
        ) * ice_heat_capacity * surface_mass

        right_hand_side_melt = copy(previous_temperature_profile)
        right_hand_side_melt[1] = melting_temperature
        a_diag[1] = a_diag₁_backup
        updated_temperature = _solve_tridiagonal_system(
            a_lower[2:end],
            a_diag,
            a_upper[1:(end - 1)],
            right_hand_side_melt;
            solver=resolved_tridiagonal_solver,
        )

        energy_to_melting += (
            melting_temperature - updated_temperature[1]
        ) * ice_heat_capacity * surface_mass
        updated_temperature[1] = melting_temperature
        _clamp_to_melt!(updated_temperature, melting_temperature)
        heating = energy_to_melting
    else
        _clamp_to_melt!(updated_temperature, melting_temperature)
        heating = dt_seconds * (
            surface_flux_constant - surface_flux_linear * updated_temperature[1]
        )
    end

    column.temperature[1:n_layers] .= updated_temperature
    column.Tsrf = column.temperature[1]

    return _energy_flux_result(
        needs_melt=needs_melt,
        energy_to_melting=energy_to_melting,
        heating = heating,
        surface_flux_constant = surface_flux_constant,
        surface_flux_linear = surface_flux_linear,
        latent_heat_linear_coefficient = latent_heat_linear_coefficient_eff,
        latent_heat_constant_term = latent_heat_constant_term_eff,
    )
end
