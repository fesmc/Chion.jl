#=
Energy-flux temperature solver for `SnowpackColumn`.

Port of BESSI `go_energy_flux_new` to the Chion column nomenclature.
=#
using LinearAlgebra

const _TINY = 1.0e-12
@inline _safe_positive(x::Float64) = x > _TINY ? x : _TINY

@inline function _clamp_to_melt!(T::Vector{Float64}, Tₘ::Float64)
    @inbounds for i in eachindex(T)
        if T[i] > Tₘ
            T[i] = Tₘ
        end
    end
    return T
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
    S_boa::Float64,
    Tₛ::Float64,
    Tₘ::Float64;
    alpha_dry::Float64=0.8,
    alpha_wet::Float64=0.6,
)
    # wet when at/near melting
    α = (Tₛ >= Tₘ) ? alpha_wet : alpha_dry
    return (1.0 - α) * S_boa
end

@inline interface_conductance(Kᵢ, Δzᵢ, Kⱼ, Δzⱼ) =
    (Kᵢ * Δzᵢ + Kⱼ * Δzⱼ) / _safe_positive((Δzᵢ + Δzⱼ)^2)

"""
    go_energy_flux!(
        column::SnowpackColumn,
        T₂m::Float64,
        S_boa::Float64,
        H_lh::Float64,
        K_lh::Float64,
        dt_sec::Float64;
        diff_model::Int = 2,
    ) -> NamedTuple

Update snow temperatures with an implicit conductive solve and surface flux terms.

Inputs:
- `S_boa`: incoming solar radiation at the bottom of the atmosphere [W m^-2]
- `H_lh`: latent-heat linear coefficient [W m^-2 K^-1]
- `K_lh`: latent-heat constant term [W m^-2]

Returns:
- `china_syndrome`: surface reached melt point and melt routine should run
- `Q_heat`: energy used to bring/cap surface at melting point [J m^-2]
- `heating`: diagnosed net heating term [J m^-2]
"""
function go_energy_flux!(
    column::SnowpackColumn,
    T₂m::Float64,
    S_boa::Float64,
    H_lh::Float64,
    K_lh::Float64,
    dt_sec::Float64;
    diff_model::Int = 2,
)
    n_layers = column.N
    if n_layers <= 0 || column.mass[1] <= 0.0
        return (china_syndrome = false, Q_heat = 0.0, heating = 0.0)
    end

    Tₘ = column.c.T0
    cᵢ = column.c.ci
    Dₛₕ = column.c.D_sh
    σ = column.c.σ
    ϵₐ = column.c.ϵ_air
    ϵₛ = column.c.ϵ_snow
    Kᵢ = column.c.Ki

    a_upper = zeros(Float64, n_layers)
    a_lower = zeros(Float64, n_layers)
    a_diag = zeros(Float64, n_layers)

    T = copy(column.temperature[1:n_layers])
    T_prev = copy(T)
    ρ = copy(column.density[1:n_layers])
    Δz = zeros(Float64, n_layers)
    for i in 1:n_layers
        Δz[i] = column.mass[i] / _safe_positive(ρ[i])
    end

    mₛ = _safe_positive(column.mass[1])
    λₛ = dt_sec / cᵢ / mₛ
    Q_sw = shortwave_absorbed(
        S_boa,
        T[1],
        Tₘ;
        alpha_dry=column.c.alpha_dry,
        alpha_wet=column.c.alpha_wet,
    )
    ### Taylor expansion of T
    F_const = (
        T₂m * Dₛₕ +
        σ * (ϵₐ * T₂m^4 + ϵₛ * 3.0 * T[1]^4) +
        Q_sw +
        K_lh
    )
    F_lin = Dₛₕ + σ * ϵₛ * 4.0 * T[1]^3 + H_lh
    surface_rhs_term = λₛ * F_const
    surface_diag_term = λₛ * F_lin

    china_syndrome = false
    Q_heat = 0.0
    heating = 0.0

    if n_layers == 1
        T_new = (T[1] + surface_rhs_term) / _safe_positive(1.0 + surface_diag_term)
        if T_new > Tₘ
            china_syndrome = true
            Q_heat = max((Tₘ - T_prev[1]) * cᵢ * mₛ, 0.0)
            T_new = Tₘ
            heating = Q_heat
        else
            heating = dt_sec * (F_const - F_lin * T_new)
        end

        column.temperature[1] = min(T_new, Tₘ)
        column.Tsrf = column.temperature[1]
        return (china_syndrome = china_syndrome, Q_heat = Q_heat, heating = heating)
    end

    Kₛ = zeros(Float64, n_layers)
    for i in 1:n_layers
        Kₛ[i] = _snow_thermal_conductivity(ρ[i], Kᵢ, diff_model)
    end

    @inbounds begin
        # i = 1 boundary
        k₁₂ = interface_conductance(Kₛ[1], Δz[1], Kₛ[2], Δz[2])
        a_upper[1] = -2.0 * dt_sec / (_safe_positive(ρ[1]) * cᵢ * _safe_positive(Δz[1])) * k₁₂
        a_diag[1] = 1.0 - a_upper[1]

        # i = n_layers boundary
        kₙ = interface_conductance(Kₛ[n_layers], Δz[n_layers], Kₛ[n_layers - 1], Δz[n_layers - 1])
        a_lower[n_layers] = -2.0 * dt_sec / (_safe_positive(ρ[n_layers]) * cᵢ * _safe_positive(Δz[n_layers])) * kₙ
        a_diag[n_layers] = 1.0 - a_lower[n_layers]

        # interior
        for i in 2:n_layers-1
            k_lower = interface_conductance(Kₛ[i], Δz[i], Kₛ[i - 1], Δz[i - 1])
            k_upper = interface_conductance(Kₛ[i], Δz[i], Kₛ[i + 1], Δz[i + 1])
            a_lower[i] = -2.0 * dt_sec / (_safe_positive(ρ[i]) * cᵢ * _safe_positive(Δz[i])) * k_lower
            a_upper[i] = -2.0 * dt_sec / (_safe_positive(ρ[i]) * cᵢ * _safe_positive(Δz[i])) * k_upper
            a_diag[i] = 1.0 - a_lower[i] - a_upper[i]
        end
    end

    a_diag₁_backup = a_diag[1]
    a_diag[1] += surface_diag_term

    rhs = copy(T)
    rhs[1] += surface_rhs_term
    A = Tridiagonal(a_lower[2:end], a_diag, a_upper[1:(end - 1)])
    T_new = A \ rhs

    if T_new[1] > Tₘ
        china_syndrome = true
        Q_heat = max((Tₘ - T_prev[1]) * cᵢ * mₛ, 0.0)

        rhs_melt = copy(T_prev)
        rhs_melt[1] = Tₘ
        a_diag[1] = a_diag₁_backup
        A = Tridiagonal(a_lower[2:end], a_diag, a_upper[1:(end - 1)])
        T_new = A \ rhs_melt

        Q_heat += max((Tₘ - T_new[1]) * cᵢ * mₛ, 0.0)
        T_new[1] = Tₘ
        _clamp_to_melt!(T_new, Tₘ)
        heating = Q_heat
    else
        _clamp_to_melt!(T_new, Tₘ)
        heating = dt_sec * (F_const - F_lin * T_new[1])
    end

    column.temperature[1:n_layers] .= T_new
    column.Tsrf = column.temperature[1]

    return (china_syndrome = china_syndrome, Q_heat = Q_heat, heating = heating)
end
