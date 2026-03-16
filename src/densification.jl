"""
Firn densification translated from BESSI `go_densification`.
"""

@inline function _overburden_pressure(snowman, mm::Int)
    columnsnow = mm > 1 ? sum(@view snowman[1:mm-1]) : 0.0
    return 9.81 * (columnsnow + snowman[mm] / 2.0)
end

@inline _bessi_low_density_rate(ρ::Float64, T::Float64, ρᵢ::Float64, At::Float64) =
    0.011 * exp(-10160.0 / 8.13 / T) * (ρᵢ - ρ) * max(At, 0.0)

@inline _htessel_snow_viscosity(T0::Float64, T::Float64, ρ::Float64) =
    3.7e7 * exp(8.1e-2 * (T0 - T) + 1.8e-2 * ρ)

@inline _htessel_thermal_metamorphism(T0::Float64, T::Float64, ρ::Float64) =
    2.8e-6 * exp(-4.2e-2 * (T0 - T) - 460.0 * max(0.0, ρ - 150.0))

@inline function _htessel_low_density_rate(
    column::SnowpackColumn,
    snowman,
    mm::Int,
    T::Float64,
    ρ::Float64,
)
    σ = _overburden_pressure(snowman, mm)
    η = _htessel_snow_viscosity(column.c.T0, T, ρ)
    ξ = _htessel_thermal_metamorphism(column.c.T0, T, ρ)
    return ρ * (σ / η + ξ)
end

function _apply_htessel_liquid_water_compaction!(
    column::SnowpackColumn,
    liquid_water_before_energy::AbstractVector{Float64},
    dt_sec::Float64,
)
    n = column.N
    if n <= 0 || dt_sec <= 0.0
        return
    end

    @views snowman = column.mass[1:n]
    @views rho = column.density[1:n]
    @views mass_w = column.mass_w[1:n]
    ρᵢ = column.c.rho_i

    for mm in 1:n
        ρ = rho[mm]
        m = snowman[mm]
        if ρ < 550.0 && m > EPS_TINY
            prev = mm <= length(liquid_water_before_energy) ? liquid_water_before_energy[mm] : 0.0
            gain = max(mass_w[mm] - prev, 0.0)
            if gain > 0.0
                # HTESSEL eq. (8): retained meltwater increases density at fixed solid mass.
                rho[mm] = min(max(ρ, ρ + ρ * gain / m), ρᵢ)
            end
        end
    end
    return
end

"""
    go_densification!(column::SnowpackColumn, At::Float64, dt_sec::Float64;
                      hl::Bool=false, rho_e::Float64=815.0, P_atm::Float64=101325.0)

Update layer densities in-place using the multi-regime densification scheme
from the original Fortran implementation.

- `At` is the accumulation proxy used by the HL branch.
- `dt_sec` is timestep in seconds.
- `hl=true` forces HL for densities above 550 kg m^-3.
"""
function go_densification!(
    column::SnowpackColumn,
    At::Float64,
    dt_sec::Float64;
    hl::Bool=false,
    rho_e::Float64=815.0,
    P_atm::Float64=101325.0,
)
    n = column.N
    if n <= 0
        return
    end

    ρᵢ = column.c.rho_i
    @views snowman = column.mass[1:n]
    @views temp = column.temperature[1:n]
    @views rho = column.density[1:n]
    low_density_scheme = column.c.low_density_densification

    for mm in 1:n
        ρ = rho[mm]
        T = temp[mm]
        m = snowman[mm]
        ddens = 0.0

        if ρ < 550.0
            if m > 0.0
                if low_density_scheme == :htessel
                    ddens = _htessel_low_density_rate(column, snowman, mm, T, ρ)
                else
                    ddens = _bessi_low_density_rate(ρ, T, ρᵢ, At)
                end
            end
        elseif hl
            if mm > 1 && m > 0.0
                ddens = 0.575 * exp(-21400.0 / 8.13 / T) * (ρᵢ - ρ) *
                        (1000.0 / 3600.0 / 24.0 / 365.0)^0.5 * max(At, 0.0)^0.5
            end
        elseif ρ < 800.0
            f = 10.0^(-29.166 * (ρ / ρᵢ)^3 + 84.422 * (ρ / ρᵢ)^2 - 87.425 * (ρ / ρᵢ) + 30.673)
            columnsnow = mm > 1 ? sum(snowman[1:mm-1]) : 0.0
            P_ice = 9.81 * (columnsnow + m / 2.0) / 1.0e6
            P_bubble = 0.0
            if ρ > rho_e
                P_bubble = P_atm * ((1.0 / rho_e - 1.0 / ρᵢ) / (1.0 / ρ - 1.0 / ρᵢ) - 1.0) / 1.0e6
            end
            dp = P_ice - P_bubble
            ddens = 25400.0 * exp(-60000.0 / 8.13 / T) * ρ * f * dp^3
        else
            denom = 1.0 - (1.0 - ρ / ρᵢ)^(1.0 / 3.0)
            if abs(denom) > EPS_TINY
                f = 3.0 / 16.0 * (1.0 - ρ / ρᵢ) / denom^3
                columnsnow = mm > 1 ? sum(snowman[1:mm-1]) : 0.0
                P_ice = 9.81 * (columnsnow + m / 2.0) / 1.0e6
                P_bubble = 0.0
                if ρ > rho_e
                    P_bubble = P_atm * ((1.0 / rho_e - 1.0 / ρᵢ) / (1.0 / ρ - 1.0 / ρᵢ) - 1.0) / 1.0e6
                end
                dp = P_ice - P_bubble
                ddens = 25400.0 * exp(-60000.0 / 8.13 / T) * ρ * f * dp^3
            end
        end

        rho_new = max(ρ, ρ + ddens * dt_sec)
        rho[mm] = min(rho_new, ρᵢ)
    end

    return
end
