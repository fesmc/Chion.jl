"""
Firn densification translated from BESSI `go_densification`.
"""

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

    for mm in 1:n
        ρ = rho[mm]
        T = temp[mm]
        m = snowman[mm]
        ddens = 0.0

        if ρ < 550.0
            if m > 0.0
                ddens = 0.011 * exp(-10160.0 / 8.13 / T) * (ρᵢ - ρ) * max(At, 0.0)
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
