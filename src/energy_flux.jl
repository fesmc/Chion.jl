#=
Energy-flux temperature solver for `SnowpackColumn`.

Port of BESSI `go_energy_flux_new` to the Chion column nomenclature.
=#
using LinearAlgebra
@inline function _safe_positive(x::Float64)
    return x > 1.0e-12 ? x : 1.0e-12
end

function _snow_thermal_conductivity(rho::Float64, K_ice::Float64, diff_model::Int)
    if diff_model == 1
        # Yen (1981)
        return K_ice * (rho / 1000.0)^1.88
    elseif diff_model == 2
        # Sturm (1997), piecewise form
        if rho > 156.0
            return 0.138 - 1.01e-3 * rho + 3.233e-6 * rho^2
        else
            return 0.023 + 0.234e-3 * rho
        end
    else
        # Van Dusen (1929)
        return 2.1e-2 + 4.2e-4 * rho + 2.2e-9 * rho^3
    end
end

@inline function shortwave_absorbed(S_boa::Float64, Ts::Float64, Tmelt::Float64;
                                   alpha_dry::Float64=0.8,
                                   alpha_wet::Float64=0.6)
    # wet when at/near melting
    α = (Ts >= Tmelt) ? alpha_wet : alpha_dry
    return (1.0 - α) * S_boa
end
@inline interfaceK(Ki, dzi, Kj, dzj) = (Ki*dzi + Kj*dzj) / _safe_positive((dzi + dzj)^2)
"""
    go_energy_flux!(
        column::SnowpackColumn,
        T2m::Float64,
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
    T2m::Float64,
    S_boa::Float64,
    H_lh::Float64,
    K_lh::Float64,
    dt_sec::Float64;
    diff_model::Int = 2,
)
    nn = column.N
    if nn <= 0 || column.mass[1] <= 0.0
        return (china_syndrome = false, Q_heat = 0.0, heating = 0.0)
    end

    kelvin = column.c.T0
    ci = column.c.ci
    D_sf = column.c.D_sh
    sigma = column.c.σ
    eps_air = column.c.ϵ_air
    eps_snow = column.c.ϵ_snow
    K_ice = column.c.Ki
   
    BB_up = zeros(Float64, nn)
    BB_down = zeros(Float64, nn)
    BB_mid = zeros(Float64, nn)

    temps = copy(column.temperature[1:nn])
    backup = copy(temps)
    rho = copy(column.density[1:nn])
    dz = zeros(Float64, nn)
    for i in 1:nn
        dz[i] = column.mass[i] / _safe_positive(rho[i])
    end

    mbox = _safe_positive(column.mass[1])
    inv = dt_sec / ci / mbox
    Qsw = shortwave_absorbed(S_boa, temps[1], kelvin; alpha_dry=column.c.alpha_dry, alpha_wet=column.c.alpha_wet)
    K1 = inv * (
        T2m * D_sf +
        sigma * (eps_air * T2m^4 + eps_snow * 3.0 * temps[1]^4) +
        Qsw +
        K_lh
    )
    H = inv * (D_sf + sigma * eps_snow * 4.0 * temps[1]^3 + H_lh)

    china_syndrome = false
    Q_heat = 0.0
    heating = 0.0

    if nn == 1
        new_temp = (temps[1] + K1) / _safe_positive(1.0 + H)
        if new_temp > kelvin
            china_syndrome = true
            Q_heat = (kelvin - backup[1]) * ci * mbox
            new_temp = kelvin
            heating = Q_heat
        else
            heating = dt_sec * (
                new_temp * (-D_sf - sigma * eps_snow * 4.0 * temps[1]^3 - H_lh) +
                (T2m * D_sf + sigma * (eps_air * T2m^4 + eps_snow * 3.0 * temps[1]^4) + Qsw + K_lh)
            )
        end

        column.temperature[1] = min(new_temp, kelvin)
        column.Tsrf = column.temperature[1]
        return (china_syndrome = china_syndrome, Q_heat = Q_heat, heating = heating)
    end

    K_snow = zeros(Float64, nn)
    for i in 1:nn
        K_snow[i] = _snow_thermal_conductivity(rho[i], K_ice, diff_model)
    end

    @inbounds begin
            # i = 1 boundary
            k12 = interfaceK(K_snow[1], dz[1], K_snow[2], dz[2])
            BB_up[1] = -2.0 * dt_sec / (_safe_positive(rho[1]) * ci * _safe_positive(dz[1])) * k12
            BB_mid[1] = 1.0 - BB_up[1]

            # i = nn boundary
            kn = interfaceK(K_snow[nn], dz[nn], K_snow[nn-1], dz[nn-1])
            BB_down[nn] = -2.0 * dt_sec / (_safe_positive(rho[nn]) * ci * _safe_positive(dz[nn])) * kn
            BB_mid[nn] = 1.0 - BB_down[nn]

            # interior
            for i in 2:nn-1
                kdn = interfaceK(K_snow[i], dz[i], K_snow[i-1], dz[i-1])
                kup = interfaceK(K_snow[i], dz[i], K_snow[i+1], dz[i+1])
                BB_down[i] = -2.0 * dt_sec / (_safe_positive(rho[i]) * ci * _safe_positive(dz[i])) * kdn
                BB_up[i]   = -2.0 * dt_sec / (_safe_positive(rho[i]) * ci * _safe_positive(dz[i])) * kup
                BB_mid[i]  = 1.0 - BB_down[i] - BB_up[i]
            end
        end

    BB_mid1_backup = BB_mid[1]
    BB_mid[1] += H

    rhs = copy(temps)
    rhs[1] += K1
    A = Tridiagonal(BB_down[2:end], BB_mid, BB_up[1:(end - 1)])
    new_temp = A \ rhs
    #new_temp = _solve_tridiagonal(BB_down, BB_mid, BB_up, rhs)

    if new_temp[1] > kelvin
        china_syndrome = true
        Q_heat = (kelvin - backup[1]) * ci * mbox

        rhs2 = copy(backup)
        rhs2[1] = kelvin
        BB_mid[1] = BB_mid1_backup
        A = Tridiagonal(BB_down[2:end], BB_mid, BB_up[1:(end - 1)])
        new_temp = A \ rhs2
        

        Q_heat += (kelvin - new_temp[1]) * ci * mbox
        new_temp[1] = kelvin
        for i in 1:nn
            if new_temp[i] > kelvin
                new_temp[i] = kelvin
            end
        end
        heating = Q_heat
    else
        for i in 1:nn
            if new_temp[i] > kelvin
                new_temp[i] = kelvin
            end
        end
        heating = dt_sec * (
            new_temp[1] * (-D_sf - sigma * eps_snow * 4.0 * temps[1]^3 - H_lh) +
            (T2m * D_sf + sigma * (eps_air * T2m^4 + eps_snow * 3.0 * temps[1]^4) + Qsw + K_lh)
        )
    end

    column.temperature[1:nn] .= new_temp
    column.Tsrf = column.temperature[1]

    return (china_syndrome = china_syndrome, Q_heat = Q_heat, heating = heating)
end
