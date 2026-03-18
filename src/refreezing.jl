"""
Liquid-water refreezing following the original BESSI-style routine.
"""

"""
    go_refreezing!(
        lwmass::AbstractVector{Float64},
        snowman::AbstractVector{Float64},
        rho_snow::AbstractVector{Float64},
        snow_temp::AbstractVector{Float64},
        kelvin::Float64,
        c_i::Float64,
        L_lh::Float64,
        rho_i::Float64,
    ) -> NamedTuple

Translate the Fortran `go_refreezing` logic to Julia.

Returns:
- `refreeze`: refrozen water mass [kg m^-2]
- `heat_fusion`: latent heat released by refreezing [J m^-2]
"""
function go_refreezing!(
    lwmass::AbstractVector{Float64},
    snowman::AbstractVector{Float64},
    rho_snow::AbstractVector{Float64},
    snow_temp::AbstractVector{Float64},
    kelvin::Float64,
    c_i::Float64,
    L_lh::Float64,
    rho_i::Float64,
)
    n_snowlayer = length(snowman)
    @assert length(lwmass) == n_snowlayer
    @assert length(rho_snow) == n_snowlayer
    @assert length(snow_temp) == n_snowlayer

    refreeze = 0.0
    heat_fusion = 0.0

    for ii in 1:n_snowlayer
        if snowman[ii] > 0.0 && lwmass[ii] > 0.0
            cold_content = (kelvin - snow_temp[ii]) * c_i * snowman[ii]
            latent_available = lwmass[ii] * L_lh

            if cold_content < latent_available
                # CASE 1: water freezes partly
                icecube = cold_content / L_lh

                heat_fusion += icecube * L_lh
                snow_temp[ii] = kelvin
                rho_snow[ii] = min(rho_snow[ii] * (icecube + snowman[ii]) / snowman[ii], rho_i)
                snowman[ii] += icecube
                lwmass[ii] -= icecube
                refreeze += icecube
            else
                # CASE 2: all water freezes
                snow_temp[ii] = (
                    lwmass[ii] * L_lh / c_i +
                    lwmass[ii] * kelvin +
                    snow_temp[ii] * snowman[ii]
                ) / (lwmass[ii] + snowman[ii])

                rho_snow[ii] = min(rho_snow[ii] * (lwmass[ii] + snowman[ii]) / snowman[ii], rho_i)
                snowman[ii] += lwmass[ii]
                refreeze += lwmass[ii]
                heat_fusion += lwmass[ii] * L_lh
                lwmass[ii] = 0.0
            end
        end
    end

    return (refreeze = refreeze, heat_fusion = heat_fusion)
end

"""
    go_refreezing!(column::SnowpackColumn) -> NamedTuple

Apply refreezing to active `SnowpackColumn` layers.
"""
function go_refreezing!(column::SnowpackColumn)
    n = column.N
    if n <= 0
        return (refreeze = 0.0, heat_fusion = 0.0)
    end

    @views return go_refreezing!(
        column.mass_w[1:n],
        column.mass[1:n],
        column.density[1:n],
        column.temperature[1:n],
        column.c.T0,
        column.c.ci,
        column.c.Lm,
        column.c.rho_i,
    )
end
