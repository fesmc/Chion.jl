"""
Liquid-water refreezing following the original BESSI-style routine.
"""

"""
    go_refreezing!(
        liquid_water_mass::AbstractVector{Float64},
        solid_mass::AbstractVector{Float64},
        snow_density::AbstractVector{Float64},
        layer_temperature::AbstractVector{Float64},
        melting_temperature::Float64,
        ice_heat_capacity::Float64,
        latent_heat_of_melting::Float64,
        ice_density::Float64,
    ) -> NamedTuple

Translate the Fortran `go_refreezing` logic to Julia.

Returns:
- `refreeze` / `refrozen_mass`: refrozen water mass [kg m^-2]
- `heat_fusion` / `released_latent_heat`: latent heat released by refreezing [J m^-2]
"""
function go_refreezing!(
    liquid_water_mass::AbstractVector{Float64},
    solid_mass::AbstractVector{Float64},
    snow_density::AbstractVector{Float64},
    layer_temperature::AbstractVector{Float64},
    melting_temperature::Float64,
    ice_heat_capacity::Float64,
    latent_heat_of_melting::Float64,
    ice_density::Float64,
)
    n_layers = length(solid_mass)
    @assert length(liquid_water_mass) == n_layers
    @assert length(snow_density) == n_layers
    @assert length(layer_temperature) == n_layers

    refrozen_mass = 0.0
    released_latent_heat = 0.0

    for layer_index in 1:n_layers
        if solid_mass[layer_index] > 0.0 && liquid_water_mass[layer_index] > 0.0
            cold_content = (
                melting_temperature - layer_temperature[layer_index]
            ) * ice_heat_capacity * solid_mass[layer_index]
            available_latent_heat = liquid_water_mass[layer_index] * latent_heat_of_melting

            if cold_content < available_latent_heat
                # CASE 1: water freezes partly
                newly_refrozen_mass = cold_content / latent_heat_of_melting

                released_latent_heat += newly_refrozen_mass * latent_heat_of_melting
                layer_temperature[layer_index] = melting_temperature
                snow_density[layer_index] = min(
                    snow_density[layer_index] *
                    (newly_refrozen_mass + solid_mass[layer_index]) /
                    solid_mass[layer_index],
                    ice_density,
                )
                solid_mass[layer_index] += newly_refrozen_mass
                liquid_water_mass[layer_index] -= newly_refrozen_mass
                refrozen_mass += newly_refrozen_mass
            else
                # CASE 2: all water freezes
                layer_temperature[layer_index] = (
                    liquid_water_mass[layer_index] * latent_heat_of_melting / ice_heat_capacity +
                    liquid_water_mass[layer_index] * melting_temperature +
                    layer_temperature[layer_index] * solid_mass[layer_index]
                ) / (liquid_water_mass[layer_index] + solid_mass[layer_index])

                snow_density[layer_index] = min(
                    snow_density[layer_index] *
                    (liquid_water_mass[layer_index] + solid_mass[layer_index]) /
                    solid_mass[layer_index],
                    ice_density,
                )
                solid_mass[layer_index] += liquid_water_mass[layer_index]
                refrozen_mass += liquid_water_mass[layer_index]
                released_latent_heat += liquid_water_mass[layer_index] * latent_heat_of_melting
                liquid_water_mass[layer_index] = 0.0
            end
        end
    end

    return (
        refreeze = refrozen_mass,
        refrozen_mass = refrozen_mass,
        heat_fusion = released_latent_heat,
        released_latent_heat = released_latent_heat,
    )
end

"""
    go_refreezing!(column::SnowpackColumn) -> NamedTuple

Apply refreezing to active `SnowpackColumn` layers.
"""
function go_refreezing!(column::SnowpackColumn)
    n = column.N
    if n <= 0
        return (
            refreeze = 0.0,
            refrozen_mass = 0.0,
            heat_fusion = 0.0,
            released_latent_heat = 0.0,
        )
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
