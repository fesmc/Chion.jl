"""
Liquid-water refreezing for array-backed snowpack states.
"""

"""
    _go_refreezing!(N_storage, mass_w, mass, density, temperature, idx, melting_temperature, ice_heat_capacity, latent_heat_of_melting, ice_density)

Refreeze liquid water in subfreezing layers of column `idx` until either cold
content or liquid water is exhausted. Mutates liquid water, solid mass,
density, and temperature in-place and returns the refrozen mass.
"""
function _go_refreezing!(
    N_storage,
    mass_w,
    mass,
    density,
    temperature,
    idx::Int,
    melting_temperature,
    ice_heat_capacity,
    latent_heat_of_melting,
    ice_density,
)
    refrozen_mass = zero(latent_heat_of_melting)

    @inbounds for layer_index in 1:_n_active(N_storage, idx)
        solid_mass = _get_layer(mass, layer_index, idx)
        liquid_water_mass = _get_layer(mass_w, layer_index, idx)
        layer_temperature = _get_layer(temperature, layer_index, idx)
        if solid_mass > zero(solid_mass) &&
           liquid_water_mass > zero(liquid_water_mass) &&
           layer_temperature < melting_temperature
            cold_content = (melting_temperature - layer_temperature) * ice_heat_capacity * solid_mass
            available_latent_heat = liquid_water_mass * latent_heat_of_melting

            if cold_content < available_latent_heat
                newly_refrozen_mass = cold_content / latent_heat_of_melting
                _set_layer!(temperature, layer_index, idx, melting_temperature)
                _set_layer!(
                    density,
                    layer_index,
                    idx,
                    min(_get_layer(density, layer_index, idx) * (newly_refrozen_mass + solid_mass) / solid_mass, ice_density),
                )
                _set_layer!(mass, layer_index, idx, solid_mass + newly_refrozen_mass)
                _set_layer!(mass_w, layer_index, idx, liquid_water_mass - newly_refrozen_mass)
                refrozen_mass += newly_refrozen_mass
            else
                updated_temperature = (
                    liquid_water_mass * latent_heat_of_melting / ice_heat_capacity +
                    liquid_water_mass * melting_temperature +
                    layer_temperature * solid_mass
                ) / (liquid_water_mass + solid_mass)
                _set_layer!(temperature, layer_index, idx, updated_temperature)
                _set_layer!(
                    density,
                    layer_index,
                    idx,
                    min(_get_layer(density, layer_index, idx) * (liquid_water_mass + solid_mass) / solid_mass, ice_density),
                )
                _set_layer!(mass, layer_index, idx, solid_mass + liquid_water_mass)
                refrozen_mass += liquid_water_mass
                _set_layer!(mass_w, layer_index, idx, zero(liquid_water_mass))
            end
        end
    end

    return refrozen_mass
end

"""
    go_refreezing!(liquid_water_mass, solid_mass, snow_density, layer_temperature, melting_temperature, ice_heat_capacity, latent_heat_of_melting, ice_density)

Run the refreezing scheme on vector-backed column data. Mutates the supplied
arrays in-place and returns the refrozen mass together with released latent
heat.
"""
function go_refreezing!(
    liquid_water_mass::AbstractVector,
    solid_mass::AbstractVector,
    snow_density::AbstractVector,
    layer_temperature::AbstractVector,
    melting_temperature,
    ice_heat_capacity,
    latent_heat_of_melting,
    ice_density,
)
    N_ref = Ref(length(solid_mass))
    refrozen_mass = _go_refreezing!(
        N_ref,
        liquid_water_mass,
        solid_mass,
        snow_density,
        layer_temperature,
        1,
        melting_temperature,
        ice_heat_capacity,
        latent_heat_of_melting,
        ice_density,
    )
    return (
        refrozen_mass=refrozen_mass,
        released_latent_heat=refrozen_mass * latent_heat_of_melting,
    )
end

"""
    go_refreezing!(domain, idx)

Run the refreezing scheme for column `idx` of `domain`. Mutates the domain
state in-place and returns the refrozen mass and released latent heat.
"""
function go_refreezing!(domain, idx::Int)
    if _n_active(domain.N, idx) <= 0
        return (
            refrozen_mass=zero(eltype(domain.mass)),
            released_latent_heat=zero(eltype(domain.mass)),
        )
    end

    refrozen_mass = _go_refreezing!(
        domain.N,
        domain.mass_w,
        domain.mass,
        domain.density,
        domain.temperature,
        idx,
        domain.c.T0,
        domain.c.ci,
        domain.c.Lm,
        domain.c.rho_i,
    )
    return (
        refrozen_mass=refrozen_mass,
        released_latent_heat=refrozen_mass * domain.c.Lm,
    )
end
