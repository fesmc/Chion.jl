"""
Surface albedo state and update rules for array-backed snow states.
"""

"""
    _constant_surface_albedo(N_storage, mass, temperature, idx, c)

Return the constant-scheme surface albedo for column `idx`, switching between
ice, dry snow, and wet snow according to surface state.
"""
@inline function _constant_surface_albedo(
    N_storage,
    mass,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
        return c.alpha_ice
    end
    return _get_layer(temperature, 1, idx) >= c.T0 ? c.alpha_wet : c.alpha_dry
end

"""
    _surface_liquid_water_content(N_storage, mass, mass_w, density, idx, c)

Estimate volumetric liquid-water content in the surface layer of column `idx`.
Returns zero when the layer is empty, ice-dense, or has no pore volume.
"""
@inline function _surface_liquid_water_content(
    N_storage,
    mass,
    mass_w,
    density,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_TINY
        return zero(eltype(mass))
    end

    surface_density = _get_layer(density, 1, idx)
    if surface_density <= EPS_TINY || surface_density >= c.rho_i - EPS_TINY
        return zero(surface_density)
    end

    pore_volume = _get_layer(mass, 1, idx) / surface_density - _get_layer(mass, 1, idx) / c.rho_i
    if pore_volume <= EPS_TINY
        return zero(pore_volume)
    end

    return max(_get_layer(mass_w, 1, idx), zero(eltype(mass))) / c.rho_w / pore_volume
end

"""
    _refresh_dynamic_albedo_from_snowfall!(albedo_dynamic, idx, c, snowfall_mass)

Refresh the dynamic surface albedo after snowfall on column `idx`. Mutates
`albedo_dynamic[idx]` and returns the updated albedo.
"""
function _refresh_dynamic_albedo_from_snowfall!(
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    snowfall_mass,
)
    snowfall_mass <= EPS_TINY && return _get_scalar(albedo_dynamic, idx)

    if _uses_constant_albedo(c)
        _set_scalar!(albedo_dynamic, idx, c.alpha_dry)
        return _get_scalar(albedo_dynamic, idx)
    end

    updated = min(
        c.alpha_dry,
        _get_scalar(albedo_dynamic, idx) +
        (c.alpha_dry - c.alpha_wet) * (one(snowfall_mass) - exp(-snowfall_mass / oftype(snowfall_mass, 3))),
    )
    _set_scalar!(albedo_dynamic, idx, updated)
    return updated
end

"""
    _update_surface_albedo_arrays!(N_storage, mass, mass_w, density, temperature, albedo_dynamic, idx, c)

Update the diagnosed surface albedo for column `idx` in-place using either the
constant or dynamic albedo scheme.
"""
function _update_surface_albedo_arrays!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        return c.alpha_ice
    end

    if _uses_constant_albedo(c)
        updated = _constant_surface_albedo(N_storage, mass, temperature, idx, c)
        _set_scalar!(albedo_dynamic, idx, updated)
        return updated
    end

    previous_albedo = clamp(_get_scalar(albedo_dynamic, idx), c.alpha_wet, c.alpha_dry)
    surface_temperature = _get_layer(temperature, 1, idx)
    updated_albedo = min(
        previous_albedo,
        previous_albedo - ((surface_temperature - c.T0) * oftype(surface_temperature, 1.35e-3) + oftype(surface_temperature, 0.0278)),
    )
    updated_albedo = max(updated_albedo, c.alpha_wet)

    liquid_water_content = _surface_liquid_water_content(N_storage, mass, mass_w, density, idx, c)
    if liquid_water_content > zero(liquid_water_content) && c.max_lwc_albedo > EPS_TINY
        wet_adjusted_albedo = updated_albedo - (
            updated_albedo - c.alpha_wet
        ) * (liquid_water_content / c.max_lwc_albedo)
        updated_albedo = max(c.alpha_wet, min(updated_albedo, wet_adjusted_albedo))
    end

    updated_albedo = clamp(updated_albedo, c.alpha_wet, c.alpha_dry)
    _set_scalar!(albedo_dynamic, idx, updated_albedo)
    return updated_albedo
end

"""
    update_surface_albedo!(state, idx)

Update the surface albedo of column `idx` in `state` and return the new
albedo.
"""
function update_surface_albedo!(state, idx::Int)
    return _update_surface_albedo_arrays!(
        state.N,
        state.mass,
        state.mass_w,
        state.density,
        state.temperature,
        state.albedo,
        idx,
        state.c,
    )
end
