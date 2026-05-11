"""
Surface melt for array-backed snowpack states.
"""

"""
    _apply_melt!(..., idx, mass_split, mass_min, melt_mass, c)

Convert up to `melt_mass` of surface snow into liquid water in column `idx`.
This mutates snow mass, liquid water, runoff routing, surface temperature, and
surface albedo as depleted layers are removed or merged.
"""
function _apply_melt!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    runoff,
    Tsrf,
    albedo_dynamic,
    idx::Int,
    mass_split,
    mass_min,
    melt_mass,
    c::SnowpackPhysicalConstants,
)
    remaining_melt = _safe_nonnegative(melt_mass)
    if remaining_melt <= zero(remaining_melt) || _n_active(N_storage, idx) <= 0
        return zero(remaining_melt)
    end

    surface_mass = _get_layer(mass, 1, idx)
    if remaining_melt < surface_mass - EPS_EMPTY_LAYER
        _set_layer!(mass, 1, idx, surface_mass - remaining_melt)
        _set_layer!(mass_w, 1, idx, _get_layer(mass_w, 1, idx) + remaining_melt)
        return remaining_melt
    end

    melted_total = zero(remaining_melt)
    while remaining_melt > zero(remaining_melt) && _n_active(N_storage, idx) > 0
        layer_mass = _get_layer(mass, 1, idx)
        if layer_mass <= EPS_TINY
            _remove_depleted_surface_and_route_water!(N_storage, mass, mass_w, density, temperature, runoff, idx, c)
            continue
        end

        dm = min(layer_mass, remaining_melt)
        _set_layer!(mass, 1, idx, layer_mass - dm)
        _set_layer!(mass_w, 1, idx, _get_layer(mass_w, 1, idx) + dm)
        remaining_melt -= dm
        melted_total += dm

        if _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
            _remove_depleted_surface_and_route_water!(N_storage, mass, mass_w, density, temperature, runoff, idx, c)
        elseif _n_active(N_storage, idx) > 1 && _get_layer(mass, 1, idx) < mass_min
            _merge_surface_layer!(N_storage, mass, mass_w, density, temperature, idx, typemax(Int), mass_split, mass_min, c)
        end
    end

    _set_scalar!(Tsrf, idx, _n_active(N_storage, idx) > 0 ? _get_layer(temperature, 1, idx) : c.T0)
    if _n_active(N_storage, idx) == 0
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end
    return melted_total
end
