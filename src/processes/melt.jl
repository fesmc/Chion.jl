"""
Surface melt for array-backed snowpack states.
"""

"""
    _melt_surface_layers!(..., idx, mass_split, mass_min, melt_mass, c, charge_sensible_heat)

Melt surface snow of column `idx` and convert it to liquid water, removing or
merging depleted layers. `melt_mass` is the melt potential in kg m⁻² of ice at
the melting point (the melt energy divided by `Lm`).

With `charge_sensible_heat = true` (melt computed from the surface energy
balance), ice taken from a layer below T₀ first has to be warmed to T₀, so each
kilogram costs `Lm + ci (T₀ − T)` and the melted mass is smaller than
`melt_mass`. This keeps the column energy budget closed. With
`charge_sensible_heat = false` (prescribed melt), exactly `melt_mass` is melted.

Returns `(melted, remaining)`: the melted snow mass and the melt potential left
when the snow is exhausted (in kg m⁻² of ice at T₀).
"""
function _melt_surface_layers!(
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
    charge_sensible_heat::Bool,
)
    remaining_melt = _safe_nonnegative(melt_mass)
    melted_total = zero(remaining_melt)
    if remaining_melt <= zero(remaining_melt) || _n_active(N_storage, idx) <= 0
        return melted_total, remaining_melt
    end

    while remaining_melt > zero(remaining_melt) && _n_active(N_storage, idx) > 0
        layer_mass = _get_layer(mass, 1, idx)
        if layer_mass <= EPS_TINY
            _remove_depleted_surface_and_route_water!(N_storage, mass, mass_w, density, temperature, runoff, idx, c)
            continue
        end

        # Melt potential consumed per kilogram of ice from this layer.
        cost = ifelse(
            charge_sensible_heat,
            one(remaining_melt) + c.ci * max(c.T0 - _get_layer(temperature, 1, idx), zero(remaining_melt)) / c.Lm,
            one(remaining_melt),
        )
        affordable = remaining_melt / cost
        if affordable < layer_mass - EPS_EMPTY_LAYER
            dm = affordable
            remaining_melt = zero(remaining_melt)
        else
            dm = layer_mass
            remaining_melt = max(remaining_melt - dm * cost, zero(remaining_melt))
        end
        _set_layer!(mass, 1, idx, layer_mass - dm)
        _set_layer!(mass_w, 1, idx, _get_layer(mass_w, 1, idx) + dm)
        melted_total += dm

        if _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
            _remove_depleted_surface_and_route_water!(N_storage, mass, mass_w, density, temperature, runoff, idx, c)
        elseif _n_active(N_storage, idx) > 1 && _get_layer(mass, 1, idx) < mass_min
            _merge_surface_layer!(N_storage, mass, mass_w, density, temperature, idx, typemax(Int), mass_split, mass_min, c)
        end
    end

    _set_scalar!(Tsrf, idx, c.T0)
    if _n_active(N_storage, idx) == 0
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end
    return melted_total, remaining_melt
end

"""
    _apply_melt!(..., idx, mass_split, mass_min, melt_mass, c)

Convert up to `melt_mass` of surface snow into liquid water in column `idx`
(prescribed-melt semantics: no sensible-heat charge) and return the melted
mass. See `_melt_surface_layers!`.
"""
_apply_melt!(N_storage, mass, mass_w, density, temperature, runoff, Tsrf, albedo_dynamic, idx::Int,
    mass_split, mass_min, melt_mass, c::SnowpackPhysicalConstants) =
    first(_melt_surface_layers!(N_storage, mass, mass_w, density, temperature, runoff, Tsrf, albedo_dynamic,
        idx, mass_split, mass_min, melt_mass, c, false))
