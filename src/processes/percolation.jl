"""
Liquid-water percolation for column and domain states.
"""

"""
    _is_lowest_active_snow_layer(N_storage, mass, idx, layer_index)

Return `true` when `layer_index` is the last active snow layer in column
`idx`.
"""
@inline function _is_lowest_active_snow_layer(
    N_storage,
    mass,
    idx::Int,
    layer_index::Int,
)
    return layer_index == _n_active(N_storage, idx) || _get_layer(mass, layer_index + 1, idx) <= zero(eltype(mass))
end

"""
    _go_percolation!(N_storage, mass, mass_w, density, idx, ice_density, water_density; max_lwc=0.1, rho_i_tol=10.0)

Route excess liquid water downward through column `idx` until each layer is at
or below the liquid-water-content threshold. Mutates `mass_w` in-place and
returns runoff leaving the bottom of the column.
"""
function _go_percolation!(
    N_storage,
    mass,
    mass_w,
    density,
    idx::Int,
    ice_density,
    water_density;
    max_lwc=oftype(water_density, 0.1),
    rho_i_tol=oftype(water_density, 10.0),
)
    runoff = zero(water_density)
    n_layers = _n_active(N_storage, idx)

    @inbounds for layer_index in 1:n_layers
        solid_mass = _get_layer(mass, layer_index, idx)
        liquid_water_mass = _get_layer(mass_w, layer_index, idx)

        if solid_mass <= zero(solid_mass)
            runoff += liquid_water_mass
            _set_layer!(mass_w, layer_index, idx, zero(liquid_water_mass))
            continue
        end

        pore_volume = solid_mass / _get_layer(density, layer_index, idx) - solid_mass / ice_density
        if pore_volume <= EPS_TINY
            excess_water = liquid_water_mass
            _set_layer!(mass_w, layer_index, idx, zero(liquid_water_mass))
            if _is_lowest_active_snow_layer(N_storage, mass, idx, layer_index)
                runoff += excess_water
            else
                _set_layer!(mass_w, layer_index + 1, idx, _get_layer(mass_w, layer_index + 1, idx) + excess_water)
            end
            continue
        end

        liquid_water_content = liquid_water_mass / water_density / pore_volume
        if liquid_water_content > max_lwc
            excess_water = (liquid_water_content - max_lwc) * pore_volume * water_density
            _set_layer!(mass_w, layer_index, idx, liquid_water_mass - excess_water)
            if _is_lowest_active_snow_layer(N_storage, mass, idx, layer_index)
                runoff += excess_water
            else
                _set_layer!(mass_w, layer_index + 1, idx, _get_layer(mass_w, layer_index + 1, idx) + excess_water)
            end
        end
    end

    return runoff
end

"""
    go_percolation!(solid_mass, liquid_water_mass, snow_density, ice_density, water_density; ...)

Run the percolation scheme on one vector-backed snow column. Mutates
`liquid_water_mass` in-place and returns the runoff mass.
"""
function go_percolation!(
    solid_mass::AbstractVector,
    liquid_water_mass::AbstractVector,
    snow_density::AbstractVector,
    ice_density,
    water_density;
    max_lwc=oftype(water_density, 0.1),
    rho_i_tol=oftype(water_density, 10.0),
)
    N_ref = Ref(length(solid_mass))
    return _go_percolation!(
        N_ref,
        solid_mass,
        liquid_water_mass,
        snow_density,
        1,
        ice_density,
        water_density;
        max_lwc=max_lwc,
        rho_i_tol=rho_i_tol,
    )
end

"""
    go_percolation!(domain, idx; ...)

Run liquid-water percolation for column `idx` of `domain`. Mutates
`domain.mass_w` and accumulates routed runoff into `domain.runoff[idx]`.
"""
function go_percolation!(
    domain,
    idx::Int;
    max_lwc=oftype(domain.c.rho_w, 0.1),
    rho_i_tol=oftype(domain.c.rho_w, 10.0),
)
    if _n_active(domain.N, idx) <= 0 || _get_layer(domain.mass, 1, idx) <= zero(eltype(domain.mass))
        return zero(eltype(domain.mass))
    end

    runoff = _go_percolation!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        idx,
        domain.c.rho_i,
        domain.c.rho_w;
        max_lwc=max_lwc,
        rho_i_tol=rho_i_tol,
    )
    _set_scalar!(domain.runoff, idx, _get_scalar(domain.runoff, idx) + runoff)
    return runoff
end
