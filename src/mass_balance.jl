"""
Mass-balance and layer-structure helpers for array-backed snowpack states.
"""

@inline _safe_nonnegative(x) = x > zero(x) ? x : zero(x)

@inline function _mass_weighted_mean(m1, x1, m2, x2)
    total_mass = m1 + m2
    return total_mass > zero(total_mass) ? (m1 * x1 + m2 * x2) / total_mass : zero(x1 + x2)
end

@inline function _fresh_snow_density(
    c::SnowpackPhysicalConstants,
    air_temperature,
    wind_speed,
)
    if _uses_constant_fresh_snow_density(c)
        return clamp(c.rho_s, oftype(c.rho_s, 50), c.rho_i)
    end

    nonnegative_wind_speed = max(wind_speed, zero(wind_speed))
    fresh_snow_density = c.rho_s_a +
                         c.rho_s_b * (air_temperature - c.T0) +
                         c.rho_s_c * sqrt(nonnegative_wind_speed)
    return clamp(fresh_snow_density, oftype(fresh_snow_density, 50), c.rho_i)
end

function _reset_layer_at_index!(
    mass,
    mass_w,
    density,
    temperature,
    idx::Int,
    layer_index::Int,
    c::SnowpackPhysicalConstants,
)
    _set_layer!(mass, layer_index, idx, zero(eltype(mass)))
    _set_layer!(mass_w, layer_index, idx, zero(eltype(mass_w)))
    _set_layer!(density, layer_index, idx, zero(eltype(density)))
    _set_layer!(temperature, layer_index, idx, c.T0)
    return nothing
end

function _split_surface_layer!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    idx::Int,
    Ntot::Int,
    mass_max,
    mass_split,
)
    surface_mass = _get_layer(mass, 1, idx)
    if !(_n_active(N_storage, idx) < Ntot && surface_mass > mass_max)
        return nothing
    end

    surface_mass_w = _get_layer(mass_w, 1, idx)
    surface_density = _get_layer(density, 1, idx)
    surface_temperature = _get_layer(temperature, 1, idx)
    new_n = _n_active(N_storage, idx) + 1
    _set_n_active!(N_storage, idx, new_n)

    @inbounds for layer_index in new_n:-1:3
        _set_layer!(mass, layer_index, idx, _get_layer(mass, layer_index - 1, idx))
        _set_layer!(mass_w, layer_index, idx, _get_layer(mass_w, layer_index - 1, idx))
        _set_layer!(density, layer_index, idx, _get_layer(density, layer_index - 1, idx))
        _set_layer!(temperature, layer_index, idx, _get_layer(temperature, layer_index - 1, idx))
    end

    _set_layer!(mass, 2, idx, mass_split)
    _set_layer!(mass, 1, idx, surface_mass - mass_split)

    water_fraction = mass_split / surface_mass
    _set_layer!(mass_w, 2, idx, surface_mass_w * water_fraction)
    _set_layer!(mass_w, 1, idx, surface_mass_w * (one(surface_mass_w) - water_fraction))
    _set_layer!(density, 1, idx, surface_density)
    _set_layer!(density, 2, idx, surface_density)
    _set_layer!(temperature, 1, idx, surface_temperature)
    _set_layer!(temperature, 2, idx, surface_temperature)
    return nothing
end

function _merge_surface_layer!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    idx::Int,
    Ntot::Int,
    mass_split,
    mass_min,
    c::SnowpackPhysicalConstants,
)
    n = _n_active(N_storage, idx)
    surface_mass = n > 0 ? _get_layer(mass, 1, idx) : zero(eltype(mass))

    if n == 1 && surface_mass < EPS_EMPTY_LAYER
        _set_n_active!(N_storage, idx, 0)
        _reset_layer_at_index!(mass, mass_w, density, temperature, idx, 1, c)
        return nothing
    elseif n <= 1 || surface_mass >= mass_min
        return nothing
    end

    subsurface_mass = _get_layer(mass, 2, idx)
    combined_mass = surface_mass + subsurface_mass

    if combined_mass > oftype(combined_mass, 2) * mass_split
        transferred_to_surface = mass_split - surface_mass
        transferred_water = transferred_to_surface / subsurface_mass * _get_layer(mass_w, 2, idx)
        _set_layer!(mass, 1, idx, mass_split)
        _set_layer!(mass, 2, idx, combined_mass - mass_split)
        _set_layer!(mass_w, 1, idx, _get_layer(mass_w, 1, idx) + transferred_water)
        _set_layer!(mass_w, 2, idx, _get_layer(mass_w, 2, idx) - transferred_water)
        _set_layer!(
            density,
            1,
            idx,
            _mass_weighted_mean(surface_mass, _get_layer(density, 1, idx), transferred_to_surface, _get_layer(density, 2, idx)),
        )
        _set_layer!(
            temperature,
            1,
            idx,
            _mass_weighted_mean(surface_mass, _get_layer(temperature, 1, idx), transferred_to_surface, _get_layer(temperature, 2, idx)),
        )
        return nothing
    end

    _set_layer!(mass, 1, idx, combined_mass)
    _set_layer!(mass_w, 1, idx, _get_layer(mass_w, 1, idx) + _get_layer(mass_w, 2, idx))
    _set_layer!(
        density,
        1,
        idx,
        _mass_weighted_mean(surface_mass, _get_layer(density, 1, idx), subsurface_mass, _get_layer(density, 2, idx)),
    )
    _set_layer!(
        temperature,
        1,
        idx,
        _mass_weighted_mean(surface_mass, _get_layer(temperature, 1, idx), subsurface_mass, _get_layer(temperature, 2, idx)),
    )

    new_n = n - 1
    _set_n_active!(N_storage, idx, new_n)
    @inbounds for layer_index in 2:new_n
        _set_layer!(mass, layer_index, idx, _get_layer(mass, layer_index + 1, idx))
        _set_layer!(mass_w, layer_index, idx, _get_layer(mass_w, layer_index + 1, idx))
        _set_layer!(density, layer_index, idx, _get_layer(density, layer_index + 1, idx))
        _set_layer!(temperature, layer_index, idx, _get_layer(temperature, layer_index + 1, idx))
    end
    _reset_layer_at_index!(mass, mass_w, density, temperature, idx, new_n + 1, c)
    return nothing
end

function _merge_bottom_layer!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    n = _n_active(N_storage, idx)
    if n < 2
        return nothing
    end

    _set_n_active!(N_storage, idx, n - 1)
    N_new = n - 1

    combined_mass = _get_layer(mass, N_new, idx) + _get_layer(mass, n, idx)
    combined_mass_w = _get_layer(mass_w, N_new, idx) + _get_layer(mass_w, n, idx)
    combined_density = _mass_weighted_mean(
        _get_layer(mass, N_new, idx),
        _get_layer(density, N_new, idx),
        _get_layer(mass, n, idx),
        _get_layer(density, n, idx),
    )
    combined_temperature = _mass_weighted_mean(
        _get_layer(mass, N_new, idx),
        _get_layer(temperature, N_new, idx),
        _get_layer(mass, n, idx),
        _get_layer(temperature, n, idx),
    )

    if combined_density > c.rho_i
        mass_limited_to_ice_density = combined_mass * (c.rho_i / combined_density)
        exported_excess = combined_mass - mass_limited_to_ice_density
        _set_scalar!(mass_base, idx, _get_scalar(mass_base, idx) + exported_excess)
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + exported_excess)
        _set_layer!(mass, N_new, idx, mass_limited_to_ice_density)
        _set_layer!(mass_w, N_new, idx, combined_mass_w)
        _set_layer!(density, N_new, idx, c.rho_i)
        _set_layer!(temperature, N_new, idx, combined_temperature)
    else
        _set_layer!(mass, N_new, idx, combined_mass)
        _set_layer!(mass_w, N_new, idx, combined_mass_w)
        _set_layer!(density, N_new, idx, combined_density)
        _set_layer!(temperature, N_new, idx, combined_temperature)
    end

    _reset_layer_at_index!(mass, mass_w, density, temperature, idx, n, c)
    return nothing
end

function _remove_surface_layer!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    n = _n_active(N_storage, idx)
    if n <= 0
        return nothing
    elseif n == 1
        _reset_layer_at_index!(mass, mass_w, density, temperature, idx, 1, c)
        _set_n_active!(N_storage, idx, 0)
        return nothing
    end

    @inbounds for layer_index in 1:(n - 1)
        _set_layer!(mass, layer_index, idx, _get_layer(mass, layer_index + 1, idx))
        _set_layer!(mass_w, layer_index, idx, _get_layer(mass_w, layer_index + 1, idx))
        _set_layer!(density, layer_index, idx, _get_layer(density, layer_index + 1, idx))
        _set_layer!(temperature, layer_index, idx, _get_layer(temperature, layer_index + 1, idx))
    end
    _reset_layer_at_index!(mass, mass_w, density, temperature, idx, n, c)
    _set_n_active!(N_storage, idx, n - 1)
    return nothing
end

@inline function _remove_depleted_surface_and_route_water!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    runoff,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    n = _n_active(N_storage, idx)
    if n > 1
        _set_layer!(mass_w, 2, idx, _get_layer(mass_w, 2, idx) + _get_layer(mass_w, 1, idx))
    else
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + _get_layer(mass_w, 1, idx))
    end
    _set_layer!(mass_w, 1, idx, zero(eltype(mass_w)))
    _remove_surface_layer!(N_storage, mass, mass_w, density, temperature, idx, c)
    return nothing
end

function _continuous_bottom_deplete!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    Tsrf,
    albedo_dynamic,
    idx::Int,
    d_m_in,
    c::SnowpackPhysicalConstants,
)
    d_m = _safe_nonnegative(d_m_in)
    ice_to_base = zero(d_m)
    runoff_total = zero(d_m)

    while d_m > EPS_TINY && _n_active(N_storage, idx) > 0
        nn = _n_active(N_storage, idx)
        layer_mass = _get_layer(mass, nn, idx)

        if layer_mass <= EPS_EMPTY_LAYER
            _reset_layer_at_index!(mass, mass_w, density, temperature, idx, nn, c)
            _set_n_active!(N_storage, idx, nn - 1)
            continue
        end

        if d_m > layer_mass
            d_m -= layer_mass
            ice_to_base += layer_mass
            runoff_mass = _get_layer(mass_w, nn, idx)
            runoff_total += runoff_mass
            _set_scalar!(mass_base, idx, _get_scalar(mass_base, idx) + layer_mass)
            _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + layer_mass)
            _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + runoff_mass)
            _reset_layer_at_index!(mass, mass_w, density, temperature, idx, nn, c)
            _set_n_active!(N_storage, idx, nn - 1)
        else
            d_lw = d_m * _get_layer(mass_w, nn, idx) / layer_mass
            _set_layer!(mass, nn, idx, layer_mass - d_m)
            _set_layer!(mass_w, nn, idx, _get_layer(mass_w, nn, idx) - d_lw)
            ice_to_base += d_m
            runoff_total += d_lw
            _set_scalar!(mass_base, idx, _get_scalar(mass_base, idx) + d_m)
            _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + d_m)
            _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + d_lw)
            d_m = zero(d_m)
        end
    end

    while _n_active(N_storage, idx) > 0 && _get_layer(mass, _n_active(N_storage, idx), idx) <= EPS_EMPTY_LAYER
        tail = _n_active(N_storage, idx)
        _reset_layer_at_index!(mass, mass_w, density, temperature, idx, tail, c)
        _set_n_active!(N_storage, idx, tail - 1)
    end

    _set_scalar!(Tsrf, idx, _n_active(N_storage, idx) > 0 ? _get_layer(temperature, 1, idx) : c.T0)
    if _n_active(N_storage, idx) == 0
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return (ice_to_base=ice_to_base, runoff=runoff_total)
end

function _free_slot_for_surface_split!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    Tsrf,
    albedo_dynamic,
    idx::Int,
    Ntot::Int,
    mass_max,
    c::SnowpackPhysicalConstants,
)
    if Ntot == 1
        overflow = max(_get_layer(mass, 1, idx) - mass_max, zero(eltype(mass)))
        if overflow > zero(overflow)
            _continuous_bottom_deplete!(
                N_storage,
                mass,
                mass_w,
                density,
                temperature,
                mass_base,
                smb_ice,
                runoff,
                Tsrf,
                albedo_dynamic,
                idx,
                overflow,
                c,
            )
        end
        return nothing
    end

    bottom_mass = _get_layer(mass, _n_active(N_storage, idx), idx)
    if bottom_mass > zero(bottom_mass)
        _continuous_bottom_deplete!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            mass_base,
            smb_ice,
            runoff,
            Tsrf,
            albedo_dynamic,
            idx,
            bottom_mass,
            c,
        )
    else
        _reset_layer_at_index!(mass, mass_w, density, temperature, idx, _n_active(N_storage, idx), c)
        _set_n_active!(N_storage, idx, _n_active(N_storage, idx) - 1)
    end
    return nothing
end

function _enforce_mass_cap!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    Tsrf,
    albedo_dynamic,
    idx::Int,
    Ntot::Int,
    mass_split,
    f_base_max,
    dt_seconds,
    c::SnowpackPhysicalConstants,
)
    total_active_solid_mass = zero(eltype(mass))
    @inbounds for layer_index in 1:_n_active(N_storage, idx)
        total_active_solid_mass += _get_layer(mass, layer_index, idx)
    end

    if Ntot <= 3
        mass_cap = oftype(total_active_solid_mass, 1.5) * mass_split
        excess = total_active_solid_mass - mass_cap
        if excess > zero(excess)
            _continuous_bottom_deplete!(
                N_storage,
                mass,
                mass_w,
                density,
                temperature,
                mass_base,
                smb_ice,
                runoff,
                Tsrf,
                albedo_dynamic,
                idx,
                excess,
                c,
            )
        end
        return nothing
    end

    reference_column_mass_cap = BESSI_REFERENCE_LAYER_COUNT * mass_split * oftype(mass_split, 1.5)
    excess_basal_mass = total_active_solid_mass - reference_column_mass_cap
    if excess_basal_mass > zero(excess_basal_mass)
        relaxation = one(excess_basal_mass) - exp(-f_base_max * dt_seconds / c.seconds_per_day)
        _continuous_bottom_deplete!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            mass_base,
            smb_ice,
            runoff,
            Tsrf,
            albedo_dynamic,
            idx,
            excess_basal_mass * relaxation,
            c,
        )
    end
    return nothing
end

function _apply_accumulation!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    Tsrf,
    snow_cover,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    Ntot::Int,
    mass_max,
    mass_split,
    mass_min,
    f_base_max,
    snowfall_rate,
    rainfall_rate,
    dt_seconds;
    air_temperature=nothing,
    T_air=nothing,
    wind_speed=oftype(dt_seconds, 5),
)
    resolved_air_temperature = _resolve_keyword_alias(air_temperature, T_air, "air_temperature", "T_air")
    resolved_air_temperature = isnothing(resolved_air_temperature) ? c.T0 : resolved_air_temperature

    if _n_active(N_storage, idx) == 0
        if snowfall_rate > zero(snowfall_rate)
            _set_n_active!(N_storage, idx, 1)
        else
            _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
            return nothing
        end
    end

    if snowfall_rate > zero(snowfall_rate)
        previous_surface_mass = _get_layer(mass, 1, idx)
        added_snow_mass = snowfall_rate * dt_seconds
        updated_surface_mass = previous_surface_mass + added_snow_mass
        fresh_snow_density = _fresh_snow_density(c, resolved_air_temperature, wind_speed)
        previous_surface_density = _get_layer(density, 1, idx) > zero(_get_layer(density, 1, idx)) ?
            _get_layer(density, 1, idx) : fresh_snow_density
        if updated_surface_mass > zero(updated_surface_mass)
            updated_density = updated_surface_mass / (
                previous_surface_mass / previous_surface_density +
                added_snow_mass / fresh_snow_density
            )
            _set_layer!(density, 1, idx, updated_density)
        end
        _set_layer!(mass, 1, idx, updated_surface_mass)
        _refresh_dynamic_albedo_from_snowfall!(albedo_dynamic, idx, c, added_snow_mass)
    end

    if _get_layer(mass, 1, idx) > zero(eltype(mass)) && rainfall_rate > zero(rainfall_rate)
        _set_layer!(mass_w, 1, idx, _get_layer(mass_w, 1, idx) + rainfall_rate * dt_seconds)
    end

    while _n_active(N_storage, idx) > 0 && _get_layer(mass, 1, idx) > mass_max
        if _n_active(N_storage, idx) == Ntot
            if Ntot <= 2
                _free_slot_for_surface_split!(
                    N_storage,
                    mass,
                    mass_w,
                    density,
                    temperature,
                    mass_base,
                    smb_ice,
                    runoff,
                    Tsrf,
                    albedo_dynamic,
                    idx,
                    Ntot,
                    mass_max,
                    c,
                )
            else
                _merge_bottom_layer!(
                    N_storage,
                    mass,
                    mass_w,
                    density,
                    temperature,
                    mass_base,
                    smb_ice,
                    idx,
                    c,
                )
            end
            if _n_active(N_storage, idx) == 0 || _get_layer(mass, 1, idx) <= mass_max
                break
            end
        end
        _split_surface_layer!(N_storage, mass, mass_w, density, temperature, idx, Ntot, mass_max, mass_split)
    end

    while _n_active(N_storage, idx) > 1 && _get_layer(mass, 1, idx) < mass_min
        _merge_surface_layer!(N_storage, mass, mass_w, density, temperature, idx, Ntot, mass_split, mass_min, c)
    end

    _enforce_mass_cap!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        albedo_dynamic,
        idx,
        Ntot,
        mass_split,
        f_base_max,
        dt_seconds,
        c,
    )

    return nothing
end

function continuous_bottom_deplete!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    d_m_in,
)
    return _continuous_bottom_deplete!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.mass_base,
        domain.smb_ice,
        domain.runoff,
        domain.Tsrf,
        domain.albedo_dynamic,
        idx,
        d_m_in,
        domain.c,
    )
end

function apply_accumulation!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    snowfall_rate,
    rainfall_rate,
    dt_seconds;
    kwargs...,
)
    return _apply_accumulation!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.mass_base,
        domain.smb_ice,
        domain.runoff,
        domain.Tsrf,
        domain.snow_cover,
        domain.albedo_dynamic,
        idx,
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        domain.f_base_max,
        snowfall_rate,
        rainfall_rate,
        dt_seconds;
        kwargs...,
    )
end

function apply_melt!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    melt_mass,
)
    return _apply_melt!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.runoff,
        domain.Tsrf,
        domain.albedo_dynamic,
        idx,
        domain.mass_split,
        domain.mass_min,
        melt_mass,
        domain.c,
    )
end

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
