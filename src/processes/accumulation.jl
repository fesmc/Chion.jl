"""
Snowfall and rainfall accumulation for array-backed snowpack states.
"""

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
