"""
Firn densification translated into array-backed kernels.
"""

@inline _bessi_low_density_rate(
    density,
    temperature,
    ice_density,
    accumulation_rate,
) = oftype(density, 0.011) * exp(-oftype(temperature, 10160.0) / oftype(temperature, 8.13) / temperature) *
    (ice_density - density) * max(accumulation_rate, zero(accumulation_rate))

@inline _htessel_snow_viscosity(
    melting_temperature,
    temperature,
    density,
) = oftype(density, 3.7e7) * exp(oftype(density, 8.1e-2) * (melting_temperature - temperature) + oftype(density, 1.8e-2) * density)

@inline _htessel_thermal_metamorphism(
    melting_temperature,
    temperature,
    density,
) = oftype(density, 2.8e-6) * exp(
    -oftype(density, 4.2e-2) * (melting_temperature - temperature) - oftype(density, 460.0) * max(zero(density), density - oftype(density, 150.0)),
)

@inline function _htessel_low_density_rate(
    c::SnowpackPhysicalConstants,
    overburden_pressure,
    temperature,
    density,
)
    snow_viscosity = _htessel_snow_viscosity(c.T0, temperature, density)
    thermal_metamorphism = _htessel_thermal_metamorphism(c.T0, temperature, density)
    return density * (overburden_pressure / snow_viscosity + thermal_metamorphism)
end

@inline _relative_porosity(density, ice_density) = clamp(one(density) - density / ice_density, zero(density), one(density))

function _apply_htessel_liquid_water_compaction!(
    N_storage,
    mass,
    mass_w,
    density,
    idx::Int,
    liquid_water_before_energy::AbstractVector,
    ice_density,
)
    n = _n_active(N_storage, idx)
    if n <= 0
        return nothing
    end

    @inbounds for layer_index in 1:n
        layer_density = _get_layer(density, layer_index, idx)
        layer_solid_mass = _get_layer(mass, layer_index, idx)
        if layer_density < oftype(layer_density, 550) && layer_solid_mass > EPS_TINY
            previous_liquid_water_mass = layer_index <= length(liquid_water_before_energy) ?
                liquid_water_before_energy[layer_index] : zero(layer_density)
            retained_liquid_water_gain = max(_get_layer(mass_w, layer_index, idx) - previous_liquid_water_mass, zero(layer_density))
            if retained_liquid_water_gain > zero(retained_liquid_water_gain)
                updated_density = min(
                    max(layer_density, layer_density + layer_density * retained_liquid_water_gain / layer_solid_mass),
                    ice_density,
                )
                _set_layer!(density, layer_index, idx, updated_density)
            end
        end
    end
    return nothing
end

@inline function _bubble_pressure_mpa(density, rho_e, ice_density, P_atm)
    if density <= rho_e
        return zero(density)
    end
    return P_atm * (
        (one(density) / rho_e - one(density) / ice_density) /
        (one(density) / density - one(density) / ice_density) - one(density)
    ) / oftype(density, 1.0e6)
end

@inline function _low_density_tendency(
    ::Val{:bessi},
    c::SnowpackPhysicalConstants,
    overburden_pressure,
    temperature,
    density,
    accumulation_rate,
)
    return _bessi_low_density_rate(density, temperature, c.rho_i, accumulation_rate)
end

@inline function _low_density_tendency(
    ::Val{:htessel},
    c::SnowpackPhysicalConstants,
    overburden_pressure,
    temperature,
    density,
    accumulation_rate,
)
    return _htessel_low_density_rate(c, overburden_pressure, temperature, density)
end

@inline function _mid_density_tendency(
    density,
    temperature,
    pressure_excess_mpa,
    ice_density,
)
    density_ratio = density / ice_density
    densification_shape_factor = oftype(density, 10.0)^(
        -oftype(density, 29.166) * density_ratio^3 +
        oftype(density, 84.422) * density_ratio^2 -
        oftype(density, 87.425) * density_ratio +
        oftype(density, 30.673)
    )
    return oftype(density, 25400.0) * exp(-oftype(density, 60000.0) / oftype(density, 8.13) / temperature) *
           density * densification_shape_factor * pressure_excess_mpa^3
end

@inline function _high_density_tendency(
    density,
    temperature,
    pressure_excess_mpa,
    ice_density,
)
    relative_porosity = _relative_porosity(density, ice_density)
    denominator = one(density) - relative_porosity^(one(density) / oftype(density, 3))
    if abs(denominator) <= EPS_TINY
        return zero(density)
    end

    densification_shape_factor = oftype(density, 3) / oftype(density, 16) *
                                 relative_porosity / denominator^3
    return oftype(density, 25400.0) * exp(-oftype(density, 60000.0) / oftype(density, 8.13) / temperature) *
           density * densification_shape_factor * pressure_excess_mpa^3
end

function _go_densification_scheme!(
    N_storage,
    mass,
    density,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
    accumulation_rate,
    dt_seconds,
    low_density_scheme::Val{scheme},
) where {scheme}
    n = _n_active(N_storage, idx)
    if n <= 0
        return nothing
    end

    ice_density = c.rho_i
    rho_e = oftype(ice_density, 815.0)
    P_atm = oftype(ice_density, 101325.0)
    solid_mass_above = zero(ice_density)

    @inbounds for layer_index in 1:n
        layer_density = _get_layer(density, layer_index, idx)
        layer_temperature = _get_layer(temperature, layer_index, idx)
        layer_solid_mass = _get_layer(mass, layer_index, idx)
        nonnegative_solid_mass = max(layer_solid_mass, zero(layer_solid_mass))

        if !isfinite(layer_density) || !isfinite(layer_temperature)
            solid_mass_above += nonnegative_solid_mass
            continue
        end

        if layer_density >= ice_density - EPS_TINY
            _set_layer!(density, layer_index, idx, ice_density)
            solid_mass_above += nonnegative_solid_mass
            continue
        end

        density_tendency = zero(layer_density)
        if layer_solid_mass > zero(layer_solid_mass)
            overburden_pressure = oftype(layer_density, 9.81) *
                                  (solid_mass_above + layer_solid_mass / oftype(layer_density, 2))

            if layer_density < oftype(layer_density, 550.0)
                density_tendency = _low_density_tendency(
                    low_density_scheme,
                    c,
                    overburden_pressure,
                    layer_temperature,
                    layer_density,
                    accumulation_rate,
                )
            else
                ice_pressure_mpa = overburden_pressure / oftype(layer_density, 1.0e6)
                bubble_pressure_mpa = _bubble_pressure_mpa(layer_density, rho_e, ice_density, P_atm)
                pressure_excess_mpa = ice_pressure_mpa - bubble_pressure_mpa
                density_tendency = if layer_density < oftype(layer_density, 800.0)
                    _mid_density_tendency(layer_density, layer_temperature, pressure_excess_mpa, ice_density)
                else
                    _high_density_tendency(layer_density, layer_temperature, pressure_excess_mpa, ice_density)
                end
            end
        end

        updated_density = max(layer_density, layer_density + density_tendency * dt_seconds)
        _set_layer!(density, layer_index, idx, min(updated_density, ice_density))
        solid_mass_above += nonnegative_solid_mass
    end

    return nothing
end

function _go_densification!(
    N_storage,
    mass,
    density,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
    accumulation_rate,
    dt_seconds,
)
    if _uses_htessel_densification(c)
        return _go_densification_scheme!(
            N_storage,
            mass,
            density,
            temperature,
            idx,
            c,
            accumulation_rate,
            dt_seconds,
            Val(:htessel),
        )
    end
    return _go_densification_scheme!(
        N_storage,
        mass,
        density,
        temperature,
        idx,
        c,
        accumulation_rate,
        dt_seconds,
        Val(:bessi),
    )
end

function go_densification!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    accumulation_rate,
    dt_seconds,
)
    return _go_densification!(
        domain.N,
        domain.mass,
        domain.density,
        domain.temperature,
        idx,
        domain.c,
        accumulation_rate,
        dt_seconds,
    )
end
