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

function _apply_htessel_liquid_water_compaction!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    liquid_water_before_energy::AbstractVector,
    dt_seconds,
)
    return _apply_htessel_liquid_water_compaction!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        idx,
        liquid_water_before_energy,
        domain.c.rho_i,
    )
end

function _go_densification!(
    N_storage,
    mass,
    density,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
    accumulation_rate,
    dt_seconds;
    hl::Bool=false,
    rho_e=oftype(c.rho_i, 815.0),
    P_atm=oftype(c.rho_i, 101325.0),
)
    n = _n_active(N_storage, idx)
    if n <= 0
        return nothing
    end

    ice_density = c.rho_i
    solid_mass_above = zero(ice_density)
    @inbounds for layer_index in 1:n
        layer_density = _get_layer(density, layer_index, idx)
        layer_temperature = _get_layer(temperature, layer_index, idx)
        layer_solid_mass = _get_layer(mass, layer_index, idx)

        if !isfinite(layer_density) || !isfinite(layer_temperature)
            solid_mass_above += max(layer_solid_mass, zero(layer_solid_mass))
            continue
        elseif layer_density >= ice_density - EPS_TINY
            _set_layer!(density, layer_index, idx, ice_density)
            solid_mass_above += max(layer_solid_mass, zero(layer_solid_mass))
            continue
        end

        density_tendency = zero(layer_density)
        overburden_pressure = oftype(layer_density, 9.81) * (solid_mass_above + layer_solid_mass / oftype(layer_density, 2))

        if layer_density < oftype(layer_density, 550.0)
            if layer_solid_mass > zero(layer_solid_mass)
                density_tendency = if c.low_density_densification == :htessel
                    _htessel_low_density_rate(c, overburden_pressure, layer_temperature, layer_density)
                else
                    _bessi_low_density_rate(layer_density, layer_temperature, ice_density, accumulation_rate)
                end
            end
        elseif hl
            if layer_index > 1 && layer_solid_mass > zero(layer_solid_mass)
                density_tendency = oftype(layer_density, 0.575) * exp(-oftype(layer_density, 21400.0) / oftype(layer_density, 8.13) / layer_temperature) *
                                   (ice_density - layer_density) *
                                   (oftype(layer_density, 1000.0) / oftype(layer_density, 3600.0) / oftype(layer_density, 24.0) / oftype(layer_density, 365.0))^oftype(layer_density, 0.5) *
                                   max(accumulation_rate, zero(accumulation_rate))^oftype(layer_density, 0.5)
            end
        elseif layer_density < oftype(layer_density, 800.0)
            density_ratio = layer_density / ice_density
            densification_shape_factor = oftype(layer_density, 10.0)^(
                -oftype(layer_density, 29.166) * density_ratio^3 +
                oftype(layer_density, 84.422) * density_ratio^2 -
                oftype(layer_density, 87.425) * density_ratio +
                oftype(layer_density, 30.673)
            )
            ice_pressure_mpa = overburden_pressure / oftype(layer_density, 1.0e6)
            bubble_pressure_mpa = zero(layer_density)
            if layer_density > rho_e
                bubble_pressure_mpa = P_atm * (
                    (one(layer_density) / rho_e - one(layer_density) / ice_density) /
                    (one(layer_density) / layer_density - one(layer_density) / ice_density) - one(layer_density)
                ) / oftype(layer_density, 1.0e6)
            end
            pressure_excess_mpa = ice_pressure_mpa - bubble_pressure_mpa
            density_tendency = oftype(layer_density, 25400.0) * exp(-oftype(layer_density, 60000.0) / oftype(layer_density, 8.13) / layer_temperature) *
                               layer_density * densification_shape_factor * pressure_excess_mpa^3
        else
            relative_porosity = _relative_porosity(layer_density, ice_density)
            denominator = one(layer_density) - relative_porosity^(one(layer_density) / oftype(layer_density, 3))
            if abs(denominator) > EPS_TINY
                densification_shape_factor = oftype(layer_density, 3) / oftype(layer_density, 16) *
                                             relative_porosity / denominator^3
                ice_pressure_mpa = overburden_pressure / oftype(layer_density, 1.0e6)
                bubble_pressure_mpa = zero(layer_density)
                if layer_density > rho_e
                    bubble_pressure_mpa = P_atm * (
                        (one(layer_density) / rho_e - one(layer_density) / ice_density) /
                        (one(layer_density) / layer_density - one(layer_density) / ice_density) - one(layer_density)
                    ) / oftype(layer_density, 1.0e6)
                end
                pressure_excess_mpa = ice_pressure_mpa - bubble_pressure_mpa
                density_tendency = oftype(layer_density, 25400.0) * exp(-oftype(layer_density, 60000.0) / oftype(layer_density, 8.13) / layer_temperature) *
                                   layer_density * densification_shape_factor * pressure_excess_mpa^3
            end
        end

        updated_density = max(layer_density, layer_density + density_tendency * dt_seconds)
        _set_layer!(density, layer_index, idx, min(updated_density, ice_density))
        solid_mass_above += max(layer_solid_mass, zero(layer_solid_mass))
    end

    return nothing
end

function go_densification!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    accumulation_rate,
    dt_seconds;
    kwargs...,
)
    return _go_densification!(
        domain.N,
        domain.mass,
        domain.density,
        domain.temperature,
        idx,
        domain.c,
        accumulation_rate,
        dt_seconds;
        kwargs...,
    )
end
