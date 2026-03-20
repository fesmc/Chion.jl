"""
Firn densification translated from BESSI `go_densification`.
"""

@inline _bessi_low_density_rate(
    density::Float64,
    temperature::Float64,
    ice_density::Float64,
    accumulation_rate::Float64,
) = 0.011 * exp(-10160.0 / 8.13 / temperature) *
    (ice_density - density) * max(accumulation_rate, 0.0)

@inline _htessel_snow_viscosity(
    melting_temperature::Float64,
    temperature::Float64,
    density::Float64,
) = 3.7e7 * exp(8.1e-2 * (melting_temperature - temperature) + 1.8e-2 * density)

@inline _htessel_thermal_metamorphism(
    melting_temperature::Float64,
    temperature::Float64,
    density::Float64,
) = 2.8e-6 * exp(
    -4.2e-2 * (melting_temperature - temperature) - 460.0 * max(0.0, density - 150.0),
)

@inline function _htessel_low_density_rate(
    column::SnowpackColumn,
    overburden_pressure::Float64,
    temperature::Float64,
    density::Float64,
)
    snow_viscosity = _htessel_snow_viscosity(column.c.T0, temperature, density)
    thermal_metamorphism = _htessel_thermal_metamorphism(column.c.T0, temperature, density)
    return density * (overburden_pressure / snow_viscosity + thermal_metamorphism)
end

@inline _relative_porosity(density::Float64, ice_density::Float64) =
    clamp(1.0 - density / ice_density, 0.0, 1.0)

function _apply_htessel_liquid_water_compaction!(
    column::SnowpackColumn,
    liquid_water_before_energy::AbstractVector{Float64},
    dt_seconds::Float64,
)
    n = column.N
    if n <= 0 || dt_seconds <= 0.0
        return
    end

    @views solid_mass = column.mass[1:n]
    @views density = column.density[1:n]
    @views liquid_water_mass = column.mass_w[1:n]
    ice_density = column.c.rho_i

    for layer_index in 1:n
        layer_density = density[layer_index]
        layer_solid_mass = solid_mass[layer_index]
        if layer_density < 550.0 && layer_solid_mass > EPS_TINY
            previous_liquid_water_mass = layer_index <= length(liquid_water_before_energy) ?
                liquid_water_before_energy[layer_index] : 0.0
            retained_liquid_water_gain = max(
                liquid_water_mass[layer_index] - previous_liquid_water_mass,
                0.0,
            )
            if retained_liquid_water_gain > 0.0
                # HTESSEL eq. (8): retained meltwater increases density at fixed solid mass.
                density[layer_index] = min(
                    max(
                        layer_density,
                        layer_density + layer_density * retained_liquid_water_gain / layer_solid_mass,
                    ),
                    ice_density,
                )
            end
        end
    end
    return
end

"""
    go_densification!(
        column::SnowpackColumn,
        accumulation_rate::Float64,
        dt_seconds::Float64;
        hl::Bool=false,
        rho_e::Float64=815.0,
        P_atm::Float64=101325.0,
    )

Update layer densities in-place using the multi-regime densification scheme
from the original Fortran implementation.

- `accumulation_rate` is the accumulation proxy used by the HL branch.
- `dt_seconds` is timestep in seconds.
- `hl=true` forces HL for densities above 550 kg m^-3.
"""
function go_densification!(
    column::SnowpackColumn,
    accumulation_rate::Float64,
    dt_seconds::Float64;
    hl::Bool=false,
    rho_e::Float64=815.0,
    P_atm::Float64=101325.0,
)
    n = column.N
    if n <= 0
        return
    end

    ice_density = column.c.rho_i
    @views solid_mass = column.mass[1:n]
    @views temperature = column.temperature[1:n]
    @views density = column.density[1:n]
    low_density_scheme = column.c.low_density_densification

    solid_mass_above = 0.0
    for layer_index in 1:n
        layer_density = density[layer_index]
        layer_temperature = temperature[layer_index]
        layer_solid_mass = solid_mass[layer_index]
        density_tendency = 0.0
        overburden_pressure = 9.81 * (solid_mass_above + layer_solid_mass / 2.0)

        if layer_density < 550.0
            if layer_solid_mass > 0.0
                if low_density_scheme == :htessel
                    density_tendency = _htessel_low_density_rate(
                        column,
                        overburden_pressure,
                        layer_temperature,
                        layer_density,
                    )
                else
                    density_tendency = _bessi_low_density_rate(
                        layer_density,
                        layer_temperature,
                        ice_density,
                        accumulation_rate,
                    )
                end
            end
        elseif hl
            if layer_index > 1 && layer_solid_mass > 0.0
                density_tendency = 0.575 * exp(-21400.0 / 8.13 / layer_temperature) *
                                   (ice_density - layer_density) *
                                   (1000.0 / 3600.0 / 24.0 / 365.0)^0.5 *
                                   max(accumulation_rate, 0.0)^0.5
            end
        elseif layer_density < 800.0
            density_ratio = layer_density / ice_density
            densification_shape_factor = 10.0^(
                -29.166 * density_ratio^3 +
                84.422 * density_ratio^2 -
                87.425 * density_ratio +
                30.673
            )
            ice_pressure_mpa = overburden_pressure / 1.0e6
            bubble_pressure_mpa = 0.0
            if layer_density > rho_e
                bubble_pressure_mpa = P_atm * (
                    (1.0 / rho_e - 1.0 / ice_density) /
                    (1.0 / layer_density - 1.0 / ice_density) - 1.0
                ) / 1.0e6
            end
            pressure_excess_mpa = ice_pressure_mpa - bubble_pressure_mpa
            density_tendency = 25400.0 * exp(-60000.0 / 8.13 / layer_temperature) *
                               layer_density * densification_shape_factor * pressure_excess_mpa^3
        else
            relative_porosity = _relative_porosity(layer_density, ice_density)
            denominator = 1.0 - relative_porosity^(1.0 / 3.0)
            if abs(denominator) > EPS_TINY
                densification_shape_factor = 3.0 / 16.0 *
                                             relative_porosity /
                                             denominator^3
                ice_pressure_mpa = overburden_pressure / 1.0e6
                bubble_pressure_mpa = 0.0
                if layer_density > rho_e
                    bubble_pressure_mpa = P_atm * (
                        (1.0 / rho_e - 1.0 / ice_density) /
                        (1.0 / layer_density - 1.0 / ice_density) - 1.0
                    ) / 1.0e6
                end
                pressure_excess_mpa = ice_pressure_mpa - bubble_pressure_mpa
                density_tendency = 25400.0 * exp(-60000.0 / 8.13 / layer_temperature) *
                                   layer_density * densification_shape_factor * pressure_excess_mpa^3
            end
        end

        updated_density = max(layer_density, layer_density + density_tendency * dt_seconds)
        density[layer_index] = min(updated_density, ice_density)
        solid_mass_above += max(layer_solid_mass, 0.0)
    end

    return
end
