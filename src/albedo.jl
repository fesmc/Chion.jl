"""
Surface albedo state and update rules for both legacy constant and
BESSI-style dynamic albedo schemes.
"""

@inline function _constant_surface_albedo(column::SnowpackColumn)
    if column.N <= 0 || column.mass[1] <= EPS_EMPTY_LAYER
        return column.c.alpha_ice
    end
    return column.temperature[1] >= column.c.T0 ? column.c.alpha_wet : column.c.alpha_dry
end

@inline function _surface_liquid_water_content(column::SnowpackColumn)
    if column.N <= 0 || column.mass[1] <= EPS_TINY
        return 0.0
    end

    density = column.density[1]
    if density <= EPS_TINY || density >= column.c.rho_i - EPS_TINY
        return 0.0
    end

    pore_volume = column.mass[1] / density - column.mass[1] / column.c.rho_i
    if pore_volume <= EPS_TINY
        return 0.0
    end

    return max(column.mass_w[1], 0.0) / column.c.rho_w / pore_volume
end

function _refresh_dynamic_albedo_from_snowfall!(
    column::SnowpackColumn,
    snowfall_mass::Float64,
)
    snowfall_mass <= EPS_TINY && return column.albedo_dynamic

    if column.c.albedo_scheme == :constant
        column.albedo_dynamic = column.c.alpha_dry
        return column.albedo_dynamic
    end

    column.albedo_dynamic = min(
        column.c.alpha_dry,
        column.albedo_dynamic +
        (column.c.alpha_dry - column.c.alpha_wet) * (1.0 - exp(-snowfall_mass / 3.0)),
    )
    return column.albedo_dynamic
end

"""
    update_surface_albedo!(column::SnowpackColumn) -> Float64

Update the current surface albedo state.

- `:constant` uses the legacy dry-snow / wet-snow / ice switch.
- `:dynamic` uses the BESSI default dynamic scheme: snowfall refreshes the
  albedo toward fresh snow, and the Aoki-style temperature/liquid-water
  reduction darkens the snow surface between snowfall events.
"""
function update_surface_albedo!(column::SnowpackColumn)
    if column.N <= 0 || column.mass[1] <= EPS_EMPTY_LAYER
        column.albedo_dynamic = column.c.alpha_ice
        return column.albedo_dynamic
    end

    if column.c.albedo_scheme == :constant
        column.albedo_dynamic = _constant_surface_albedo(column)
        return column.albedo_dynamic
    end

    previous_albedo = clamp(column.albedo_dynamic, column.c.alpha_wet, column.c.alpha_dry)
    surface_temperature = column.temperature[1]

    # BESSI radiation.f90, albedo_module == 4 (Aoki-style reduction).
    updated_albedo = min(
        previous_albedo,
        previous_albedo - ((surface_temperature - column.c.T0) * 1.35e-3 + 0.0278),
    )
    updated_albedo = max(updated_albedo, column.c.alpha_wet)

    liquid_water_content = _surface_liquid_water_content(column)
    if liquid_water_content > 0.0 && column.c.max_lwc_albedo > EPS_TINY
        wet_adjusted_albedo = updated_albedo - (
            updated_albedo - column.c.alpha_wet
        ) * (liquid_water_content / column.c.max_lwc_albedo)
        updated_albedo = max(
            column.c.alpha_wet,
            min(updated_albedo, wet_adjusted_albedo),
        )
    end

    column.albedo_dynamic = clamp(updated_albedo, column.c.alpha_wet, column.c.alpha_dry)
    return column.albedo_dynamic
end
