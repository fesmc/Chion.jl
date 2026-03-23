"""
Column-level helper utilities shared across process modules and diagnostics.
"""

@inline function _resolve_keyword_alias(
    preferred_value,
    legacy_value,
    preferred_name::AbstractString,
    legacy_name::AbstractString,
)
    if !isnothing(preferred_value) && !isnothing(legacy_value) &&
       !isequal(preferred_value, legacy_value)
        error(
            "Received both `$preferred_name` and `$legacy_name` with different values. " *
            "Use one or provide matching values.",
        )
    end
    return isnothing(preferred_value) ? legacy_value : preferred_value
end

@inline function _bulk_snow_density(column::SnowpackColumn)
    if column.N <= 0
        return 0.0
    end

    total_mass = 0.0
    total_thickness = 0.0
    @inbounds for layer_index in 1:column.N
        layer_mass = column.mass[layer_index]
        layer_density = column.density[layer_index]
        if layer_mass > 0.0 && layer_density > EPS_TINY
            total_mass += layer_mass
            total_thickness += layer_mass / layer_density
        end
    end

    if total_mass <= 0.0 || total_thickness <= EPS_TINY
        return 0.0
    end
    return total_mass / total_thickness
end

@inline function _total_snow_water_mass(column::SnowpackColumn)
    if column.N <= 0
        return 0.0
    end

    total_wet_mass = 0.0
    @inbounds for layer_index in 1:column.N
        total_wet_mass += max(column.mass[layer_index], 0.0) +
                          max(column.mass_w[layer_index], 0.0)
    end
    return total_wet_mass
end

@inline function _snow_cover_fraction(column::SnowpackColumn)
    if column.N <= 0
        return 0.0
    end

    total_wet_mass = _total_snow_water_mass(column)
    total_wet_mass <= 0.0 && return 0.0

    bulk_density = _bulk_snow_density(column)
    bulk_density <= EPS_TINY && return 0.0

    return min(1.0, (total_wet_mass / bulk_density) / 0.1)
end

@inline function update_snow_cover!(column::SnowpackColumn)
    column.snow_cover = _snow_cover_fraction(column)
    return column.snow_cover
end

@inline function _column_has_liquid_water(column::SnowpackColumn)
    @inbounds for layer_index in 1:column.N
        if column.mass_w[layer_index] > EPS_TINY
            return true
        end
    end
    return false
end
