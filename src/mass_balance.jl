"""
Mass-balance helpers for `SnowpackColumn`.
"""

@inline _safe_nonnegative(x::Float64) = x > 0.0 ? x : 0.0

function _remove_surface_layer!(column::SnowpackColumn)
    if column.N <= 0
        return
    elseif column.N == 1
        reset_column_at_index!(column, 1)
        column.N = 0
        return
    end

    @inbounds for i in 1:(column.N - 1)
        column.mass[i] = column.mass[i + 1]
        column.mass_w[i] = column.mass_w[i + 1]
        column.density[i] = column.density[i + 1]
        column.temperature[i] = column.temperature[i + 1]
    end

    reset_column_at_index!(column, column.N)
    column.N -= 1
    return
end

"""
    apply_melt!(column::SnowpackColumn, melt_mass::Float64) -> Float64

Remove melt mass from the active snowpack (surface down), returning melted snow mass [kg m^-2].
Melted snow is converted to liquid water in-pack; runoff is only produced when no receiving snow layer exists.
"""
function apply_melt!(column::SnowpackColumn, melt_mass::Float64)
    remaining_melt = _safe_nonnegative(melt_mass)
    if remaining_melt <= 0.0 || column.N <= 0
        return 0.0
    end

    melted_total = 0.0
    runoff_from_melt = 0.0

    while remaining_melt > 0.0 && column.N > 0
        m_layer = column.mass[1]
        if m_layer <= 1.0e-12
            if column.N > 1
                column.mass_w[2] += column.mass_w[1]
                column.mass_w[1] = 0.0
            else
                runoff_from_melt += column.mass_w[1]
            end
            _remove_surface_layer!(column)
            continue
        end

        dm = min(m_layer, remaining_melt)
        column.mass[1] -= dm
        column.mass_w[1] += dm

        remaining_melt -= dm
        melted_total += dm

        if column.mass[1] <= 1.0e-10
            if column.N > 1
                column.mass_w[2] += column.mass_w[1]
                column.mass_w[1] = 0.0
            else
                runoff_from_melt += column.mass_w[1]
            end
            _remove_surface_layer!(column)
        elseif column.N > 1 && column.mass[1] < column.mass_min
            merge_surface_layer!(column)
        end
    end

    column.runoff += runoff_from_melt
    column.Tsrf = column.N > 0 ? column.temperature[1] : column.c.T0
    return melted_total
end
