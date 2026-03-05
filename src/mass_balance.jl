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

Remove melt mass from the active snowpack (surface down), returning melted mass [kg m^-2].
Melted mass is added to `column.runoff`.
"""
function apply_melt!(column::SnowpackColumn, melt_mass::Float64)
    remaining_melt = _safe_nonnegative(melt_mass)
    if remaining_melt <= 0.0 || column.N <= 0
        return 0.0
    end

    melted_total = 0.0

    while remaining_melt > 0.0 && column.N > 0
        m_layer = column.mass[1]
        if m_layer <= 1.0e-12
            _remove_surface_layer!(column)
            continue
        end

        dm = min(m_layer, remaining_melt)
        wfrac = clamp(column.mass_w[1] / m_layer, 0.0, 1.0)
        column.mass[1] -= dm
        column.mass_w[1] = _safe_nonnegative(column.mass_w[1] - dm * wfrac)

        remaining_melt -= dm
        melted_total += dm

        if column.mass[1] <= 1.0e-10
            _remove_surface_layer!(column)
        elseif column.N > 1 && column.mass[1] < column.mass_min
            merge_surface_layer!(column)
        else
            break
        end
    end

    column.runoff += melted_total
    column.Tsrf = column.N > 0 ? column.temperature[1] : column.c.T0
    return melted_total
end
