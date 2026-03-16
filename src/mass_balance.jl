"""
Mass-balance and layer-structure helpers for `SnowpackColumn`.
"""

@inline _safe_nonnegative(x::Float64) = x > 0.0 ? x : 0.0
@inline _mass_weighted_mean(m1::Float64, x1::Float64, m2::Float64, x2::Float64) =
    (m1 * x1 + m2 * x2) / (m1 + m2)
@inline function _fresh_snow_density(column::SnowpackColumn, T_air::Float64, wind_speed::Float64)
    c = column.c
    V = max(wind_speed, 0.0)
    ρfresh = c.rho_s_a + c.rho_s_b * (T_air - c.T0) + c.rho_s_c * sqrt(V)
    return clamp(ρfresh, 50.0, c.rho_i)
end

function reset_column_at_index!(column::SnowpackColumn, i::Int)
    column.mass[i] = 0.0
    column.mass_w[i] = 0.0
    column.density[i] = 0.0
    column.temperature[i] = column.c.T0
    return
end

"""
    split_surface_layer!(column::SnowpackColumn)

Split the surface layer when it exceeds `mass_max`.
Lower part contains `mass_split`, upper part contains remainder.
"""
function split_surface_layer!(column::SnowpackColumn)
    @assert column.mass[1] > column.mass_max
    @assert column.N < column.Ntot

    surface_mass = column.mass[1]
    surface_mass_w = column.mass_w[1]
    surface_density = column.density[1]
    surface_temperature = column.temperature[1]

    column.N += 1

    # Shift all active layers down, keeping indices 1 and 2 for split products.
    @inbounds for i in column.N:-1:3
        column.mass[i] = column.mass[i - 1]
        column.mass_w[i] = column.mass_w[i - 1]
        column.density[i] = column.density[i - 1]
        column.temperature[i] = column.temperature[i - 1]
    end

    column.mass[2] = column.mass_split
    column.mass[1] = surface_mass - column.mass_split

    water_fraction = column.mass_split / surface_mass
    column.mass_w[2] = surface_mass_w * water_fraction
    column.mass_w[1] = surface_mass_w * (1.0 - water_fraction)

    column.density[1] = surface_density
    column.density[2] = surface_density
    column.temperature[1] = surface_temperature
    column.temperature[2] = surface_temperature
    return
end

"""
    merge_surface_layer!(column::SnowpackColumn)

Merge surface layer with the second layer when `mass[1] < mass_min`.
If combined mass exceeds `2*mass_split`, only transfer enough mass to
restore the surface to `mass_split`.
"""
function merge_surface_layer!(column::SnowpackColumn)
    if column.N == 1 && column.mass[1] < EPS_EMPTY_LAYER
        column.N = 0
        reset_column_at_index!(column, 1)
        return
    elseif column.N == 1
        return
    end

    surface_mass = column.mass[1]
    subsurface_mass = column.mass[2]
    combined_mass = surface_mass + subsurface_mass

    if combined_mass > 2.0 * column.mass_split
        transferred_to_surface = column.mass_split - surface_mass
        transferred_water = transferred_to_surface / column.mass[2] * column.mass_w[2]

        column.mass[1] = column.mass_split
        column.mass[2] = combined_mass - column.mass_split
        column.mass_w[1] += transferred_water
        column.mass_w[2] -= transferred_water

        column.density[1] = _mass_weighted_mean(
            surface_mass,
            column.density[1],
            transferred_to_surface,
            column.density[2],
        )
        column.temperature[1] = _mass_weighted_mean(
            surface_mass,
            column.temperature[1],
            transferred_to_surface,
            column.temperature[2],
        )
    else
        column.mass[1] = combined_mass
        column.mass_w[1] += column.mass_w[2]
        column.density[1] = _mass_weighted_mean(
            surface_mass,
            column.density[1],
            subsurface_mass,
            column.density[2],
        )
        column.temperature[1] = _mass_weighted_mean(
            surface_mass,
            column.temperature[1],
            subsurface_mass,
            column.temperature[2],
        )

        column.N -= 1
        @inbounds for i in 2:column.N
            column.mass[i] = column.mass[i + 1]
            column.mass_w[i] = column.mass_w[i + 1]
            column.density[i] = column.density[i + 1]
            column.temperature[i] = column.temperature[i + 1]
        end
        reset_column_at_index!(column, column.N + 1)
    end

    return
end

"""
    merge_bottom_layer!(column::SnowpackColumn)

Merge the two lowest layers to make room for a new surface layer split.
"""
function merge_bottom_layer!(column::SnowpackColumn)
    @assert column.N == column.Ntot

    column.N -= 1
    N = column.N
    Np1 = column.N + 1

    combined_mass = column.mass[N] + column.mass[Np1]
    combined_mass_w = column.mass_w[N] + column.mass_w[Np1]
    combined_density = _mass_weighted_mean(
        column.mass[N],
        column.density[N],
        column.mass[Np1],
        column.density[Np1],
    )
    combined_temperature = _mass_weighted_mean(
        column.mass[N],
        column.temperature[N],
        column.mass[Np1],
        column.temperature[Np1],
    )

    if combined_density > column.c.rho_i
        mass_limited_to_ice_density = combined_mass * (column.c.rho_i / combined_density)
        column.mass_base += combined_mass - mass_limited_to_ice_density
        column.mass[N] = mass_limited_to_ice_density
        column.mass_w[N] = combined_mass_w
        column.density[N] = column.c.rho_i
        column.temperature[N] = combined_temperature
    else
        column.mass[N] = combined_mass
        column.mass_w[N] = combined_mass_w
        column.density[N] = combined_density
        column.temperature[N] = combined_temperature
    end

    reset_column_at_index!(column, Np1)
    return
end

"""
    apply_accumulation!(column::SnowpackColumn, P_snow::Float64, P_rain::Float64, dt_sec::Float64;
                        T_air::Float64=column.c.T0, wind_speed::Float64=5.0)

Apply snowfall/rainfall forcing and enforce dynamic layer bounds.
Fresh-snow density is parameterized as
`a + b*(T_air - T0) + c*sqrt(wind_speed)` using `column.c.rho_s_a/b/c`.
"""
function apply_accumulation!(
    column::SnowpackColumn,
    P_snow::Float64,
    P_rain::Float64,
    dt_sec::Float64,
    ;
    T_air::Float64=column.c.T0,
    wind_speed::Float64=5.0,
)
    # Fortran behavior: rain alone does not create a new snow layer.
    if column.N == 0
        if P_snow > 0.0
            column.N = 1
        else
            return
        end
    end

    if P_snow > 0.0
        old_mass = column.mass[1]
        add_snow = P_snow * dt_sec
        masssum = old_mass + add_snow
        fresh_rho = _fresh_snow_density(column, T_air, wind_speed)
        old_rho = column.density[1] > 0.0 ? column.density[1] : fresh_rho
        if masssum > 0.0
            column.density[1] = masssum / (old_mass / old_rho + add_snow / fresh_rho)
        end
        column.mass[1] = masssum
    end

    if column.mass[1] > 0.0 && P_rain > 0.0
        column.mass_w[1] += P_rain * dt_sec
    end

    while column.mass[1] > column.mass_max
        if column.N == column.Ntot
            merge_bottom_layer!(column)
        end
        split_surface_layer!(column)
    end

    while column.N > 1 && column.mass[1] < column.mass_min
        merge_surface_layer!(column)
    end

    if column.N == column.Ntot && column.mass[column.N] > column.mass_max
        mass_to_base = column.f_base_max * column.mass[column.N]
        column.mass[column.N] -= mass_to_base
        column.mass_base += mass_to_base
    end

    return
end

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

@inline function _remove_depleted_surface_and_route_water!(column::SnowpackColumn)
    if column.N > 1
        column.mass_w[2] += column.mass_w[1]
    else
        column.runoff += column.mass_w[1]
    end
    column.mass_w[1] = 0.0
    _remove_surface_layer!(column)
    return
end

"""
    apply_melt!(column::SnowpackColumn, melt_mass::Float64) -> Float64

Remove melt mass from the active snowpack (surface down), returning melted snow mass [kg m^-2].
Melted snow is converted to liquid water in-pack; runoff is only produced when
no receiving snow layer exists.
"""
function apply_melt!(column::SnowpackColumn, melt_mass::Float64)
    remaining_melt = _safe_nonnegative(melt_mass)
    if remaining_melt <= 0.0 || column.N <= 0
        return 0.0
    end

    melted_total = 0.0

    while remaining_melt > 0.0 && column.N > 0
        m_layer = column.mass[1]
        if m_layer <= EPS_TINY
            _remove_depleted_surface_and_route_water!(column)
            continue
        end

        dm = min(m_layer, remaining_melt)
        column.mass[1] -= dm
        column.mass_w[1] += dm
        remaining_melt -= dm
        melted_total += dm

        if column.mass[1] <= EPS_EMPTY_LAYER
            _remove_depleted_surface_and_route_water!(column)
        elseif column.N > 1 && column.mass[1] < column.mass_min
            merge_surface_layer!(column)
        end
    end

    column.Tsrf = column.N > 0 ? column.temperature[1] : column.c.T0
    return melted_total
end
