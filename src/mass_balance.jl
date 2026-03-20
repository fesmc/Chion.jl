"""
Mass-balance and layer-structure helpers for `SnowpackColumn`.
"""

@inline _safe_nonnegative(x::Float64) = x > 0.0 ? x : 0.0
@inline _mass_weighted_mean(m1::Float64, x1::Float64, m2::Float64, x2::Float64) =
    (m1 * x1 + m2 * x2) / (m1 + m2)
@inline function _fresh_snow_density(
    column::SnowpackColumn,
    air_temperature::Float64,
    wind_speed::Float64,
)
    c = column.c
    if c.fresh_snow_density_scheme == :constant
        return clamp(c.rho_s, 50.0, c.rho_i)
    end

    nonnegative_wind_speed = max(wind_speed, 0.0)
    fresh_snow_density = c.rho_s_a +
                         c.rho_s_b * (air_temperature - c.T0) +
                         c.rho_s_c * sqrt(nonnegative_wind_speed)
    return clamp(fresh_snow_density, 50.0, c.rho_i)
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

function _free_slot_for_surface_split!(column::SnowpackColumn)
    @assert column.N == column.Ntot

    if column.Ntot == 1
        overflow = max(column.mass[1] - column.mass_max, 0.0)
        if overflow > 0.0
            continuous_bottom_deplete!(column, overflow)
        end
        return
    end

    # When the column is already full, retire the deepest active layer to the
    # base so a new surface split can occur without cycling merge/split forever.
    bottom_mass = column.mass[column.N]
    if bottom_mass > 0.0
        continuous_bottom_deplete!(column, bottom_mass)
    else
        reset_column_at_index!(column, column.N)
        column.N -= 1
    end
    return
end

"""
    apply_accumulation!(
        column::SnowpackColumn,
        snowfall_rate::Float64,
        rainfall_rate::Float64,
        dt_seconds::Float64;
        air_temperature::Float64=column.c.T0,
        wind_speed::Float64=5.0,
    )

Apply snowfall/rainfall forcing and enforce dynamic layer bounds.
Fresh-snow density is parameterized as
`a + b*(air_temperature - T0) + c*sqrt(wind_speed)` using `column.c.rho_s_a/b/c`.

Legacy keyword alias `T_air` is still accepted.
"""
function apply_accumulation!(
    column::SnowpackColumn,
    snowfall_rate::Float64,
    rainfall_rate::Float64,
    dt_seconds::Float64,
    ;
    air_temperature::Union{Nothing, Float64}=nothing,
    T_air::Union{Nothing, Float64}=nothing,
    wind_speed::Float64=5.0,
)
    resolved_air_temperature = _resolve_keyword_alias(air_temperature, T_air, "air_temperature", "T_air")
    resolved_air_temperature = isnothing(resolved_air_temperature) ? column.c.T0 : resolved_air_temperature

    # Fortran behavior: rain alone does not create a new snow layer.
    if column.N == 0
        if snowfall_rate > 0.0
            column.N = 1
        else
            return
        end
    end

    if snowfall_rate > 0.0
        previous_surface_mass = column.mass[1]
        added_snow_mass = snowfall_rate * dt_seconds
        updated_surface_mass = previous_surface_mass + added_snow_mass
        fresh_snow_density = _fresh_snow_density(column, resolved_air_temperature, wind_speed)
        previous_surface_density = column.density[1] > 0.0 ? column.density[1] : fresh_snow_density
        if updated_surface_mass > 0.0
            column.density[1] = updated_surface_mass / (
                previous_surface_mass / previous_surface_density +
                added_snow_mass / fresh_snow_density
            )
        end
        column.mass[1] = updated_surface_mass
    end

    if column.mass[1] > 0.0 && rainfall_rate > 0.0
        column.mass_w[1] += rainfall_rate * dt_seconds
    end

    while column.mass[1] > column.mass_max
        if column.N == column.Ntot
            _free_slot_for_surface_split!(column)
            column.N == 0 && break
        end
        split_surface_layer!(column)
    end

    while column.N > 1 && column.mass[1] < column.mass_min
        merge_surface_layer!(column)
    end
    total_active_solid_mass = column.N > 0 ? sum(@view column.mass[1:column.N]) : 0.0
    reference_column_mass_cap = BESSI_REFERENCE_LAYER_COUNT * column.mass_split * 1.5
    excess_basal_mass = total_active_solid_mass - reference_column_mass_cap
    if excess_basal_mass > 0.0
        continuous_bottom_deplete!(column, excess_basal_mass)
    end



    return
end

"""
    continuous_bottom_deplete!(column::SnowpackColumn, d_m_in::Float64) -> NamedTuple

Julia translation of the Fortran `continous_snowman_depleet` routine.
Remove `d_m_in` solid mass [kg m^-2] from the bottom upward, adding the solid
part to `mass_base` and routing removed liquid water to `runoff`.
"""
function continuous_bottom_deplete!(column::SnowpackColumn, d_m_in::Float64)
    d_m = _safe_nonnegative(d_m_in)
    ice_to_base = 0.0
    runoff = 0.0

    while d_m > EPS_TINY && column.N > 0
        nn = column.N
        m = column.mass[nn]

        if m <= EPS_EMPTY_LAYER
            reset_column_at_index!(column, nn)
            column.N -= 1
            continue
        end

        if d_m > m
            d_m -= m
            ice_to_base += m
            runoff += column.mass_w[nn]

            column.mass_base += m
            column.runoff += column.mass_w[nn]

            reset_column_at_index!(column, nn)
            column.N -= 1
        else
            d_lw = d_m * column.mass_w[nn] / m
            column.mass[nn] -= d_m
            column.mass_w[nn] -= d_lw

            ice_to_base += d_m
            runoff += d_lw

            column.mass_base += d_m
            column.runoff += d_lw
            d_m = 0.0
        end
    end

    while column.N > 0 && column.mass[column.N] <= EPS_EMPTY_LAYER
        reset_column_at_index!(column, column.N)
        column.N -= 1
    end

    column.Tsrf = column.N > 0 ? column.temperature[1] : column.c.T0
    return (ice_to_base = ice_to_base, runoff = runoff)
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
