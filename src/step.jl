"""
Top-level model step orchestration.
"""

@inline function _default_snow_fraction(column::SnowpackColumn, air_temperature::Float64)
    return air_temperature > column.c.T0 ? 0.0 : 1.0
end

function _resolve_step_forcing(
    column::SnowpackColumn,
    air_temperature::Float64,
    precipitation_rate::Float64;
    snow_fraction,
    f_s,
    snowfall_rate::Union{Nothing, Float64},
    rainfall_rate::Union{Nothing, Float64},
    shortwave_down::Union{Nothing, Float64},
    p_snow::Union{Nothing, Float64},
    p_rain::Union{Nothing, Float64},
    s_boa::Union{Nothing, Float64},
)
    resolved_snow_fraction = _resolve_keyword_alias(snow_fraction, f_s, "snow_fraction", "f_s")
    resolved_snowfall_rate = _resolve_keyword_alias(snowfall_rate, p_snow, "snowfall_rate", "p_snow")
    resolved_rainfall_rate = _resolve_keyword_alias(rainfall_rate, p_rain, "rainfall_rate", "p_rain")
    resolved_shortwave_down = _resolve_keyword_alias(shortwave_down, s_boa, "shortwave_down", "s_boa")

    if !isnothing(resolved_snowfall_rate) || !isnothing(resolved_rainfall_rate)
        return (
            snowfall_rate=isnothing(resolved_snowfall_rate) ? 0.0 : resolved_snowfall_rate,
            rainfall_rate=isnothing(resolved_rainfall_rate) ? 0.0 : resolved_rainfall_rate,
            shortwave_down=resolved_shortwave_down,
        )
    end

    snowfall_fraction = isnothing(resolved_snow_fraction) ?
        _default_snow_fraction(column, air_temperature) :
        resolved_snow_fraction
    rainfall = precipitation_rate * (1.0 - snowfall_fraction)
    snowfall = precipitation_rate - rainfall
    return (
        snowfall_rate=snowfall,
        rainfall_rate=rainfall,
        shortwave_down=resolved_shortwave_down,
    )
end

@inline function _surface_has_snow(column::SnowpackColumn)
    return column.N > 0 && column.mass[1] > EPS_EMPTY_LAYER
end

@inline function _diagnosed_shortwave_down(shortwave_down::Union{Nothing, Float64})
    return isnothing(shortwave_down) ? 400.0 : max(shortwave_down, 0.0)
end

function _run_liquid_water_processes!(
    column::SnowpackColumn,
    liquid_water_before_energy::AbstractVector{Float64},
    dt_seconds::Float64;
    timings::Union{Nothing, StepTimingStats}=nothing,
)
    has_liquid_water = _column_has_liquid_water(column)
    if has_liquid_water
        _time_block!(timings, :percolation) do
            go_percolation!(column)
        end
        has_liquid_water = _column_has_liquid_water(column)
    end

    if column.c.low_density_densification == :htessel &&
       !isempty(liquid_water_before_energy) &&
       has_liquid_water
        _time_block!(timings, :liquid_water_compaction) do
            _apply_htessel_liquid_water_compaction!(
                column,
                liquid_water_before_energy,
                dt_seconds,
            )
        end
        has_liquid_water = _column_has_liquid_water(column)
    end

    if has_liquid_water
        _time_block!(timings, :refreezing) do
            go_refreezing!(column)
        end
    end

    return nothing
end

"""
    step!(column::SnowpackColumn, air_temperature::Float64, precipitation_rate::Float64, dt_days::Float64)

Advance the snowpack column by one time step.

# Arguments
- `column`: The snowpack column to update
- `air_temperature`: Near-surface air temperature [K]
- `precipitation_rate`: Total precipitation rate at surface [kg/m^2/s]
- `dt_days`: Time step [d]
- `snow_fraction`: Optional fraction of precipitation that falls as snow [1]
- `snowfall_rate`: Optional direct snowfall rate [kg/m^2/s]
- `rainfall_rate`: Optional direct rainfall rate [kg/m^2/s]
- `shortwave_down`: Optional downward shortwave forcing [W m^-2], used when
  `q_sw_net` is not provided
- `wind_speed`: Optional near-surface wind speed [m/s], default = `5.0`

Legacy keyword aliases `f_s`, `p_snow`, `p_rain`, and `s_boa` are still accepted.

# Process
1. Apply surface mass flux
2. Handle layer splitting/merging
3. Solve temperature evolution
4. Apply melt, percolation, and refreezing when needed
"""
function step!(
    column::SnowpackColumn,
    air_temperature::Float64,
    precipitation_rate::Float64,
    dt_days::Float64;
    snow_fraction=nothing,
    f_s=nothing,
    P_ave=precipitation_rate,
    snowfall_rate::Union{Nothing, Float64}=nothing,
    rainfall_rate::Union{Nothing, Float64}=nothing,
    shortwave_down::Union{Nothing, Float64}=nothing,
    p_snow::Union{Nothing, Float64}=nothing,
    p_rain::Union{Nothing, Float64}=nothing,
    s_boa::Union{Nothing, Float64}=nothing,
    wind_speed::Float64=10.0,
    q_sw_net::Union{Nothing, Float64}=nothing,
    q_lw_down::Union{Nothing, Float64}=nothing,
    q_sh::Union{Nothing, Float64}=nothing,
    q_lh::Union{Nothing, Float64}=nothing,
    timings::Union{Nothing, StepTimingStats}=nothing,
)
    forcing = _resolve_step_forcing(
        column,
        air_temperature,
        precipitation_rate;
        snow_fraction=snow_fraction,
        f_s=f_s,
        snowfall_rate=snowfall_rate,
        rainfall_rate=rainfall_rate,
        shortwave_down=shortwave_down,
        p_snow=p_snow,
        p_rain=p_rain,
        s_boa=s_boa,
    )

    dt_seconds = dt_days * column.c.seconds_per_day
    started_without_surface_snow = !_surface_has_snow(column)

    _time_block!(timings, :accumulation) do
        apply_accumulation!(
            column,
            forcing.snowfall_rate,
            forcing.rainfall_rate,
            dt_seconds;
            air_temperature=air_temperature,
            wind_speed=wind_speed,
        )
    end

    if forcing.snowfall_rate > 0.0 && started_without_surface_snow && column.N > 0
        column.temperature[1] = air_temperature
    end

    _time_block!(timings, :snow_cover) do
        update_snow_cover!(column)
    end

    if !_surface_has_snow(column)
        column.albedo_dynamic = column.c.alpha_ice
        bare_ice_ablation = _time_block!(timings, :bare_ice_ablation) do
            bare_ice_ablation_mass(
                column,
                air_temperature,
                forcing.rainfall_rate,
                dt_seconds;
                shortwave_down=forcing.shortwave_down,
                q_sw_net=q_sw_net,
                q_lw_down=q_lw_down,
                q_sh=q_sh,
                q_lh=q_lh,
            )
        end
        column.smb_ice -= bare_ice_ablation
        return nothing
    end

    liquid_water_before_energy = column.c.low_density_densification == :htessel ?
        copy(@view column.mass_w[1:column.N]) : Float64[]
    accumulation_rate = max(forcing.snowfall_rate, 0.0) +
                        (_surface_has_snow(column) ? forcing.rainfall_rate : 0.0)

    if _surface_has_snow(column)
        _time_block!(timings, :densification) do
            go_densification!(column, accumulation_rate, dt_seconds)
        end
    end

    diagnosed_shortwave_down = _diagnosed_shortwave_down(forcing.shortwave_down)
    energy = _time_block!(timings, :energy_flux) do
        go_energy_flux!(
            column,
            air_temperature,
            diagnosed_shortwave_down,
            nothing,
            nothing,
            dt_seconds;
            snowfall_rate=forcing.snowfall_rate,
            rainfall_rate=forcing.rainfall_rate,
            diffusion_model=1,
            q_sw_net=q_sw_net,
            q_lw_down=q_lw_down,
            q_sh=q_sh,
            q_lh=q_lh,
            tridiagonal_solver=:thomas,
        )
    end

    if energy.needs_melt
        melt_mass = energy.melt_energy_available / column.c.Lm
        melted_snow = _time_block!(timings, :melt) do
            apply_melt!(column, melt_mass)
        end
        if melted_snow < melt_mass && column.N == 0
            column.smb_ice -= (melt_mass - melted_snow)
        end
    end

    _run_liquid_water_processes!(
        column,
        liquid_water_before_energy,
        dt_seconds;
        timings=timings,
    )

    _time_block!(timings, :snow_cover) do
        update_snow_cover!(column)
    end
    if !_surface_has_snow(column)
        column.albedo_dynamic = column.c.alpha_ice
    end
    return nothing
end
