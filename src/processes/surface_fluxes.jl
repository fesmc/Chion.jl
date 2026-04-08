"""
Bare-ice and diurnal surface-flux helpers used by `step.jl`.
"""

"""
    _resolved_nonshortwave_surface_flux_components(c, air_temperature, rainfall_rate, dt_seconds, surface_temperature, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Resolve longwave, sensible, latent, and rain heat flux components for a known
surface temperature. Returned fluxes are instantaneous energy fluxes in model
surface-flux units.
"""
@inline function _resolved_nonshortwave_surface_flux_components(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    surface_temperature,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
)
    longwave_flux = use_q_lw_down ?
        q_lw_down_value - c.σ * c.ϵ_snow * surface_temperature^4 :
        c.σ * (c.ϵ_air * air_temperature^4 - c.ϵ_snow * surface_temperature^4)
    sensible_heat_flux = use_q_sh ? q_sh_value : c.D_sh * (air_temperature - surface_temperature)
    latent_heat_flux = use_q_lh ? q_lh_value : zero(dt_seconds)
    rain_heat_flux = rainfall_rate * c.cw * (air_temperature - c.T0)
    return longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux
end

"""
    _resolved_bare_ice_surface_flux_components(c, air_temperature, rainfall_rate, dt_seconds, shortwave_down, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Resolve all bare-ice surface-flux components, including absorbed shortwave
energy.
"""
@inline function _resolved_bare_ice_surface_flux_components(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    shortwave_down,
    use_q_sw_net::Bool,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
)
    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        max(shortwave_down, zero(dt_seconds)) * (one(dt_seconds) - c.alpha_ice)
    return (
        absorbed_shortwave,
        _resolved_nonshortwave_surface_flux_components(
            c,
            air_temperature,
            rainfall_rate,
            dt_seconds,
            c.T0,
            use_q_lw_down,
            q_lw_down_value,
            use_q_sh,
            q_sh_value,
            use_q_lh,
            q_lh_value,
        )...,
    )
end

"""
    _bare_ice_ablation_mass_resolved(c, air_temperature, rainfall_rate, dt_seconds, shortwave_down, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value)

Convert net positive bare-ice surface energy into melt mass over `dt_seconds`.
Returns zero when the surface energy balance is negative.
"""
function _bare_ice_ablation_mass_resolved(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    shortwave_down,
    use_q_sw_net::Bool,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
)
    absorbed_shortwave, longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux =
        _resolved_bare_ice_surface_flux_components(
            c,
            air_temperature,
            rainfall_rate,
            dt_seconds,
            shortwave_down,
            use_q_sw_net,
            q_sw_net_value,
            use_q_lw_down,
            q_lw_down_value,
            use_q_sh,
            q_sh_value,
            use_q_lh,
            q_lh_value,
        )
    net_surface_flux = absorbed_shortwave + longwave_flux + sensible_heat_flux + latent_heat_flux + rain_heat_flux
    return max(net_surface_flux, zero(net_surface_flux)) * dt_seconds / c.Lm
end

"""
    _bare_ice_ablation_mass_diurnal_resolved(c, air_temperature, rainfall_rate, dt_seconds, shortwave_down, use_q_sw_net, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value, latitude, day_of_year)

Estimate bare-ice ablation using the dEBM-style diurnal melt-window
partitioning.
"""
function _bare_ice_ablation_mass_diurnal_resolved(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    shortwave_down,
    use_q_sw_net::Bool,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
    latitude,
    day_of_year,
)
    absorbed_shortwave, longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux =
        _resolved_bare_ice_surface_flux_components(
            c,
            air_temperature,
            rainfall_rate,
            dt_seconds,
            shortwave_down,
            use_q_sw_net,
            q_sw_net_value,
            use_q_lw_down,
            q_lw_down_value,
            use_q_sh,
            q_sh_value,
            use_q_lh,
            q_lh_value,
        )
    nonshortwave_flux = longwave_flux + sensible_heat_flux + latent_heat_flux + rain_heat_flux
    partition = _debm_melt_window_fluxes(absorbed_shortwave, nonshortwave_flux, latitude, day_of_year)
    return max(partition.melt_window_daily_flux, zero(dt_seconds)) * dt_seconds / c.Lm
end

"""
    _bare_ice_ablation_mass(c, forcing, dt_seconds)

Dispatch bare-ice ablation to the standard or diurnal formulation according to
`forcing.diurnal_shortwave`.
"""
function _bare_ice_ablation_mass(
    c::SnowpackPhysicalConstants,
    forcing::SnowpackStepForcing,
    dt_seconds,
)
    return if forcing.diurnal_shortwave
        _bare_ice_ablation_mass_diurnal_resolved(
            c,
            forcing.air_temperature,
            forcing.rainfall_rate,
            dt_seconds,
            forcing.shortwave_down,
            forcing.has_q_sw_net,
            forcing.q_sw_net,
            forcing.has_q_lw_down,
            forcing.q_lw_down,
            forcing.has_q_sh,
            forcing.q_sh,
            forcing.has_q_lh,
            forcing.q_lh,
            forcing.latitude,
            forcing.day_of_year,
        )
    else
        _bare_ice_ablation_mass_resolved(
            c,
            forcing.air_temperature,
            forcing.rainfall_rate,
            dt_seconds,
            forcing.shortwave_down,
            forcing.has_q_sw_net,
            forcing.q_sw_net,
            forcing.has_q_lw_down,
            forcing.q_lw_down,
            forcing.has_q_sh,
            forcing.q_sh,
            forcing.has_q_lh,
            forcing.q_lh,
        )
    end
end

"""
    _diagnose_debm_diurnal_adjustment_resolved(c, air_temperature, rainfall_rate, dt_seconds, surface_temperature, q_sw_net_value, use_q_lw_down, q_lw_down_value, use_q_sh, q_sh_value, use_q_lh, q_lh_value, latitude, day_of_year)

Diagnose the extra melt energy and refreezing recharge implied by the dEBM
diurnal melt-window partition for an existing snow surface.
"""
function _diagnose_debm_diurnal_adjustment_resolved(
    c::SnowpackPhysicalConstants,
    air_temperature,
    rainfall_rate,
    dt_seconds,
    surface_temperature,
    q_sw_net_value,
    use_q_lw_down::Bool,
    q_lw_down_value,
    use_q_sh::Bool,
    q_sh_value,
    use_q_lh::Bool,
    q_lh_value,
    latitude,
    day_of_year,
)
    longwave_component, sensible_component, latent_component, rain_component =
        _resolved_nonshortwave_surface_flux_components(
            c,
            air_temperature,
            rainfall_rate,
            dt_seconds,
            surface_temperature,
            use_q_lw_down,
            q_lw_down_value,
            use_q_sh,
            q_sh_value,
            use_q_lh,
            q_lh_value,
        )
    baseline_nonshortwave_flux = longwave_component + sensible_component + latent_component + rain_component
    partition = _debm_melt_window_fluxes(q_sw_net_value, baseline_nonshortwave_flux, latitude, day_of_year)
    baseline_positive_flux = max(q_sw_net_value + baseline_nonshortwave_flux, zero(dt_seconds))
    corrected_positive_flux = partition.melt_window_daily_flux
    extra_melt_energy = max(corrected_positive_flux - baseline_positive_flux, zero(dt_seconds)) * dt_seconds
    return (
        baseline_positive_flux=baseline_positive_flux,
        corrected_positive_flux=corrected_positive_flux,
        extra_melt_energy=extra_melt_energy,
        refreezing_recharge_energy=partition.refreezing_daily_flux * dt_seconds,
        refreezing_period_seconds=(one(dt_seconds) - partition.melt_period_fraction) * dt_seconds,
        melt_period_seconds=partition.melt_period_fraction * dt_seconds,
        partition=partition,
    )
end
