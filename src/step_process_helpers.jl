"""
Internal helpers shared by the core stepping routine.
"""

function _copy_liquid_water_before_energy!(
    liquid_water_before_energy,
    N_storage,
    mass_w,
    idx::Int,
)
    n = _n_active(N_storage, idx)
    @inbounds for layer_index in 1:n
        _set_layer!(liquid_water_before_energy, layer_index, idx, _get_layer(mass_w, layer_index, idx))
    end
    return n
end

@inline function _apply_accumulation_resolved!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    Tsrf,
    snow_cover,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    Ntot::Int,
    mass_max,
    mass_split,
    mass_min,
    f_base_max,
    snowfall_rate,
    rainfall_rate,
    dt_seconds,
    air_temperature,
    wind_speed,
)
    return _apply_accumulation!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        snow_cover,
        albedo_dynamic,
        idx,
        c,
        Ntot,
        mass_max,
        mass_split,
        mass_min,
        f_base_max,
        snowfall_rate,
        rainfall_rate,
        dt_seconds;
        air_temperature=air_temperature,
        wind_speed=wind_speed,
    )
end

@inline function _bare_ice_ablation_mass_from_forcing(
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

function _run_liquid_water_processes!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    runoff,
    idx::Int,
    c::SnowpackPhysicalConstants,
    liquid_water_before_energy,
    n_liquid_water_before_energy::Int,
    dt_seconds;
    timings=nothing,
)
    has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    if has_liquid_water
        routed_runoff = _time_call!(
            timings,
            :percolation,
            _go_percolation!,
            N_storage,
            mass,
            mass_w,
            density,
            idx,
            c.rho_i,
            c.rho_w,
        )
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + routed_runoff)
        has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    end

    if _uses_htessel_densification(c) &&
       n_liquid_water_before_energy > 0 &&
       has_liquid_water
        _time_call!(
            timings,
            :liquid_water_compaction,
            _apply_htessel_liquid_water_compaction!,
            N_storage,
            mass,
            mass_w,
            density,
            idx,
            liquid_water_before_energy,
            c.rho_i,
        )
    end

    if has_liquid_water
        _time_call!(
            timings,
            :refreezing,
            _go_refreezing!,
            N_storage,
            mass_w,
            mass,
            density,
            temperature,
            idx,
            c.T0,
            c.ci,
            c.Lm,
            c.rho_i,
        )
    end

    return nothing
end

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
