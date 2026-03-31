"""
Top-level model-step orchestration for domain and raw-array entrypoints.
"""

@inline _default_snow_fraction(c::SnowpackPhysicalConstants, air_temperature) =
    air_temperature > c.T0 ? zero(air_temperature) : one(air_temperature)

function _resolve_step_forcing(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate;
    snow_fraction=nothing,
    f_s=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    shortwave_down=nothing,
    p_snow=nothing,
    p_rain=nothing,
    s_boa=nothing,
)
    resolved_snow_fraction = _resolve_keyword_alias(snow_fraction, f_s, "snow_fraction", "f_s")
    resolved_snowfall_rate = _resolve_keyword_alias(snowfall_rate, p_snow, "snowfall_rate", "p_snow")
    resolved_rainfall_rate = _resolve_keyword_alias(rainfall_rate, p_rain, "rainfall_rate", "p_rain")
    resolved_shortwave_down = _resolve_keyword_alias(shortwave_down, s_boa, "shortwave_down", "s_boa")

    if !isnothing(resolved_snowfall_rate) || !isnothing(resolved_rainfall_rate)
        return (
            snowfall_rate=isnothing(resolved_snowfall_rate) ? zero(precipitation_rate) : resolved_snowfall_rate,
            rainfall_rate=isnothing(resolved_rainfall_rate) ? zero(precipitation_rate) : resolved_rainfall_rate,
            shortwave_down=resolved_shortwave_down,
        )
    end

    snowfall_fraction = isnothing(resolved_snow_fraction) ?
        _default_snow_fraction(c, air_temperature) :
        resolved_snow_fraction
    rainfall = precipitation_rate * (one(precipitation_rate) - snowfall_fraction)
    snowfall = precipitation_rate - rainfall
    return (
        snowfall_rate=snowfall,
        rainfall_rate=rainfall,
        shortwave_down=resolved_shortwave_down,
    )
end

@inline function _diagnosed_shortwave_down(shortwave_down)
    return isnothing(shortwave_down) ? 400.0 : max(shortwave_down, zero(shortwave_down))
end

function _validate_diurnal_configuration(diurnal_shortwave::Bool, latitude, day_of_year)
    if diurnal_shortwave && (isnothing(latitude) || isnothing(day_of_year))
        error("`diurnal_shortwave=true` requires both `latitude` and `day_of_year`.")
    end
    return nothing
end

function _resolved_step_forcing(
    c::SnowpackPhysicalConstants,
    air_temperature,
    precipitation_rate,
    dt_days;
    snow_fraction=nothing,
    f_s=nothing,
    snowfall_rate=nothing,
    rainfall_rate=nothing,
    shortwave_down=nothing,
    p_snow=nothing,
    p_rain=nothing,
    s_boa=nothing,
    wind_speed=oftype(air_temperature, 10.0),
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
    diurnal_shortwave::Bool=false,
    latitude=nothing,
    day_of_year=nothing,
)
    _validate_diurnal_configuration(diurnal_shortwave, latitude, day_of_year)
    resolved = _resolve_step_forcing(
        c,
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
    return SnowpackStepForcing(
        c,
        air_temperature,
        precipitation_rate,
        dt_days;
        snowfall_rate=resolved.snowfall_rate,
        rainfall_rate=resolved.rainfall_rate,
        shortwave_down=_diagnosed_shortwave_down(resolved.shortwave_down),
        wind_speed=wind_speed,
        q_sw_net=q_sw_net,
        q_lw_down=q_lw_down,
        q_sh=q_sh,
        q_lh=q_lh,
        diurnal_shortwave=diurnal_shortwave,
        latitude=isnothing(latitude) ? zero(air_temperature) : latitude,
        day_of_year=isnothing(day_of_year) ? zero(air_temperature) : day_of_year,
    )
end

function _copy_liquid_water_before_energy!(
    liquid_water_before_energy::AbstractVector,
    N_storage,
    mass_w,
    idx::Int,
)
    n = _n_active(N_storage, idx)
    fill!(liquid_water_before_energy, zero(eltype(liquid_water_before_energy)))
    @inbounds for layer_index in 1:n
        liquid_water_before_energy[layer_index] = _get_layer(mass_w, layer_index, idx)
    end
    return n
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
    liquid_water_before_energy::AbstractVector,
    n_liquid_water_before_energy::Int,
    dt_seconds;
    timings=nothing,
)
    has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    if has_liquid_water
        routed_runoff = _time_block!(timings, :percolation) do
            _go_percolation!(N_storage, mass, mass_w, density, idx, c.rho_i, c.rho_w)
        end
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + routed_runoff)
        has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    end

    if c.low_density_densification == :htessel &&
       n_liquid_water_before_energy > 0 &&
       has_liquid_water
        _time_block!(timings, :liquid_water_compaction) do
            _apply_htessel_liquid_water_compaction!(
                N_storage,
                mass,
                mass_w,
                density,
                idx,
                liquid_water_before_energy,
                c.rho_i,
            )
        end
        has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    end

    if has_liquid_water
        _time_block!(timings, :refreezing) do
            _go_refreezing!(N_storage, mass_w, mass, density, temperature, idx, c.T0, c.ci, c.Lm, c.rho_i)
        end
    end

    return nothing
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
    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        max(shortwave_down, zero(dt_seconds)) * (one(dt_seconds) - c.alpha_ice)
    longwave_flux = use_q_lw_down ?
        q_lw_down_value - c.σ * c.ϵ_snow * c.T0^4 :
        c.σ * (c.ϵ_air * air_temperature^4 - c.ϵ_snow * c.T0^4)
    sensible_heat_flux = use_q_sh ? q_sh_value : c.D_sh * (air_temperature - c.T0)
    latent_heat_flux = use_q_lh ? q_lh_value : zero(dt_seconds)
    rain_heat_flux = rainfall_rate * c.cw * (air_temperature - c.T0)
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
    absorbed_shortwave = use_q_sw_net ?
        q_sw_net_value :
        max(shortwave_down, zero(dt_seconds)) * (one(dt_seconds) - c.alpha_ice)
    longwave_flux = use_q_lw_down ?
        q_lw_down_value - c.σ * c.ϵ_snow * c.T0^4 :
        c.σ * (c.ϵ_air * air_temperature^4 - c.ϵ_snow * c.T0^4)
    sensible_heat_flux = use_q_sh ? q_sh_value : c.D_sh * (air_temperature - c.T0)
    latent_heat_flux = use_q_lh ? q_lh_value : zero(dt_seconds)
    rain_heat_flux = rainfall_rate * c.cw * (air_temperature - c.T0)
    partition = _debm_melt_window_fluxes(
        absorbed_shortwave,
        longwave_flux + sensible_heat_flux + latent_heat_flux + rain_heat_flux,
        latitude,
        day_of_year,
    )
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
    longwave_component = use_q_lw_down ?
        q_lw_down_value - c.σ * c.ϵ_snow * surface_temperature^4 :
        c.σ * (c.ϵ_air * air_temperature^4 - c.ϵ_snow * surface_temperature^4)
    sensible_component = use_q_sh ? q_sh_value : c.D_sh * (air_temperature - surface_temperature)
    latent_component = use_q_lh ? q_lh_value : zero(dt_seconds)
    rain_component = rainfall_rate * c.cw * (air_temperature - c.T0)
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

function _step_state_resolved!(
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
    forcing::SnowpackStepForcing,
    workspace::StepWorkspace;
    timings=nothing,
)
    dt_seconds = forcing.dt_days * c.seconds_per_day
    started_without_surface_snow = !_surface_has_snow(N_storage, mass, idx)

    _time_block!(timings, :accumulation) do
        _apply_accumulation!(
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
            forcing.snowfall_rate,
            forcing.rainfall_rate,
            dt_seconds;
            air_temperature=forcing.air_temperature,
            wind_speed=forcing.wind_speed,
        )
    end

    if forcing.snowfall_rate > zero(dt_seconds) &&
       started_without_surface_snow &&
       _n_active(N_storage, idx) > 0
        _set_layer!(temperature, 1, idx, forcing.air_temperature)
    end

    _time_block!(timings, :snow_cover) do
        _update_snow_cover_arrays!(N_storage, mass, mass_w, density, snow_cover, idx)
    end
    _time_block!(timings, :surface_albedo) do
        _update_surface_albedo_arrays!(N_storage, mass, mass_w, density, temperature, albedo_dynamic, idx, c)
    end

    has_surface_snow = _surface_has_snow(N_storage, mass, idx)
    if !has_surface_snow
        bare_ice_ablation = _time_block!(timings, :bare_ice_ablation) do
            if forcing.diurnal_shortwave
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
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - bare_ice_ablation)
        return nothing
    end

    n_liquid_water_before_energy = if c.low_density_densification == :htessel
        _copy_liquid_water_before_energy!(workspace.liquid_water_before_energy, N_storage, mass_w, idx)
    else
        0
    end

    accumulation_rate = max(forcing.snowfall_rate, zero(dt_seconds)) +
                        (has_surface_snow ? forcing.rainfall_rate : zero(dt_seconds))

    if has_surface_snow
        _time_block!(timings, :densification) do
            _go_densification!(N_storage, mass, density, temperature, idx, c, accumulation_rate, dt_seconds)
        end
    end

    latent_heat_linear, latent_heat_constant = _diagnose_latent_heat_flux_coefficients(
        has_surface_snow,
        c,
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
    )
    energy = _time_block!(timings, :energy_flux) do
        _go_energy_flux_resolved!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            Tsrf,
            albedo_dynamic,
            idx,
            c,
            workspace.energy,
            forcing.air_temperature,
            forcing.shortwave_down,
            latent_heat_linear,
            latent_heat_constant,
            dt_seconds,
            1,
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

    extra_melt_energy = zero(dt_seconds)
    if forcing.diurnal_shortwave
        q_sw_effective = forcing.has_q_sw_net ?
            forcing.q_sw_net :
            shortwave_absorbed(forcing.shortwave_down; surface_albedo=_get_scalar(albedo_dynamic, idx))
        adjustment = _diagnose_debm_diurnal_adjustment_resolved(
            c,
            forcing.air_temperature,
            forcing.rainfall_rate,
            dt_seconds,
            _get_layer(temperature, 1, idx),
            q_sw_effective,
            forcing.has_q_lw_down,
            forcing.q_lw_down,
            forcing.has_q_sh,
            forcing.q_sh,
            forcing.has_q_lh,
            forcing.q_lh,
            forcing.latitude,
            forcing.day_of_year,
        )
        extra_melt_energy = adjustment.extra_melt_energy
        if adjustment.refreezing_recharge_energy > zero(dt_seconds)
            _apply_diurnal_refreezing_recharge!(
                N_storage,
                mass,
                mass_w,
                temperature,
                idx,
                c,
                adjustment.refreezing_recharge_energy,
                adjustment.refreezing_period_seconds,
            )
        end
    end

    if energy.needs_melt || extra_melt_energy > zero(dt_seconds)
        melt_mass = (energy.melt_energy_available + extra_melt_energy) / c.Lm
        melted_snow = _time_block!(timings, :melt) do
            _apply_melt!(
                N_storage,
                mass,
                mass_w,
                density,
                temperature,
                runoff,
                Tsrf,
                albedo_dynamic,
                idx,
                mass_split,
                mass_min,
                melt_mass,
                c,
            )
        end
        if melted_snow < melt_mass && _n_active(N_storage, idx) == 0
            _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - (melt_mass - melted_snow))
        end
    end

    _run_liquid_water_processes!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        runoff,
        idx,
        c,
        workspace.liquid_water_before_energy,
        n_liquid_water_before_energy,
        dt_seconds;
        timings=timings,
    )

    _time_block!(timings, :snow_cover) do
        _update_snow_cover_arrays!(N_storage, mass, mass_w, density, snow_cover, idx)
    end
    if !_surface_has_snow(N_storage, mass, idx)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return nothing
end

function step!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    forcing::SnowpackStepForcing,
)
    return step!(domain, idx, forcing, StepWorkspace(domain))
end

function step!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace::StepWorkspace,
)
    return _step_state_resolved!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.mass_base,
        domain.smb_ice,
        domain.runoff,
        domain.Tsrf,
        domain.snow_cover,
        domain.albedo_dynamic,
        idx,
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        domain.f_base_max,
        forcing,
        workspace,
    )
end

function step!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace::StepWorkspace,
    timings::StepTimingStats,
)
    return _step_state_resolved!(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.mass_base,
        domain.smb_ice,
        domain.runoff,
        domain.Tsrf,
        domain.snow_cover,
        domain.albedo_dynamic,
        idx,
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        domain.f_base_max,
        forcing,
        workspace;
        timings=timings,
    )
end

function step!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    air_temperature,
    precipitation_rate,
    dt_days;
    workspace::StepWorkspace=StepWorkspace(domain),
    timings=nothing,
    kwargs...,
)
    forcing = _resolved_step_forcing(
        domain.c,
        air_temperature,
        precipitation_rate,
        dt_days;
        kwargs...,
    )
    return isnothing(timings) ?
        step!(domain, idx, forcing, workspace) :
        step!(domain, idx, forcing, workspace, timings)
end

function step!(
    N::AbstractVector{<:Integer},
    mass::AbstractMatrix,
    mass_w::AbstractMatrix,
    density::AbstractMatrix,
    temperature::AbstractMatrix,
    mass_base::AbstractVector,
    smb_ice::AbstractVector,
    runoff::AbstractVector,
    Tsrf::AbstractVector,
    snow_cover::AbstractVector,
    albedo_dynamic::AbstractVector,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace::StepWorkspace;
    c::SnowpackPhysicalConstants=SnowpackPhysicalConstants(eltype(mass)),
    mass_max=DEFAULT_MASS_MAX,
    mass_split=DEFAULT_MASS_SPLIT,
    mass_min=DEFAULT_MASS_MIN,
    f_base_max=DEFAULT_F_BASE_MAX,
    timings=nothing,
)
    return _step_state_resolved!(
        N,
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
        size(mass, 1),
        mass_max,
        mass_split,
        mass_min,
        f_base_max,
        forcing,
        workspace;
        timings=timings,
    )
end

function step!(
    N::AbstractVector{<:Integer},
    mass::AbstractMatrix,
    mass_w::AbstractMatrix,
    density::AbstractMatrix,
    temperature::AbstractMatrix,
    mass_base::AbstractVector,
    smb_ice::AbstractVector,
    runoff::AbstractVector,
    Tsrf::AbstractVector,
    snow_cover::AbstractVector,
    albedo_dynamic::AbstractVector,
    idx::Int,
    air_temperature,
    precipitation_rate,
    dt_days;
    c::SnowpackPhysicalConstants=SnowpackPhysicalConstants(eltype(mass)),
    mass_max=DEFAULT_MASS_MAX,
    mass_split=DEFAULT_MASS_SPLIT,
    mass_min=DEFAULT_MASS_MIN,
    f_base_max=DEFAULT_F_BASE_MAX,
    workspace::StepWorkspace=StepWorkspace(eltype(mass), size(mass, 1)),
    timings=nothing,
    kwargs...,
)
    forcing = _resolved_step_forcing(c, air_temperature, precipitation_rate, dt_days; kwargs...)
    return step!(
        N,
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
        forcing,
        workspace;
        c=c,
        mass_max=mass_max,
        mass_split=mass_split,
        mass_min=mass_min,
        f_base_max=f_base_max,
        timings=timings,
    )
end
