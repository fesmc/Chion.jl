"""
Core stepping flow and public single-column entrypoints.
"""

"""
    _step_state_resolved!(..., forcing, workspace, update_snow_cover=true; timings=nothing)

Advance one snowpack column by one forcing step using already-resolved arrays,
constants, and scratch storage. This mutates the supplied state arrays
in-place, may update runoff and SMB diagnostics, and optionally records stage
timings.
"""
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
    forcing::SnowpackStepForcing,
    workspace,
    update_snow_cover::Bool=true;
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

    has_surface_snow = _surface_has_snow(N_storage, mass, idx)
    if !has_surface_snow
        _time_call!(timings, :surface_albedo, _set_scalar!, albedo_dynamic, idx, c.alpha_ice)
        bare_ice_ablation = _time_call!(timings, :bare_ice_ablation, _bare_ice_ablation_mass, c, forcing, dt_seconds)
        if update_snow_cover
            _time_call!(timings, :snow_cover, _update_snow_cover_arrays!, N_storage, mass, mass_w, density, snow_cover, idx)
        end
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - bare_ice_ablation)
        return nothing
    end

    _time_call!(
        timings,
        :surface_albedo,
        _update_surface_albedo_arrays!,
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        albedo_dynamic,
        idx,
        c,
    )

    n_liquid_water_before_energy = 0
    if _uses_htessel_densification(c)
        n_liquid_water_before_energy = _n_active(N_storage, idx)
        @inbounds for layer_index in 1:n_liquid_water_before_energy
            _set_layer!(
                workspace.liquid_water_before_energy,
                layer_index,
                idx,
                _get_layer(mass_w, layer_index, idx),
            )
        end
    end

    accumulation_rate = max(forcing.snowfall_rate, zero(dt_seconds)) +
                        (has_surface_snow ? forcing.rainfall_rate : zero(dt_seconds))
    _time_call!(
        timings,
        :densification,
        _go_densification!,
        N_storage,
        mass,
        density,
        temperature,
        idx,
        c,
        accumulation_rate,
        dt_seconds,
    )

    latent_heat_linear, latent_heat_constant = _diagnose_latent_heat_flux_coefficients(
        has_surface_snow,
        c,
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
    )
    energy = _time_call!(
        timings,
        :energy_flux,
        _go_energy_flux_resolved!,
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
        melted_snow = _time_call!(
            timings,
            :melt,
            _apply_melt!,
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
        if melted_snow < melt_mass && _n_active(N_storage, idx) == 0
            _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - (melt_mass - melted_snow))
        end
    end

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
            workspace.liquid_water_before_energy,
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

    if update_snow_cover
        _time_call!(timings, :snow_cover, _update_snow_cover_arrays!, N_storage, mass, mass_w, density, snow_cover, idx)
    end
    if !_surface_has_snow(N_storage, mass, idx)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return nothing
end

"""
    step!(domain, idx, forcing, workspace=StepWorkspace(domain); timings=nothing, update_snow_cover=true)

Advance column `idx` of `domain` by one step using a prebuilt
`SnowpackStepForcing`. Mutates `domain` in-place and reuses `workspace` for
temporary storage.
"""
function step!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace=StepWorkspace(domain);
    timings=nothing,
    update_snow_cover::Bool=true,
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
        forcing,
        workspace,
        update_snow_cover;
        timings=timings,
    )
end

"""
    step!(domain, idx, air_temperature, precipitation_rate, dt_days; workspace=StepWorkspace(domain), timings=nothing, kwargs...)

Advance column `idx` of `domain` by one step from scalar meteorological input.
Keyword arguments are normalized into a `SnowpackStepForcing` before the core
step routine is called.
"""
function step!(
    domain::AbstractSnowpackDomain,
    idx::Int,
    air_temperature,
    precipitation_rate,
    dt_days;
    workspace=StepWorkspace(domain),
    timings=nothing,
    kwargs...,
)
    return step!(
        domain,
        idx,
        _resolved_step_forcing(
            domain.c,
            air_temperature,
            precipitation_rate,
            dt_days;
            kwargs...,
        ),
        workspace;
        timings=timings,
    )
end
