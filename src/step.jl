"""BESSI column stepping and CPU/GPU schedulers."""

@inline _workspace_array(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} =
    similar(storage, NF, dims...)
@inline _workspace_bool_array(storage, dims::Vararg{Int,N}) where {N} =
    similar(storage, Bool, dims...)
@inline _workspace_index_array(storage, dims::Vararg{Int,N}) where {N} =
    similar(storage, Int, dims...)

struct EnergyWorkspace{LT,DT,UT,RT,IT,PT,TT,KT}
    lower::LT
    diag::DT
    upper::UT
    rhs::RT
    interface_conductance::IT
    previous_temperature::PT
    layer_thickness::TT
    thermal_conductivity::KT
end

function EnergyWorkspace(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(storage, NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

EnergyWorkspace(domain::AbstractSnowpackDomain) =
    EnergyWorkspace(domain.mass, number_type(domain.c), domain.Ntot)

struct EnergyPhaseWorkspace{BT,IT,VT}
    post_active::BT
    solve_active::BT
    needs_second_solve::BT
    needs_melt::BT
    n_liquid_water_before_energy::IT
    surface_mass::VT
    previous_surface_temperature::VT
    surface_diag_term::VT
    surface_flux_constant::VT
    surface_flux_linear::VT
    latent_heat_linear_coefficient::VT
    latent_heat_constant_term::VT
    energy_to_melting::VT
    melt_energy_available::VT
    heating::VT
end

function EnergyPhaseWorkspace(storage, ::Type{NF}, ncol::Int) where {NF <: AbstractFloat}
    bools() = _workspace_bool_array(storage, ncol)
    ints() = _workspace_index_array(storage, ncol)
    values() = _workspace_array(storage, NF, ncol)
    return EnergyPhaseWorkspace(
        bools(),
        bools(),
        bools(),
        bools(),
        ints(),
        values(),
        values(),
        values(),
        values(),
        values(),
        values(),
        values(),
        values(),
        values(),
        values(),
    )
end

struct ColumnarStepWorkspace{LWT,ET,EPT}
    liquid_water_before_energy::LWT
    energy::ET
    energy_phase::EPT
end

function ColumnarStepWorkspace(storage, ::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return ColumnarStepWorkspace(
        _workspace_array(storage, NF, Ntot, ncol),
        EnergyWorkspace(storage, NF, Ntot, ncol),
        EnergyPhaseWorkspace(storage, NF, ncol),
    )
end

function ColumnarStepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    return ColumnarStepWorkspace(domain.mass, NF, domain.Ntot, column_count(domain))
end

Adapt.@adapt_structure EnergyWorkspace
Adapt.@adapt_structure EnergyPhaseWorkspace
Adapt.@adapt_structure ColumnarStepWorkspace

const DEFAULT_BESSI_STEP_OPTIONS = (
    diurnal_shortwave_substeps=false,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps=3,
    diurnal_shortwave_min_air_temperature=265.15,
    diurnal_temperature_cycle=false,
    diurnal_temperature_amplitude=0.0,
)

@inline _step_config_from_keywords(; kwargs...) = merge(DEFAULT_BESSI_STEP_OPTIONS, (; kwargs...))

@inline _step_config_from_keywords(kwargs::NamedTuple) = _step_config_from_keywords(; kwargs...)

"""
Core stepping flow shared by batch stepping kernels.
"""

@inline function _prescribed_surface_albedo(forcing::SnowpackStepForcing)
    return clamp(forcing.prescribed_albedo, zero(forcing.prescribed_albedo), one(forcing.prescribed_albedo))
end

@inline function _set_prescribed_surface_albedo!(albedo_dynamic, idx::Int, forcing::SnowpackStepForcing)
    _set_scalar!(albedo_dynamic, idx, _prescribed_surface_albedo(forcing))
    return nothing
end

function _step_diurnal_shortwave_interval!(
    domain,
    idx::Int,
    forcing::SnowpackStepForcing,
    config,
    workspace,
    hour_angle_start,
    hour_angle_end,
)
    fraction = (hour_angle_end - hour_angle_start) / oftype(forcing.dt_days, 2π)
    fraction <= zero(fraction) && return nothing
    shortwave_down = _diurnal_shortwave_interval_average(
        forcing.shortwave_down,
        forcing.latitude_deg,
        forcing.solar_longitude_deg,
        hour_angle_start,
        hour_angle_end,
    )
    air_temperature = config.diurnal_temperature_cycle ?
        _diurnal_temperature_interval_average(
            forcing.air_temperature,
            config.diurnal_temperature_amplitude,
            hour_angle_start,
            hour_angle_end,
        ) :
        forcing.air_temperature
    q_sw_net = forcing.has_q_sw_net ?
        _diurnal_shortwave_interval_average(
            forcing.q_sw_net,
            forcing.latitude_deg,
            forcing.solar_longitude_deg,
            hour_angle_start,
            hour_angle_end,
        ) :
        forcing.q_sw_net
    subforcing = _diurnal_substep_forcing(forcing, fraction, air_temperature, shortwave_down, q_sw_net)
    return column_step_core!(domain, idx, subforcing, workspace)
end

"""
    column_step!(domain, idx, forcing, config, workspace)

Advance one snowpack column by one forcing step using already-resolved arrays,
constants, and scratch storage. This mutates the supplied state arrays
in-place and may update runoff and SMB diagnostics.
"""
function column_step!(
    domain,
    idx::Int,
    forcing::SnowpackStepForcing,
    config,
    workspace,
)
    if config.diurnal_shortwave_substeps
        shortwave_for_criterion = forcing.has_q_sw_net ? forcing.q_sw_net : forcing.shortwave_down
        n_substeps = _diurnal_shortwave_substep_count(
            forcing.dt_days,
            shortwave_for_criterion,
            forcing.air_temperature,
            config.diurnal_shortwave_min_air_temperature,
            forcing.latitude_deg,
            forcing.solar_longitude_deg,
            config.diurnal_shortwave_threshold,
            config.diurnal_shortwave_max_substeps,
        )
        if n_substeps > 1
            day_start = -oftype(forcing.dt_days, π)
            day_end = oftype(forcing.dt_days, π)
            substep_width = (day_end - day_start) / n_substeps
            for substep_index in 1:n_substeps
                hour_angle_start = day_start + (substep_index - 1) * substep_width
                hour_angle_end = substep_index == n_substeps ? day_end : hour_angle_start + substep_width
                _step_diurnal_shortwave_interval!(
                    domain,
                    idx,
                    forcing,
                    config,
                    workspace,
                    hour_angle_start,
                    hour_angle_end,
                )
            end
            return nothing
        end
    end

    return column_step_core!(domain, idx, forcing, workspace)
end

function _step_column_from_fields!(
    domain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace,
    idx::Int,
    _update_diagnostics::Bool,
    diurnal_shortwave_substeps::Bool,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
    diurnal_temperature_cycle::Bool,
    diurnal_temperature_amplitude,
)
    config = _step_config_from_keywords((
        diurnal_shortwave_substeps=diurnal_shortwave_substeps,
        diurnal_shortwave_threshold=diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle=diurnal_temperature_cycle,
        diurnal_temperature_amplitude=diurnal_temperature_amplitude,
    ))
    return column_step!(domain, idx, _step_forcing_at(forcing, idx, time_index), config, workspace)
end

function column_step_core!(
    domain,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace,
)
    N_storage = domain.N
    mass = domain.mass
    mass_w = domain.mass_w
    density = domain.density
    temperature = domain.temperature
    mass_base = domain.mass_base
    smb_ice = domain.smb_ice
    runoff = domain.runoff
    melt = domain.melt
    refreezing = domain.refreezing
    vapor_mass = domain.vapor_mass
    sublimation = domain.sublimation
    latent_heat_flux_sum = domain.latent_heat_flux_sum
    Tsrf = domain.Tsrf
    albedo_dynamic = domain.albedo_dynamic
    c = domain.c

    dt_seconds = forcing.dt_days * c.seconds_per_day
    started_without_surface_snow = !_surface_has_snow(N_storage, mass, idx)
    use_prescribed_albedo = _uses_prescribed_albedo(c) && forcing.has_prescribed_albedo

    _apply_accumulation_resolved!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        albedo_dynamic,
        idx,
        c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        dt_seconds,
        forcing.air_temperature,
        forcing.wind_speed,
    )

    if forcing.snowfall_rate > zero(dt_seconds) &&
       started_without_surface_snow &&
       _n_active(N_storage, idx) > 0
        _set_layer!(temperature, 1, idx, forcing.air_temperature)
    end

    use_prescribed_albedo && _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)

    has_surface_snow = _surface_has_snow(N_storage, mass, idx)
    if !has_surface_snow
        use_prescribed_albedo || _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        bare_ice_fluxes = _bare_ice_ablation_mass(c, forcing, dt_seconds)
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + bare_ice_fluxes.net_mass_change)
        _set_scalar!(melt, idx, _get_scalar(melt, idx) + bare_ice_fluxes.melt_mass)
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + bare_ice_fluxes.melt_mass)
        _set_scalar!(vapor_mass, idx, _get_scalar(vapor_mass, idx) + bare_ice_fluxes.vapor_mass)
        _set_scalar!(sublimation, idx, _get_scalar(sublimation, idx) + bare_ice_fluxes.sublimation_mass)
        _set_scalar!(latent_heat_flux_sum, idx, _get_scalar(latent_heat_flux_sum, idx) + bare_ice_fluxes.latent_heat_flux * forcing.dt_days)
        return nothing
    end

    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    else
        _update_surface_albedo_arrays!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            albedo_dynamic,
            idx,
            c,
        )
    end

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
    _go_densification!(
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
    energy = _go_energy_flux_resolved!(
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
        forcing.has_q_sw_net,
        forcing.q_sw_net,
        forcing.has_q_lw_down,
        forcing.q_lw_down,
        forcing.has_q_sh,
        forcing.q_sh,
        forcing.has_q_lh,
        forcing.q_lh,
        forcing.has_relative_humidity,
        forcing.relative_humidity,
        forcing.air_pressure,
    )

    snow_vapor_fluxes = _apply_snow_surface_vapor_mass_flux!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        runoff,
        Tsrf,
        albedo_dynamic,
        idx,
        c,
        forcing,
        dt_seconds,
        domain.mass_split,
        domain.mass_min,
    )
    _set_scalar!(vapor_mass, idx, _get_scalar(vapor_mass, idx) + snow_vapor_fluxes.vapor_mass)
    _set_scalar!(sublimation, idx, _get_scalar(sublimation, idx) + snow_vapor_fluxes.sublimation_mass)
    _set_scalar!(latent_heat_flux_sum, idx, _get_scalar(latent_heat_flux_sum, idx) + snow_vapor_fluxes.latent_heat_flux * forcing.dt_days)

    if energy.needs_melt
        melt_mass = energy.melt_energy_available / c.Lm
        melted_snow = _apply_melt!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            runoff,
            Tsrf,
            albedo_dynamic,
            idx,
            domain.mass_split,
            domain.mass_min,
            melt_mass,
            c,
        )
        if melted_snow < melt_mass && _n_active(N_storage, idx) == 0
            ice_melt = melt_mass - melted_snow
            _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - ice_melt)
            _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + ice_melt)
        end
        _set_scalar!(melt, idx, _get_scalar(melt, idx) + melt_mass)
    end

    has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    if has_liquid_water
        routed_runoff = _go_percolation!(
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
        _apply_htessel_liquid_water_compaction!(
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
        refrozen_mass = _go_refreezing!(
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
        _set_scalar!(refreezing, idx, _get_scalar(refreezing, idx) + refrozen_mass)
        has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    end

    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    elseif !_surface_has_snow(N_storage, mass, idx)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return nothing
end

function _pre_energy_column_step!(
    domain,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace,
)
    N_storage = domain.N
    mass = domain.mass
    mass_w = domain.mass_w
    density = domain.density
    temperature = domain.temperature
    mass_base = domain.mass_base
    smb_ice = domain.smb_ice
    runoff = domain.runoff
    melt = domain.melt
    vapor_mass = domain.vapor_mass
    sublimation = domain.sublimation
    latent_heat_flux_sum = domain.latent_heat_flux_sum
    Tsrf = domain.Tsrf
    albedo_dynamic = domain.albedo_dynamic
    c = domain.c
    phase = workspace.energy_phase

    dt_seconds = forcing.dt_days * c.seconds_per_day
    _reset_energy_phase_column!(phase, idx, dt_seconds, zero(dt_seconds), zero(dt_seconds))

    started_without_surface_snow = !_surface_has_snow(N_storage, mass, idx)
    use_prescribed_albedo = _uses_prescribed_albedo(c) && forcing.has_prescribed_albedo

    _apply_accumulation_resolved!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        albedo_dynamic,
        idx,
        c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        dt_seconds,
        forcing.air_temperature,
        forcing.wind_speed,
    )

    if forcing.snowfall_rate > zero(dt_seconds) &&
       started_without_surface_snow &&
       _n_active(N_storage, idx) > 0
        _set_layer!(temperature, 1, idx, forcing.air_temperature)
    end

    use_prescribed_albedo && _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)

    has_surface_snow = _surface_has_snow(N_storage, mass, idx)
    if !has_surface_snow
        use_prescribed_albedo || _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        bare_ice_fluxes = _bare_ice_ablation_mass(c, forcing, dt_seconds)
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + bare_ice_fluxes.net_mass_change)
        _set_scalar!(melt, idx, _get_scalar(melt, idx) + bare_ice_fluxes.melt_mass)
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + bare_ice_fluxes.melt_mass)
        _set_scalar!(vapor_mass, idx, _get_scalar(vapor_mass, idx) + bare_ice_fluxes.vapor_mass)
        _set_scalar!(sublimation, idx, _get_scalar(sublimation, idx) + bare_ice_fluxes.sublimation_mass)
        _set_scalar!(latent_heat_flux_sum, idx, _get_scalar(latent_heat_flux_sum, idx) + bare_ice_fluxes.latent_heat_flux * forcing.dt_days)
        return nothing
    end

    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    else
        _update_surface_albedo_arrays!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            albedo_dynamic,
            idx,
            c,
        )
    end

    if _uses_htessel_densification(c)
        n_liquid_water_before_energy = _n_active(N_storage, idx)
        _set_scalar!(phase.n_liquid_water_before_energy, idx, n_liquid_water_before_energy)
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
    _go_densification!(
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

    return _prepare_energy_flux_system!(
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
        phase,
        forcing.air_temperature,
        forcing.shortwave_down,
        latent_heat_linear,
        latent_heat_constant,
        dt_seconds,
        forcing.has_q_sw_net,
        forcing.q_sw_net,
        forcing.has_q_lw_down,
        forcing.q_lw_down,
        forcing.has_q_sh,
        forcing.q_sh,
        forcing.has_q_lh,
        forcing.q_lh,
        forcing.has_relative_humidity,
        forcing.relative_humidity,
        forcing.air_pressure,
    )
end

function _post_energy_column_step!(
    domain,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace,
)
    phase = workspace.energy_phase
    _get_scalar(phase.post_active, idx) || return nothing

    N_storage = domain.N
    mass = domain.mass
    mass_w = domain.mass_w
    density = domain.density
    temperature = domain.temperature
    smb_ice = domain.smb_ice
    runoff = domain.runoff
    melt = domain.melt
    refreezing = domain.refreezing
    vapor_mass = domain.vapor_mass
    sublimation = domain.sublimation
    latent_heat_flux_sum = domain.latent_heat_flux_sum
    Tsrf = domain.Tsrf
    albedo_dynamic = domain.albedo_dynamic
    c = domain.c

    dt_seconds = forcing.dt_days * c.seconds_per_day
    use_prescribed_albedo = _uses_prescribed_albedo(c) && forcing.has_prescribed_albedo

    snow_vapor_fluxes = _apply_snow_surface_vapor_mass_flux!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        runoff,
        Tsrf,
        albedo_dynamic,
        idx,
        c,
        forcing,
        dt_seconds,
        domain.mass_split,
        domain.mass_min,
    )
    _set_scalar!(vapor_mass, idx, _get_scalar(vapor_mass, idx) + snow_vapor_fluxes.vapor_mass)
    _set_scalar!(sublimation, idx, _get_scalar(sublimation, idx) + snow_vapor_fluxes.sublimation_mass)
    _set_scalar!(latent_heat_flux_sum, idx, _get_scalar(latent_heat_flux_sum, idx) + snow_vapor_fluxes.latent_heat_flux * forcing.dt_days)

    if _get_scalar(phase.needs_melt, idx)
        melt_mass = _get_scalar(phase.melt_energy_available, idx) / c.Lm
        melted_snow = _apply_melt!(
            N_storage,
            mass,
            mass_w,
            density,
            temperature,
            runoff,
            Tsrf,
            albedo_dynamic,
            idx,
            domain.mass_split,
            domain.mass_min,
            melt_mass,
            c,
        )
        if melted_snow < melt_mass && _n_active(N_storage, idx) == 0
            ice_melt = melt_mass - melted_snow
            _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - ice_melt)
            _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + ice_melt)
        end
        _set_scalar!(melt, idx, _get_scalar(melt, idx) + melt_mass)
    end

    has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    if has_liquid_water
        routed_runoff = _go_percolation!(
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
       _get_scalar(phase.n_liquid_water_before_energy, idx) > 0 &&
       has_liquid_water
        _apply_htessel_liquid_water_compaction!(
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
        refrozen_mass = _go_refreezing!(
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
        _set_scalar!(refreezing, idx, _get_scalar(refreezing, idx) + refrozen_mass)
    end

    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    elseif !_surface_has_snow(N_storage, mass, idx)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return nothing
end

"""Batch stepping over forcing fields with one KernelAbstractions path for CPU and GPU."""

"""
    _step_columns_kernel!(...)

KernelAbstractions kernel that extracts one time slice of the full forcing
and advances each column independently in-place.
"""
@kernel function _step_columns_kernel!(
    domain,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing::SnowpackForcing,
    config,
    time_start::Int,
    time_stop::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        for time_index in time_start:time_stop
            step_forcing = _step_forcing_at(forcing, idx, time_index)
            column_step!(domain, idx, step_forcing, config, workspace)
        end
    end
end

@kernel function _step_pre_energy_kernel!(
    domain,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing::SnowpackForcing,
    time_index::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        _pre_energy_column_step!(domain, idx, _step_forcing_at(forcing, idx, time_index), workspace)
    end
end

@kernel function _step_finish_energy_first_kernel!(
    domain,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing::SnowpackForcing,
    time_index::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        dt_seconds = _step_dt(forcing.dt_days, time_index) * domain.c.seconds_per_day
        _finish_energy_flux_first_solve!(
            domain.N,
            domain.mass,
            domain.temperature,
            domain.Tsrf,
            idx,
            domain.c,
            workspace.energy,
            workspace.energy_phase,
            dt_seconds,
        )
    end
end

@kernel function _step_finish_energy_second_kernel!(
    domain,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing::SnowpackForcing,
    time_index::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        dt_seconds = _step_dt(forcing.dt_days, time_index) * domain.c.seconds_per_day
        _finish_energy_flux_second_solve!(
            domain.N,
            domain.temperature,
            domain.Tsrf,
            idx,
            domain.c,
            workspace.energy,
            workspace.energy_phase,
            dt_seconds,
        )
    end
end

@kernel function _step_post_energy_kernel!(
    domain,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing::SnowpackForcing,
    time_index::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        _post_energy_column_step!(domain, idx, _step_forcing_at(forcing, idx, time_index), workspace)
    end
end

"""
    _launch_step_columns_kernel!(domain, forcing, time_start, time_stop, workspace, active_indices, config)

Launch the backend-specific batch stepping kernel for one contiguous forcing
range and return the KernelAbstractions event.
"""
@inline _step_kernel_workgroupsize(backend) =
    backend isa KernelAbstractions.CPU ? 512 : 256

@inline function _cpu_time_block_steps()
    value = get(ENV, "CHION_CPU_TIME_BLOCK_STEPS", "30")
    parsed = tryparse(Int, value)
    return isnothing(parsed) || parsed < 1 ? 1 : parsed
end

@inline _step_time_block_steps(backend) =
    backend isa KernelAbstractions.CPU ? _cpu_time_block_steps() : 1

@inline function _launch_step_columns_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_start::Int,
    time_stop::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
    config,
)
    time_stop >= time_start || return nothing
    backend = _ka_backend(domain.mass)
    kernel! = _step_columns_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(
        domain,
        workspace,
        active_indices,
        forcing,
        config,
        time_start,
        time_stop,
        ndrange=length(active_indices),
    )
end

@inline function _launch_step_pre_energy_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
)
    backend = _ka_backend(domain.mass)
    kernel! = _step_pre_energy_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(domain, workspace, active_indices, forcing, time_index, ndrange=length(active_indices))
end

@inline function _launch_step_finish_energy_first_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
)
    backend = _ka_backend(domain.mass)
    kernel! = _step_finish_energy_first_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(domain, workspace, active_indices, forcing, time_index, ndrange=length(active_indices))
end

@inline function _launch_step_finish_energy_second_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
)
    backend = _ka_backend(domain.mass)
    kernel! = _step_finish_energy_second_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(domain, workspace, active_indices, forcing, time_index, ndrange=length(active_indices))
end

@inline function _launch_step_post_energy_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
)
    backend = _ka_backend(domain.mass)
    kernel! = _step_post_energy_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(domain, workspace, active_indices, forcing, time_index, ndrange=length(active_indices))
end

@inline function _time_range_bounds(time_range)
    time_range isa AbstractUnitRange ||
        error("BESSI KA stepping requires a contiguous unit-step time range.")
    first_time = Int(first(time_range))
    last_time = Int(last(time_range))
    return first_time, last_time
end

function _step_range!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_range,
    workspace::ColumnarStepWorkspace,
    active_indices,
    config,
)
    first_time, last_time = _time_range_bounds(time_range)
    first_time <= last_time || return nothing
    if config.diurnal_shortwave_substeps
        return _step_range_column_kernel!(domain, forcing, time_range, workspace, active_indices, config)
    end

    for time_index in first_time:last_time
        _step_phased_energy_time!(domain, forcing, time_index, workspace, active_indices)
    end
    return nothing
end

function _step_range_column_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_range,
    workspace::ColumnarStepWorkspace,
    active_indices,
    config,
)
    first_time, last_time = _time_range_bounds(time_range)
    first_time <= last_time || return nothing
    backend = _ka_backend(domain.mass)
    block_steps = _step_time_block_steps(backend)
    time_index = first_time
    while time_index <= last_time
        time_stop = min(last_time, time_index + block_steps - 1)
        _wait_kernel(_launch_step_columns_kernel!(
            domain,
            forcing,
            time_index,
            time_stop,
            workspace,
            active_indices,
            config,
        ))
        time_index = time_stop + 1
    end
    return nothing
end

function _step_phased_energy_time!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
)
    energy = workspace.energy
    phase = workspace.energy_phase

    _wait_kernel(_launch_step_pre_energy_kernel!(domain, forcing, time_index, workspace, active_indices))
    _solve_batched_tridiagonal_prefix!(
        energy.rhs,
        energy.lower,
        energy.diag,
        energy.upper,
        energy.rhs,
        energy.previous_temperature,
        domain.N,
        active_indices;
        max_layers=domain.Ntot,
        active_mask=phase.solve_active,
    )
    _wait_kernel(_launch_step_finish_energy_first_kernel!(domain, forcing, time_index, workspace, active_indices))
    _solve_batched_tridiagonal_prefix!(
        energy.rhs,
        energy.lower,
        energy.diag,
        energy.upper,
        energy.rhs,
        energy.previous_temperature,
        domain.N,
        active_indices;
        max_layers=domain.Ntot,
        active_mask=phase.needs_second_solve,
    )
    _wait_kernel(_launch_step_finish_energy_second_kernel!(domain, forcing, time_index, workspace, active_indices))
    _wait_kernel(_launch_step_post_energy_kernel!(domain, forcing, time_index, workspace, active_indices))
    return nothing
end

"""
    step!(domain, forcing, time_index, workspace::ColumnarStepWorkspace)

Advance all columns for one time step using the KernelAbstractions backend
associated with `domain.mass`. The same kernel runs on CPU or GPU depending on
the storage backend of the domain and workspace arrays.
"""
function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace;
    kwargs...,
)
    return step!(domain, forcing, time_index, workspace, 1:column_count(domain); kwargs...)
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices;
    kwargs...,
)
    config = _step_config_from_keywords(; kwargs...)
    return _step_range!(domain, forcing, time_index:time_index, workspace, active_indices, config)
end

"""
    step!(domain, forcing, workspace::ColumnarStepWorkspace)

Advance all columns through the full forcing sequence in `forcing`. The outer
time loop runs in Julia, while each time step is advanced by the same
KernelAbstractions batch kernel on the storage backend of `workspace`.
"""
function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace;
    kwargs...,
)
    return step!(domain, forcing, workspace, 1:column_count(domain); kwargs...)
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace,
    active_indices;
    kwargs...,
)
    config = _step_config_from_keywords(; kwargs...)
    return _step_range!(domain, forcing, 1:_step_time_count(forcing), workspace, active_indices, config)
end

function step_interval_threads!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_range,
    workspace::ColumnarStepWorkspace,
    active_indices;
    kwargs...,
)
    config = _step_config_from_keywords(; kwargs...)
    return _step_range!(domain, forcing, time_range, workspace, active_indices, config)
end
