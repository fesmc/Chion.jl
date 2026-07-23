"""BESSI column stepping and CPU/GPU schedulers."""

@inline _workspace_array(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} =
    similar(storage, NF, dims...)

struct EnergyWorkspace{LT,DT,UT,RT,IT,PT}
    lower::LT
    diag::DT
    upper::UT
    rhs::RT
    interface_conductance::IT
    previous_temperature::PT
end

function EnergyWorkspace(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(storage, NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

EnergyWorkspace(state) =
    EnergyWorkspace(state.mass, number_type(state.c), state.Ntot)

struct ColumnarStepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

function ColumnarStepWorkspace(storage, ::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return ColumnarStepWorkspace(
        _workspace_array(storage, NF, Ntot, ncol),
        EnergyWorkspace(storage, NF, Ntot, ncol),
    )
end

function ColumnarStepWorkspace(state)
    NF = number_type(state.c)
    return ColumnarStepWorkspace(state.mass, NF, state.Ntot, ncols(state))
end

Adapt.@adapt_structure EnergyWorkspace
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
    state,
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
    return column_step_core!(state, idx, subforcing, workspace)
end

"""
    column_step!(state, idx, forcing, config, workspace)

Advance one snowpack column by one forcing step using already-resolved arrays,
constants, and scratch storage. This mutates the supplied state arrays
in-place and may update runoff and SMB diagnostics.
"""
function column_step!(
    state,
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
                    state,
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

    return column_step_core!(state, idx, forcing, workspace)
end

function column_step_core!(
    state,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace,
)
    N_storage = state.N
    mass = state.mass
    mass_w = state.mass_w
    density = state.density
    temperature = state.temperature
    mass_base = state.mass_base
    smb_ice = state.smb_ice
    runoff = state.runoff
    melt = state.melt
    refreezing = state.refreezing
    vapor_mass = state.vapor_mass
    sublimation = state.sublimation
    latent_heat_flux_sum = state.latent_heat_flux_sum
    Tsrf = state.Tsrf
    albedo_dynamic = state.albedo
    c = state.c

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
        state.Ntot,
        state.mass_max,
        state.mass_split,
        state.mass_min,
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
        rainfall_mass = max(forcing.rainfall_rate, zero(forcing.rainfall_rate)) * dt_seconds
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + bare_ice_fluxes.net_mass_change)
        _set_scalar!(melt, idx, _get_scalar(melt, idx) + bare_ice_fluxes.melt_mass)
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + rainfall_mass + bare_ice_fluxes.melt_mass)
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
        state.mass_split,
        state.mass_min,
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
            state.mass_split,
            state.mass_min,
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
    state,
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
            column_step!(state, idx, step_forcing, config, workspace)
        end
    end
end

"""
    _launch_step_columns_kernel!(state, forcing, time_start, time_stop, workspace, active_indices, config)

Launch the backend-specific batch stepping kernel for one contiguous forcing
range and return the KernelAbstractions event.
"""
@inline _step_kernel_workgroupsize(backend) =
    backend isa KernelAbstractions.CPU ? 16 : 256

@inline function _step_time_block_steps(backend)
    name = backend isa KernelAbstractions.CPU ? "CHION_CPU_TIME_BLOCK_STEPS" : "CHION_GPU_TIME_BLOCK_STEPS"
    value = get(ENV, name, "1")
    parsed = tryparse(Int, value)
    return isnothing(parsed) || parsed < 1 ? 1 : parsed
end

@inline function _launch_step_columns_kernel!(
    state,
    forcing::SnowpackForcing,
    time_start::Int,
    time_stop::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
    config,
)
    time_stop >= time_start || return nothing
    backend = _ka_backend(state.mass)
    kernel! = _step_columns_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(
        state,
        workspace,
        active_indices,
        forcing,
        config,
        time_start,
        time_stop,
        ndrange=length(active_indices),
    )
end

@inline function _time_range_bounds(time_range)
    time_range isa AbstractUnitRange ||
        error("BESSI KA stepping requires a contiguous unit-step time range.")
    first_time = Int(first(time_range))
    last_time = Int(last(time_range))
    return first_time, last_time
end

function _step_range!(
    state,
    forcing::SnowpackForcing,
    time_range,
    workspace::ColumnarStepWorkspace,
    active_indices,
    config,
)
    first_time, last_time = _time_range_bounds(time_range)
    first_time <= last_time || return nothing
    backend = _ka_backend(state.mass)
    block_steps = _step_time_block_steps(backend)
    time_index = first_time
    while time_index <= last_time
        time_stop = min(last_time, time_index + block_steps - 1)
        _wait_kernel(_launch_step_columns_kernel!(
            state,
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

"""
    step!(state, forcing, time_index, workspace::ColumnarStepWorkspace)

Advance all columns for one time step using the KernelAbstractions backend
associated with `state.mass`. The same kernel runs on CPU or GPU depending on
the storage backend of the state and workspace arrays.
"""
function step!(
    state,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace;
    kwargs...,
)
    return step!(state, forcing, time_index, workspace, 1:ncols(state); kwargs...)
end

function step!(
    state,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices;
    kwargs...,
)
    config = _step_config_from_keywords(; kwargs...)
    return _step_range!(state, forcing, time_index:time_index, workspace, active_indices, config)
end

"""
    step!(state, forcing, workspace::ColumnarStepWorkspace)

Advance all columns through the full forcing sequence in `forcing`. The outer
time loop runs in Julia, while each time step is advanced by the same
KernelAbstractions batch kernel on the storage backend of `workspace`.
"""
function step!(
    state,
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace;
    kwargs...,
)
    return step!(state, forcing, workspace, 1:ncols(state); kwargs...)
end

function step!(
    state,
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace,
    active_indices;
    kwargs...,
)
    config = _step_config_from_keywords(; kwargs...)
    return _step_range!(state, forcing, 1:_step_time_count(forcing), workspace, active_indices, config)
end
