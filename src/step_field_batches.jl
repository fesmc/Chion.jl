"""
Batch stepping over forcing fields through a single KernelAbstractions path.
"""

"""
    _step_columns_kernel!(...)

KernelAbstractions kernel that extracts one time slice of `SnowpackStepFields`
and advances each column independently in-place.
"""
@kernel function _step_columns_kernel!(
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
    c::SnowpackPhysicalConstants,
    Ntot::Int,
    mass_max,
    mass_split,
    mass_min,
    workspace::ColumnarStepWorkspace,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    shortwave_down,
    wind_speed,
    q_lw_down,
    has_q_lw_down,
    q_sh,
    has_q_sh,
    q_lh,
    has_q_lh,
    time_index::Int,
    dt_days,
    update_snow_cover::Bool,
)
    idx = @index(Global)
    if idx <= length(N_storage)
        forcing = _step_forcing_from_fields(
            air_temperature[idx, time_index],
            snowfall_rate[idx, time_index],
            rainfall_rate[idx, time_index],
            dt_days,
            shortwave_down[idx, time_index],
            wind_speed[idx, time_index],
            q_lw_down[idx, time_index],
            has_q_lw_down[idx, time_index],
            q_sh[idx, time_index],
            has_q_sh[idx, time_index],
            q_lh[idx, time_index],
            has_q_lh[idx, time_index],
        )
        _step_state_resolved!(
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
            forcing,
            workspace,
            update_snow_cover,
        )
    end
end

"""
    _launch_step_columns_kernel!(domain, forcing, time_index, workspace, update_snow_cover)

Launch the backend-specific batch stepping kernel for one forcing time step and
return the KernelAbstractions event.
"""
@inline function _launch_step_columns_kernel!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    update_snow_cover::Bool,
)
    kernel! = _step_columns_kernel!(_ka_backend(domain.mass))
    return kernel!(
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
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        workspace,
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        forcing.shortwave_down,
        forcing.wind_speed,
        forcing.q_lw_down,
        forcing.has_q_lw_down,
        forcing.q_sh,
        forcing.has_q_sh,
        forcing.q_lh,
        forcing.has_q_lh,
        time_index,
        _step_dt(forcing.dt_days, time_index),
        update_snow_cover;
        ndrange=column_count(domain),
    )
end

"""
    step!(domain, forcing, time_index, workspace::ColumnarStepWorkspace; update_snow_cover=true)

Advance all columns for one time step using the KernelAbstractions backend
associated with `domain.mass`. The same kernel runs on CPU or GPU depending on
the storage backend of the domain and workspace arrays.
"""
function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    time_index::Int,
    workspace::ColumnarStepWorkspace;
    update_snow_cover::Bool=true,
)
    _wait_kernel(_launch_step_columns_kernel!(domain, forcing, time_index, workspace, update_snow_cover))
    return nothing
end

"""
    step!(domain, forcing, workspace::ColumnarStepWorkspace; update_snow_cover=true)

Advance all columns through the full forcing sequence in `forcing`. The outer
time loop runs in Julia, while each time step is advanced by the same
KernelAbstractions batch kernel on the storage backend of `workspace`.
"""
function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    workspace::ColumnarStepWorkspace;
    update_snow_cover::Bool=true,
)
    for time_index in 1:_step_time_count(forcing)
        step!(domain, forcing, time_index, workspace; update_snow_cover=update_snow_cover)
    end
    return nothing
end
