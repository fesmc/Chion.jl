"""
Batch stepping kernels and workspace-specialized `step!` methods.
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
    f_base_max,
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
            f_base_max,
            forcing,
            workspace,
            update_snow_cover,
        )
    end
end

@kernel function _step_cycle_columns_kernel!(
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
    f_base_max,
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
    dt_days,
    ntime::Int,
    update_snow_cover::Bool,
)
    idx = @index(Global)
    if idx <= length(N_storage)
        for time_index in 1:ntime
            forcing = _step_forcing_from_fields(
                air_temperature[idx, time_index],
                snowfall_rate[idx, time_index],
                rainfall_rate[idx, time_index],
                dt_days[time_index],
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
                f_base_max,
                forcing,
                workspace,
                update_snow_cover,
            )
        end
    end
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    time_index::Int,
    workspaces::AbstractVector{<:StepWorkspace};
    update_snow_cover::Bool=true,
)
    Threads.@threads for idx in 1:column_count(domain)
        step!(
            domain,
            idx,
            _step_forcing_from_fields(forcing, idx, time_index),
            workspaces[Threads.threadid()];
            update_snow_cover=update_snow_cover,
        )
    end
    return nothing
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    time_index::Int,
    workspace::ColumnarStepWorkspace;
    update_snow_cover::Bool=true,
)
    kernel! = _step_columns_kernel!(_ka_backend(domain.mass))
    event = _launch_step_columns_kernel!(
        kernel!,
        domain,
        workspace,
        forcing,
        time_index,
        _step_dt(forcing.dt_days, time_index),
        update_snow_cover,
    )
    _wait_kernel(event)
    return nothing
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    workspaces::AbstractVector{<:StepWorkspace};
    update_snow_cover::Bool=true,
)
    for time_index in 1:_step_time_count(forcing)
        step!(
            domain,
            forcing,
            time_index,
            workspaces;
            update_snow_cover=update_snow_cover,
        )
    end
    return nothing
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackStepFields,
    workspace::ColumnarStepWorkspace;
    update_snow_cover::Bool=true,
)
    kernel! = _step_cycle_columns_kernel!(_ka_backend(domain.mass))
    event = _launch_step_columns_kernel!(
        kernel!,
        domain,
        workspace,
        forcing,
        forcing.dt_days,
        _step_time_count(forcing),
        update_snow_cover,
    )
    _wait_kernel(event)
    return nothing
end
