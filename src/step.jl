"""
Scratch storage for stepping and energy-flux solves.
"""

"""
    _workspace_array(storage, ::Type{NF}, dims...)

Allocate scratch storage compatible with `storage`, preserving the active
backend while changing element type and shape.
"""
@inline _workspace_array(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} =
    similar(storage, NF, dims...)

"""
    EnergyWorkspace

Scratch arrays reused by the implicit temperature solver in
[`go_energy_flux!`](@ref).
"""
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

"""
    EnergyWorkspace(storage, ::Type{NF}, dims...)

Allocate energy-solver scratch arrays on the same backend as `storage`.
Returns an `EnergyWorkspace` whose fields are mutated by the energy-flux
solver.
"""
function EnergyWorkspace(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(storage, NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

"""
    EnergyWorkspace(domain)

Allocate energy-flux scratch storage sized for `domain`.
"""
EnergyWorkspace(domain::AbstractSnowpackDomain) =
    EnergyWorkspace(domain.mass, number_type(domain.c), domain.Ntot)

"""
    ColumnarStepWorkspace

Column-major scratch storage that holds temporary state for every column in a
batch run, regardless of whether the backing arrays live on CPU or GPU.
"""
struct ColumnarStepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

"""
    ColumnarStepWorkspace(storage, ::Type{NF}, Ntot, ncol)

Allocate column-major scratch arrays that hold per-layer temporary state for
every column in a batch on the same backend as `storage`.
"""
function ColumnarStepWorkspace(storage, ::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return ColumnarStepWorkspace(
        _workspace_array(storage, NF, Ntot, ncol),
        EnergyWorkspace(storage, NF, Ntot, ncol),
    )
end

"""
    ColumnarStepWorkspace(domain)

Allocate batch stepping scratch compatible with `domain`'s backend and sized
for all columns.
"""
function ColumnarStepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    return ColumnarStepWorkspace(domain.mass, NF, domain.Ntot, column_count(domain))
end

Adapt.@adapt_structure EnergyWorkspace
Adapt.@adapt_structure ColumnarStepWorkspace

"""
Step-forcing container types and field-to-forcing conversion helpers.
"""

"""
    SnowpackStepForcing{NF}

Per-column forcing bundle consumed by [`step!`](@ref). It stores temperatures,
mass fluxes, optional prescribed surface-flux terms, and metadata needed by
the optional diurnal shortwave adjustment.
"""
struct SnowpackStepForcing{NF <: AbstractFloat}
    air_temperature::NF
    precipitation_rate::NF
    dt_days::NF
    snowfall_rate::NF
    rainfall_rate::NF
    shortwave_down::NF
    wind_speed::NF
    q_sw_net::NF
    q_lw_down::NF
    q_sh::NF
    q_lh::NF
    has_q_sw_net::Bool
    has_q_lw_down::Bool
    has_q_sh::Bool
    has_q_lh::Bool
    diurnal_shortwave_substeps::Bool
    latitude_deg::NF
    day_of_year::NF
    solar_longitude_deg::NF
    diurnal_shortwave_threshold::NF
    diurnal_shortwave_max_substeps::Int
    diurnal_shortwave_min_air_temperature::NF
end

function SnowpackStepForcing(
    air_temperature,
    precipitation_rate,
    dt_days,
    snowfall_rate,
    rainfall_rate,
    shortwave_down,
    wind_speed,
    q_sw_net,
    q_lw_down,
    q_sh,
    q_lh,
    has_q_sw_net::Bool,
    has_q_lw_down::Bool,
    has_q_sh::Bool,
    has_q_lh::Bool,
    diurnal_shortwave_substeps::Bool,
    latitude_deg,
    day_of_year,
    solar_longitude_deg,
)
    return SnowpackStepForcing(
        air_temperature,
        precipitation_rate,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        q_sw_net,
        q_lw_down,
        q_sh,
        q_lh,
        has_q_sw_net,
        has_q_lw_down,
        has_q_sh,
        has_q_lh,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        solar_longitude_deg,
        zero(air_temperature),
        3,
        oftype(air_temperature, 265.15),
    )
end

function SnowpackStepForcing(
    air_temperature,
    precipitation_rate,
    dt_days,
    snowfall_rate,
    rainfall_rate,
    shortwave_down,
    wind_speed,
    q_sw_net,
    q_lw_down,
    q_sh,
    q_lh,
    has_q_sw_net::Bool,
    has_q_lw_down::Bool,
    has_q_sh::Bool,
    has_q_lh::Bool,
    diurnal_shortwave_substeps::Bool,
    latitude_deg,
    day_of_year,
)
    return SnowpackStepForcing(
        air_temperature,
        precipitation_rate,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        q_sw_net,
        q_lw_down,
        q_sh,
        q_lh,
        has_q_sw_net,
        has_q_lw_down,
        has_q_sh,
        has_q_lh,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        _solar_longitude_deg_from_calendar_day(day_of_year),
    )
end

"""
    _step_time_count(fields)

Return the number of forcing time steps stored in `fields`.
"""
@inline _step_time_count(fields) = size(fields.air_temperature, 2)

"""
    _step_dt(dt_days, time_index)

Resolve the step duration in days for `time_index`, supporting both scalar and
vector-valued `dt_days` storage.
"""
@inline _step_dt(dt_days::Number, ::Int) = dt_days
@inline _step_dt(dt_days::AbstractVector, time_index::Int) = @inbounds dt_days[time_index]

"""
    _step_forcing_from_fields(air_temperature, snowfall_rate, rainfall_rate, dt_days, shortwave_down, wind_speed, q_lw_down, has_q_lw_down, q_sh, has_q_sh, q_lh, has_q_lh)

Build a single-column `SnowpackStepForcing` from already-indexed forcing
values. The returned forcing disables optional fluxes that are not present and
sets precipitation rate to snowfall plus rainfall.
"""
@inline function _step_forcing_from_fields(
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    dt_days,
    shortwave_down,
    wind_speed,
    q_lw_down,
    has_q_lw_down::Bool,
    q_sh,
    has_q_sh::Bool,
    q_lh,
    has_q_lh::Bool,
    diurnal_shortwave_substeps::Bool,
    latitude_deg,
    day_of_year,
    solar_longitude_deg,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
)
    return SnowpackStepForcing(
        air_temperature,
        snowfall_rate + rainfall_rate,
        dt_days,
        snowfall_rate,
        rainfall_rate,
        shortwave_down,
        wind_speed,
        zero(air_temperature),
        q_lw_down,
        q_sh,
        q_lh,
        false,
        has_q_lw_down,
        has_q_sh,
        has_q_lh,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        solar_longitude_deg,
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
    )
end

@adapt_structure SnowpackStepForcing

"""
Core stepping flow shared by batch stepping kernels.
"""

function _step_diurnal_shortwave_interval_resolved!(
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
    hour_angle_start,
    hour_angle_end,
    update_snow_cover::Bool,
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
    q_sw_net = forcing.has_q_sw_net ?
        _diurnal_shortwave_interval_average(
            forcing.q_sw_net,
            forcing.latitude_deg,
            forcing.solar_longitude_deg,
            hour_angle_start,
            hour_angle_end,
        ) :
        forcing.q_sw_net
    subforcing = _diurnal_substep_forcing(forcing, fraction, shortwave_down, q_sw_net)
    return _step_state_core_resolved!(
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
        subforcing,
        workspace,
        update_snow_cover,
    )
end

"""
    _step_state_resolved!(..., forcing, workspace, update_snow_cover=true)

Advance one snowpack column by one forcing step using already-resolved arrays,
constants, and scratch storage. This mutates the supplied state arrays
in-place and may update runoff and SMB diagnostics.
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
    update_snow_cover::Bool=true,
)
    if forcing.diurnal_shortwave_substeps
        shortwave_for_criterion = forcing.has_q_sw_net ? forcing.q_sw_net : forcing.shortwave_down
        n_substeps = _diurnal_shortwave_substep_count(
            forcing.dt_days,
            shortwave_for_criterion,
            forcing.air_temperature,
            forcing.diurnal_shortwave_min_air_temperature,
            forcing.latitude_deg,
            forcing.solar_longitude_deg,
            forcing.diurnal_shortwave_threshold,
            forcing.diurnal_shortwave_max_substeps,
        )
        if n_substeps > 1
            day_start = -oftype(forcing.dt_days, π)
            day_end = oftype(forcing.dt_days, π)
            substep_width = (day_end - day_start) / n_substeps
            for substep_index in 1:n_substeps
                hour_angle_start = day_start + (substep_index - 1) * substep_width
                hour_angle_end = substep_index == n_substeps ? day_end : hour_angle_start + substep_width
                _step_diurnal_shortwave_interval_resolved!(
                    N_storage, mass, mass_w, density, temperature,
                    mass_base, smb_ice, runoff, Tsrf, snow_cover, albedo_dynamic,
                    idx, c, Ntot, mass_max, mass_split, mass_min, forcing,
                    workspace, hour_angle_start, hour_angle_end, update_snow_cover && substep_index == n_substeps,
                )
            end
            return nothing
        end
    end

    return _step_state_core_resolved!(
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

function _step_state_core_resolved!(
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
    update_snow_cover::Bool=true,
)
    dt_seconds = forcing.dt_days * c.seconds_per_day
    started_without_surface_snow = !_surface_has_snow(N_storage, mass, idx)

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

    if forcing.snowfall_rate > zero(dt_seconds) &&
       started_without_surface_snow &&
       _n_active(N_storage, idx) > 0
        _set_layer!(temperature, 1, idx, forcing.air_temperature)
    end

    has_surface_snow = _surface_has_snow(N_storage, mass, idx)
    if !has_surface_snow
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        bare_ice_ablation = _bare_ice_ablation_mass(c, forcing, dt_seconds)
        if update_snow_cover
            _update_snow_cover_arrays!(N_storage, mass, mass_w, density, snow_cover, idx)
        end
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - bare_ice_ablation)
        return nothing
    end

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
    )

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
        _go_refreezing!(
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
        has_liquid_water = _column_has_liquid_water(N_storage, mass_w, idx)
    end

    if update_snow_cover
        _update_snow_cover_arrays!(N_storage, mass, mass_w, density, snow_cover, idx)
    end
    if !_surface_has_snow(N_storage, mass, idx)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return nothing
end

"""
Batch stepping over forcing fields through a single KernelAbstractions path.
"""

"""
    _step_columns_kernel!(...)

KernelAbstractions kernel that extracts one time slice of the full forcing
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
    latitude_deg,
    day_of_year,
    solar_longitude_deg,
    diurnal_shortwave_substeps::Bool,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
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
            diurnal_shortwave_substeps,
            latitude_deg[idx, time_index],
            day_of_year,
            solar_longitude_deg,
            diurnal_shortwave_threshold,
            diurnal_shortwave_max_substeps,
            diurnal_shortwave_min_air_temperature,
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
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    update_snow_cover::Bool,
    diurnal_shortwave_substeps::Bool,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
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
        forcing.latitude_deg,
        forcing.day_of_year[time_index],
        forcing.solar_longitude_deg[time_index],
        diurnal_shortwave_substeps,
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
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
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace;
    update_snow_cover::Bool=true,
    diurnal_shortwave_substeps::Bool=false,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps::Int=3,
    diurnal_shortwave_min_air_temperature=265.15,
)
    _wait_kernel(_launch_step_columns_kernel!(
        domain,
        forcing,
        time_index,
        workspace,
        update_snow_cover,
        diurnal_shortwave_substeps,
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
    ))
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
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace;
    update_snow_cover::Bool=true,
    diurnal_shortwave_substeps::Bool=false,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps::Int=3,
    diurnal_shortwave_min_air_temperature=265.15,
)
    for time_index in 1:_step_time_count(forcing)
        step!(
            domain,
            forcing,
            time_index,
            workspace;
            update_snow_cover=update_snow_cover,
            diurnal_shortwave_substeps=diurnal_shortwave_substeps,
            diurnal_shortwave_threshold=diurnal_shortwave_threshold,
            diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
            diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
        )
    end
    return nothing
end
