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
    relative_humidity::NF
    has_relative_humidity::Bool
    air_pressure::NF
    prescribed_albedo::NF
    has_prescribed_albedo::Bool
    diurnal_shortwave_substeps::Bool
    latitude_deg::NF
    day_of_year::NF
    solar_longitude_deg::NF
    diurnal_shortwave_threshold::NF
    diurnal_shortwave_max_substeps::Int
    diurnal_shortwave_min_air_temperature::NF
    diurnal_temperature_cycle::Bool
    diurnal_temperature_amplitude::NF
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
    relative_humidity,
    has_relative_humidity::Bool,
    air_pressure,
    prescribed_albedo,
    has_prescribed_albedo::Bool,
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
        relative_humidity,
        has_relative_humidity,
        air_pressure,
        prescribed_albedo,
        has_prescribed_albedo,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        solar_longitude_deg,
        zero(air_temperature),
        3,
        oftype(air_temperature, 265.15),
        false,
        zero(air_temperature),
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
    prescribed_albedo,
    has_prescribed_albedo::Bool,
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
        zero(air_temperature),
        false,
        oftype(air_temperature, 101_325.0),
        prescribed_albedo,
        has_prescribed_albedo,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        solar_longitude_deg,
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
    relative_humidity,
    has_relative_humidity::Bool,
    air_pressure,
    prescribed_albedo,
    has_prescribed_albedo::Bool,
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
        relative_humidity,
        has_relative_humidity,
        air_pressure,
        prescribed_albedo,
        has_prescribed_albedo,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        _solar_longitude_deg_from_calendar_day(day_of_year),
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
    prescribed_albedo,
    has_prescribed_albedo::Bool,
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
        zero(air_temperature),
        false,
        oftype(air_temperature, 101_325.0),
        prescribed_albedo,
        has_prescribed_albedo,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
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
    _step_forcing_from_fields(air_temperature, snowfall_rate, rainfall_rate, dt_days, shortwave_down, wind_speed, q_lw_down, has_q_lw_down, q_sh, has_q_sh, q_lh, has_q_lh, prescribed_albedo, has_prescribed_albedo)

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
    relative_humidity,
    has_relative_humidity::Bool,
    air_pressure,
    prescribed_albedo,
    has_prescribed_albedo::Bool,
    diurnal_shortwave_substeps::Bool,
    latitude_deg,
    day_of_year,
    solar_longitude_deg,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
    diurnal_temperature_cycle::Bool,
    diurnal_temperature_amplitude,
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
        relative_humidity,
        has_relative_humidity,
        air_pressure,
        prescribed_albedo,
        has_prescribed_albedo,
        diurnal_shortwave_substeps,
        latitude_deg,
        day_of_year,
        solar_longitude_deg,
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle,
        diurnal_temperature_amplitude,
    )
end

@adapt_structure SnowpackStepForcing

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

function _step_diurnal_shortwave_interval_resolved!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    melt,
    refreezing,
    vapor_mass,
    sublimation,
    latent_heat_flux_sum,
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
    air_temperature = forcing.diurnal_temperature_cycle ?
        _diurnal_temperature_interval_average(
            forcing.air_temperature,
            forcing.diurnal_temperature_amplitude,
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
    return _step_state_core_resolved!(
        N_storage,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        melt,
        refreezing,
        vapor_mass,
        sublimation,
        latent_heat_flux_sum,
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
    melt,
    refreezing,
    vapor_mass,
    sublimation,
    latent_heat_flux_sum,
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
                    mass_base, smb_ice, runoff, melt, refreezing, vapor_mass, sublimation, latent_heat_flux_sum, Tsrf, snow_cover, albedo_dynamic,
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
        melt,
        refreezing,
        vapor_mass,
        sublimation,
        latent_heat_flux_sum,
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
    melt,
    refreezing,
    vapor_mass,
    sublimation,
    latent_heat_flux_sum,
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
        if update_snow_cover
            _update_snow_cover_arrays!(N_storage, mass, mass_w, density, snow_cover, idx)
        end
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
        mass_split,
        mass_min,
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
            mass_split,
            mass_min,
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

    if update_snow_cover
        _update_snow_cover_arrays!(N_storage, mass, mass_w, density, snow_cover, idx)
    end
    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    elseif !_surface_has_snow(N_storage, mass, idx)
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
    end

    return nothing
end

"""
Batch stepping over forcing fields.

GPU-backed state uses KernelAbstractions kernels. CPU-backed BESSI state uses
plain Julia threaded loops over chunks of columns to avoid per-step KA launch
overhead and keep the branch-heavy scalar physics on the normal CPU compiler
path.
"""

Base.@propagate_inbounds @inline function _step_column_from_fields!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    idx::Int,
    update_snow_cover::Bool,
    diurnal_shortwave_substeps::Bool,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
    diurnal_temperature_cycle::Bool,
    diurnal_temperature_amplitude,
)
    step_forcing = _step_forcing_from_fields(
        forcing.air_temperature[idx, time_index],
        forcing.snowfall_rate[idx, time_index],
        forcing.rainfall_rate[idx, time_index],
        _step_dt(forcing.dt_days, time_index),
        forcing.shortwave_down[idx, time_index],
        forcing.wind_speed[idx, time_index],
        forcing.q_lw_down[idx, time_index],
        forcing.has_q_lw_down[idx, time_index],
        forcing.q_sh[idx, time_index],
        forcing.has_q_sh[idx, time_index],
        forcing.q_lh[idx, time_index],
        forcing.has_q_lh[idx, time_index],
        forcing.relative_humidity[idx, time_index],
        forcing.has_relative_humidity[idx, time_index],
        forcing.air_pressure[idx, time_index],
        forcing.prescribed_albedo[idx, time_index],
        forcing.has_prescribed_albedo[idx, time_index],
        diurnal_shortwave_substeps,
        forcing.latitude_deg[idx, time_index],
        forcing.day_of_year[time_index],
        forcing.solar_longitude_deg[time_index],
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle,
        diurnal_temperature_amplitude,
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
        domain.melt,
        domain.refreezing,
        domain.vapor_mass,
        domain.sublimation,
        domain.latent_heat_flux_sum,
        domain.Tsrf,
        domain.snow_cover,
        domain.albedo_dynamic,
        idx,
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        step_forcing,
        workspace,
        update_snow_cover,
    )
end

@inline function _default_step_threads_chunk_size(nactive::Int)
    nactive <= 0 && return 1
    nchunks = max(1, Base.Threads.nthreads() * 4)
    return max(1, min(nactive, max(16, cld(nactive, nchunks))))
end

"""
    step_interval_threads!(domain, forcing, time_range, workspace[, active_indices]; chunk_size)

Advance BESSI columns over a forcing interval using Base.Threads on CPU arrays.
The loop order is time-within-column-chunk, so forcing reads stay contiguous
for the `ncol x ntime` forcing layout while each column mutates its contiguous
`Ntot x ncol` state slice.
"""
function step_interval_threads!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_range,
    workspace::ColumnarStepWorkspace,
    active_indices=1:column_count(domain);
    update_snow_cover::Bool=true,
    diurnal_shortwave_substeps::Bool=false,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps::Int=3,
    diurnal_shortwave_min_air_temperature=265.15,
    diurnal_temperature_cycle::Bool=false,
    diurnal_temperature_amplitude=0.0,
    chunk_size::Union{Nothing,Integer}=nothing,
)
    domain.mass isa Array || error("step_interval_threads! is CPU-only; use step! for GPU-backed BESSI state.")
    nactive = length(active_indices)
    nactive == 0 && return nothing
    resolved_chunk_size = isnothing(chunk_size) ? _default_step_threads_chunk_size(nactive) : Int(chunk_size)
    resolved_chunk_size > 0 || error("`chunk_size` must be positive.")
    time_indices = time_range isa AbstractUnitRange ? time_range : collect(time_range)
    chunks = Iterators.partition(active_indices, resolved_chunk_size)
    chunk_list = collect(chunks)
    @threads :dynamic for chunk_index in eachindex(chunk_list)
        cols = chunk_list[chunk_index]
        for time_index in time_indices
            @inbounds for idx in cols
                _step_column_from_fields!(
                    domain,
                    forcing,
                    Int(time_index),
                    workspace,
                    Int(idx),
                    update_snow_cover,
                    diurnal_shortwave_substeps,
                    diurnal_shortwave_threshold,
                    diurnal_shortwave_max_substeps,
                    diurnal_shortwave_min_air_temperature,
                    diurnal_temperature_cycle,
                    diurnal_temperature_amplitude,
                )
            end
        end
    end
    return nothing
end

step_year_threads!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace,
    active_indices=1:column_count(domain);
    kwargs...,
) = step_interval_threads!(domain, forcing, 1:_step_time_count(forcing), workspace, active_indices; kwargs...)

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
    melt,
    refreezing,
    vapor_mass,
    sublimation,
    latent_heat_flux_sum,
    Tsrf,
    snow_cover,
    albedo_dynamic,
    c::SnowpackPhysicalConstants,
    Ntot::Int,
    mass_max,
    mass_split,
    mass_min,
    workspace::ColumnarStepWorkspace,
    active_indices,
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
    relative_humidity,
    has_relative_humidity,
    air_pressure,
    prescribed_albedo,
    has_prescribed_albedo,
    latitude_deg,
    day_of_year,
    solar_longitude_deg,
    diurnal_shortwave_substeps::Bool,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
    diurnal_temperature_cycle::Bool,
    diurnal_temperature_amplitude,
    time_index::Int,
    dt_days,
    update_snow_cover::Bool,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        idx = active_indices[active_idx]
        forcing = _step_forcing_from_fields(
            air_temperature[idx, time_index],
            snowfall_rate[idx, time_index],
            rainfall_rate[idx, time_index],
            _step_dt(dt_days, time_index),
            shortwave_down[idx, time_index],
            wind_speed[idx, time_index],
            q_lw_down[idx, time_index],
            has_q_lw_down[idx, time_index],
            q_sh[idx, time_index],
            has_q_sh[idx, time_index],
            q_lh[idx, time_index],
            has_q_lh[idx, time_index],
            relative_humidity[idx, time_index],
            has_relative_humidity[idx, time_index],
            air_pressure[idx, time_index],
            prescribed_albedo[idx, time_index],
            has_prescribed_albedo[idx, time_index],
            diurnal_shortwave_substeps,
            latitude_deg[idx, time_index],
            day_of_year[time_index],
            solar_longitude_deg[time_index],
            diurnal_shortwave_threshold,
            diurnal_shortwave_max_substeps,
            diurnal_shortwave_min_air_temperature,
            diurnal_temperature_cycle,
            diurnal_temperature_amplitude,
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
            melt,
            refreezing,
            vapor_mass,
            sublimation,
            latent_heat_flux_sum,
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
    active_indices,
    update_snow_cover::Bool,
    diurnal_shortwave_substeps::Bool,
    diurnal_shortwave_threshold,
    diurnal_shortwave_max_substeps::Int,
    diurnal_shortwave_min_air_temperature,
    diurnal_temperature_cycle::Bool,
    diurnal_temperature_amplitude,
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
        domain.melt,
        domain.refreezing,
        domain.vapor_mass,
        domain.sublimation,
        domain.latent_heat_flux_sum,
        domain.Tsrf,
        domain.snow_cover,
        domain.albedo_dynamic,
        domain.c,
        domain.Ntot,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        workspace,
        active_indices,
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
        forcing.relative_humidity,
        forcing.has_relative_humidity,
        forcing.air_pressure,
        forcing.prescribed_albedo,
        forcing.has_prescribed_albedo,
        forcing.latitude_deg,
        forcing.day_of_year,
        forcing.solar_longitude_deg,
        diurnal_shortwave_substeps,
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle,
        diurnal_temperature_amplitude,
        time_index,
        forcing.dt_days,
        update_snow_cover;
        ndrange=length(active_indices),
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
    diurnal_temperature_cycle::Bool=false,
    diurnal_temperature_amplitude=0.0,
)
    active_indices = 1:column_count(domain)
    return step!(
        domain,
        forcing,
        time_index,
        workspace,
        active_indices;
        update_snow_cover=update_snow_cover,
        diurnal_shortwave_substeps=diurnal_shortwave_substeps,
        diurnal_shortwave_threshold=diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle=diurnal_temperature_cycle,
        diurnal_temperature_amplitude=diurnal_temperature_amplitude,
    )
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    time_index::Int,
    workspace::ColumnarStepWorkspace,
    active_indices;
    update_snow_cover::Bool=true,
    diurnal_shortwave_substeps::Bool=false,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps::Int=3,
    diurnal_shortwave_min_air_temperature=265.15,
    diurnal_temperature_cycle::Bool=false,
    diurnal_temperature_amplitude=0.0,
)
    if domain.mass isa Array
        return step_interval_threads!(
            domain,
            forcing,
            time_index:time_index,
            workspace,
            active_indices;
            update_snow_cover=update_snow_cover,
            diurnal_shortwave_substeps=diurnal_shortwave_substeps,
            diurnal_shortwave_threshold=diurnal_shortwave_threshold,
            diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
            diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
            diurnal_temperature_cycle=diurnal_temperature_cycle,
            diurnal_temperature_amplitude=diurnal_temperature_amplitude,
        )
    end
    _wait_kernel(_launch_step_columns_kernel!(
        domain,
        forcing,
        time_index,
        workspace,
        active_indices,
        update_snow_cover,
        diurnal_shortwave_substeps,
        diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle,
        diurnal_temperature_amplitude,
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
    diurnal_temperature_cycle::Bool=false,
    diurnal_temperature_amplitude=0.0,
)
    active_indices = 1:column_count(domain)
    return step!(
        domain,
        forcing,
        workspace,
        active_indices;
        update_snow_cover=update_snow_cover,
        diurnal_shortwave_substeps=diurnal_shortwave_substeps,
        diurnal_shortwave_threshold=diurnal_shortwave_threshold,
        diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
        diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
        diurnal_temperature_cycle=diurnal_temperature_cycle,
        diurnal_temperature_amplitude=diurnal_temperature_amplitude,
    )
end

function step!(
    domain::AbstractSnowpackDomain,
    forcing::SnowpackForcing,
    workspace::ColumnarStepWorkspace,
    active_indices;
    update_snow_cover::Bool=true,
    diurnal_shortwave_substeps::Bool=false,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps::Int=3,
    diurnal_shortwave_min_air_temperature=265.15,
    diurnal_temperature_cycle::Bool=false,
    diurnal_temperature_amplitude=0.0,
)
    if domain.mass isa Array
        return step_interval_threads!(
            domain,
            forcing,
            1:_step_time_count(forcing),
            workspace,
            active_indices;
            update_snow_cover=update_snow_cover,
            diurnal_shortwave_substeps=diurnal_shortwave_substeps,
            diurnal_shortwave_threshold=diurnal_shortwave_threshold,
            diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
            diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
            diurnal_temperature_cycle=diurnal_temperature_cycle,
            diurnal_temperature_amplitude=diurnal_temperature_amplitude,
        )
    end
    for time_index in 1:_step_time_count(forcing)
        step!(
            domain,
            forcing,
            time_index,
            workspace,
            active_indices;
            update_snow_cover=update_snow_cover,
            diurnal_shortwave_substeps=diurnal_shortwave_substeps,
            diurnal_shortwave_threshold=diurnal_shortwave_threshold,
            diurnal_shortwave_max_substeps=diurnal_shortwave_max_substeps,
            diurnal_shortwave_min_air_temperature=diurnal_shortwave_min_air_temperature,
            diurnal_temperature_cycle=diurnal_temperature_cycle,
            diurnal_temperature_amplitude=diurnal_temperature_amplitude,
        )
    end
    return nothing
end
