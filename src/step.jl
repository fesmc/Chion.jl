"""BESSI column stepping and CPU/GPU schedulers."""

@inline _workspace_array(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} =
    similar(storage, NF, dims...)

struct EnergyWorkspace{A}
    lower::A
    diag::A
    upper::A
    rhs::A
    interface_conductance::A
    previous_temperature::A
end

function EnergyWorkspace(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(storage, NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

# Thermal rows cover the snow/firn layers followed by the ice substrate.
@inline _thermal_row_capacity(state) = state.Ntot + size(state.ice_temperature, 1)

EnergyWorkspace(state) =
    EnergyWorkspace(state.mass, number_type(state.c), _thermal_row_capacity(state))

struct ColumnarStepWorkspace{LWT,ET,LWE,SWE,SHE,RHE}
    liquid_water_before_energy::LWT
    energy::ET
    # Cumulative net-longwave energy (J m⁻²), recorded from the same
    # linearized boundary condition used by the energy solver.
    net_longwave_energy::LWE
    absorbed_shortwave_energy::SWE
    sensible_heat_energy::SHE
    rain_heat_energy::RHE
end

function ColumnarStepWorkspace(storage, ::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    net_longwave_energy = _workspace_array(storage, NF, ncol)
    absorbed_shortwave_energy = _workspace_array(storage, NF, ncol)
    sensible_heat_energy = _workspace_array(storage, NF, ncol)
    rain_heat_energy = _workspace_array(storage, NF, ncol)
    fill!(net_longwave_energy, zero(NF))
    fill!(absorbed_shortwave_energy, zero(NF))
    fill!(sensible_heat_energy, zero(NF))
    fill!(rain_heat_energy, zero(NF))
    return ColumnarStepWorkspace(
        _workspace_array(storage, NF, Ntot, ncol),
        EnergyWorkspace(storage, NF, Ntot, ncol),
        net_longwave_energy,
        absorbed_shortwave_energy,
        sensible_heat_energy,
        rain_heat_energy,
    )
end

function ColumnarStepWorkspace(state)
    NF = number_type(state.c)
    return ColumnarStepWorkspace(state.mass, NF, _thermal_row_capacity(state), ncols(state))
end

Adapt.@adapt_structure EnergyWorkspace
Adapt.@adapt_structure ColumnarStepWorkspace

const DEFAULT_BESSI_STEP_OPTIONS = (
    diurnal_shortwave_substeps=Val(false),
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps=3,
    diurnal_shortwave_min_air_temperature=265.15,
    diurnal_temperature_cycle=Val(false),
    diurnal_temperature_amplitude=0.0,
    diurnal_temperature_amplitude_gradient=0.0,
    diurnal_temperature_amplitude_reference_height=0.0,
    diurnal_temperature_amplitude_max=Inf,
)

@inline function _step_config_from_keywords(; kwargs...)
    config = merge(DEFAULT_BESSI_STEP_OPTIONS, (; kwargs...))
    return merge(
        config,
        (
            diurnal_shortwave_substeps=Val(_config_enabled(config.diurnal_shortwave_substeps)),
            diurnal_temperature_cycle=Val(_config_enabled(config.diurnal_temperature_cycle)),
        ),
    )
end

@inline _step_config_from_keywords(kwargs::NamedTuple) = _step_config_from_keywords(; kwargs...)

@inline _config_enabled(::Val{Enabled}) where {Enabled} = Enabled
@inline _config_enabled(enabled::Bool) = enabled
@inline _diurnal_air_temperature(::Val{false}, forcing, config, hour_angle_start, hour_angle_end) =
    forcing.air_temperature
@inline function _diurnal_air_temperature(::Val{true}, forcing, config, hour_angle_start, hour_angle_end)
    # Forcing without surface height (NaN) gets no elevation adjustment.
    height_excess = ifelse(
        isfinite(forcing.surface_height),
        max(forcing.surface_height - config.diurnal_temperature_amplitude_reference_height, zero(forcing.surface_height)),
        zero(forcing.surface_height),
    )
    return _diurnal_temperature_interval_average(
        forcing.air_temperature,
        clamp(
            config.diurnal_temperature_amplitude + config.diurnal_temperature_amplitude_gradient * height_excess,
            zero(forcing.air_temperature),
            config.diurnal_temperature_amplitude_max,
        ),
        hour_angle_start,
        hour_angle_end,
    )
end

@inline _copy_liquid_water_for_compaction!(::Val{:bessi}, args...) = 0
@inline function _copy_liquid_water_for_compaction!(
    ::Val{:htessel},
    N_storage,
    mass_w,
    liquid_water_before_energy,
    idx,
)
    n = _n_active(N_storage, idx)
    @inbounds for layer_index in 1:n
        _set_layer!(
            liquid_water_before_energy,
            layer_index,
            idx,
            _get_layer(mass_w, layer_index, idx),
        )
    end
    return n
end

@inline _apply_liquid_water_compaction!(::Val{:bessi}, args...) = nothing
@inline function _apply_liquid_water_compaction!(
    ::Val{:htessel},
    N_storage,
    mass,
    mass_w,
    density,
    idx,
    liquid_water_before_energy,
    ice_density,
)
    return _apply_htessel_liquid_water_compaction!(
        N_storage,
        mass,
        mass_w,
        density,
        idx,
        liquid_water_before_energy,
        ice_density,
    )
end

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

@inline function _step_diurnal_shortwave_interval!(
    fields,
    parameters::BESSIParameters,
    idx::Int,
    forcing::SnowpackStepForcing,
    config,
    workspace,
    hour_angle_start,
    hour_angle_end,
)
    fraction = (hour_angle_end - hour_angle_start) / oftype(forcing.dt_days, 2π)
    shortwave_down = _diurnal_shortwave_interval_average(
        forcing.shortwave_down, forcing.latitude_deg,
        forcing.solar_longitude_deg, hour_angle_start, hour_angle_end,
    )
    air_temperature = _diurnal_air_temperature(
        config.diurnal_temperature_cycle,
        forcing,
        config,
        hour_angle_start,
        hour_angle_end,
    )
    reconstructed_q_sw_net = _diurnal_shortwave_interval_average(
        forcing.q_sw_net, forcing.latitude_deg,
        forcing.solar_longitude_deg, hour_angle_start, hour_angle_end,
    )
    q_sw_net = ifelse(forcing.has_q_sw_net, reconstructed_q_sw_net, forcing.q_sw_net)
    # Preserve the daily prescribed melt integral while allocating it by the
    # reconstructed absorbed shortwave energy.  This is deliberately a
    # hydrological forcing experiment: it does not make the SEB energy budget
    # close to the imposed MAR melt.
    melt_fraction = ifelse(
        forcing.has_prescribed_melt & forcing.has_q_sw_net & (forcing.q_sw_net > zero(forcing.q_sw_net)),
        max(q_sw_net, zero(q_sw_net)) * fraction / forcing.q_sw_net,
        fraction,
    )
    prescribed_melt = ifelse(
        forcing.has_prescribed_melt,
        forcing.prescribed_melt * melt_fraction / fraction,
        forcing.prescribed_melt,
    )
    subforcing = _diurnal_substep_forcing(forcing, fraction, air_temperature, shortwave_down, q_sw_net, prescribed_melt)
    return column_step_core!(fields, parameters, idx, subforcing, workspace)
end

"""
    column_step!(state, idx, forcing, config, workspace)

Advance one snowpack column by one forcing step using already-resolved arrays,
constants, and scratch storage. This mutates the supplied state arrays
in-place and may update runoff and SMB diagnostics.
"""
@inline function column_step!(
    state::BESSIState,
    idx::Int,
    forcing::SnowpackStepForcing,
    config,
    workspace,
)
    return column_step!(get_fields(state), state.parameters, idx, forcing, config, workspace)
end

@inline function column_step!(
    fields,
    parameters::BESSIParameters,
    idx::Int,
    forcing::SnowpackStepForcing,
    config,
    workspace,
)
    return _column_step_diurnal!(
        config.diurnal_shortwave_substeps,
        fields,
        parameters,
        idx,
        _with_parameterized_longwave(parameters.c, forcing),
        config,
        workspace,
    )
end

@inline function _column_step_diurnal!(
    ::Val{false},
    fields,
    parameters,
    idx,
    forcing,
    config,
    workspace,
)
    return column_step_core!(fields, parameters, idx, forcing, workspace)
end

@inline function _column_step_diurnal!(
    ::Val{true},
    fields,
    parameters,
    idx,
    forcing,
    config,
    workspace,
)
    shortwave_for_criterion = ifelse(forcing.has_q_sw_net, forcing.q_sw_net, forcing.shortwave_down)
    n_substeps = _diurnal_shortwave_substep_count(
        forcing.dt_days, shortwave_for_criterion, forcing.air_temperature,
        config.diurnal_shortwave_min_air_temperature, forcing.latitude_deg,
        forcing.solar_longitude_deg, config.diurnal_shortwave_threshold,
        config.diurnal_shortwave_max_substeps,
    )
    day_start = -oftype(forcing.dt_days, π)
    day_end = oftype(forcing.dt_days, π)
    substep_width = (day_end - day_start) / n_substeps
    for substep_index in 1:n_substeps
        hour_angle_start = day_start + (substep_index - 1) * substep_width
        hour_angle_end = ifelse(
            substep_index == n_substeps,
            day_end,
            hour_angle_start + substep_width,
        )
        _step_diurnal_shortwave_interval!(
            fields,
            parameters,
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

"""
    _step_bare_ice_substrate!(fields, parameters, idx, forcing, workspace, dt_seconds)

Bare-ice step with a thermal ice substrate. The top substrate layer is the
surface: it cools when the surface energy balance is negative and must be
re-warmed before melting, instead of being held at the melting point with the
negative energy discarded. Melt and rain leave as runoff from the impermeable
surface.
"""
@inline function _step_bare_ice_substrate!(fields, parameters::BESSIParameters, idx::Int, forcing::SnowpackStepForcing, workspace, dt_seconds)
    c = parameters.c
    Tsrf = fields.Tsrf
    _set_scalar!(fields.albedo, idx, _bare_ice_albedo(c, forcing))
    # Rain enthalpy enters the surface balance as on a snow surface.
    latent_heat_linear, latent_heat_constant = _diagnose_latent_heat_flux_coefficients(
        true, c, forcing.air_temperature, forcing.snowfall_rate, forcing.rainfall_rate,
    )
    energy = _go_energy_flux_resolved!(
        fields.N, fields.mass, fields.mass_w, fields.density, fields.temperature,
        Tsrf, fields.albedo, idx, c, workspace.energy,
        forcing.air_temperature, forcing.shortwave_down,
        latent_heat_linear, latent_heat_constant, dt_seconds,
        forcing.has_q_sw_net, forcing.q_sw_net,
        forcing.has_q_lw_down, forcing.q_lw_down,
        forcing.has_q_sh, forcing.q_sh,
        forcing.has_q_lh, forcing.q_lh,
        forcing.has_relative_humidity, forcing.relative_humidity,
        forcing.air_pressure, forcing.wind_speed,
        fields.ice_temperature, parameters.ice_substrate_layers,
        parameters.ice_substrate_top_thickness_m,
    )
    surface_temperature = _get_scalar(Tsrf, idx)
    longwave_flux, sensible_heat_flux, latent_heat_flux, rain_heat_flux =
        _resolved_nonshortwave_surface_flux_components(
            c, forcing.air_temperature, forcing.rainfall_rate, dt_seconds, surface_temperature,
            forcing.has_q_lw_down, forcing.q_lw_down, forcing.has_q_sh, forcing.q_sh,
            forcing.has_q_lh, forcing.q_lh, forcing.has_relative_humidity,
            forcing.relative_humidity, forcing.air_pressure, forcing.wind_speed,
            c.semix_z0m_ice, c.eps_ice,
        )
    _set_scalar!(workspace.net_longwave_energy, idx, _get_scalar(workspace.net_longwave_energy, idx) + longwave_flux * dt_seconds)
    _set_scalar!(workspace.absorbed_shortwave_energy, idx, _get_scalar(workspace.absorbed_shortwave_energy, idx) + energy.absorbed_shortwave * dt_seconds)
    _set_scalar!(workspace.sensible_heat_energy, idx, _get_scalar(workspace.sensible_heat_energy, idx) + sensible_heat_flux * dt_seconds)
    _set_scalar!(workspace.rain_heat_energy, idx, _get_scalar(workspace.rain_heat_energy, idx) + rain_heat_flux * dt_seconds)

    # Bare ice is a solid surface: vapour exchange carries Lᵥ+Lₘ.
    vapor_mass = latent_heat_flux * dt_seconds / (c.Lv + c.Lm)
    vapor_fluxes = _surface_vapor_fluxes(vapor_mass, latent_heat_flux)
    endogenous_melt_mass = ifelse(energy.needs_melt, energy.melt_energy_available / c.Lm, zero(dt_seconds))
    melt_mass = ifelse(
        forcing.has_prescribed_melt,
        max(forcing.prescribed_melt, zero(forcing.prescribed_melt)) * dt_seconds,
        endogenous_melt_mass,
    )
    liquid_mass = max(forcing.rainfall_rate, zero(forcing.rainfall_rate)) * dt_seconds + melt_mass

    _set_scalar!(fields.smb_ice, idx, _get_scalar(fields.smb_ice, idx) + vapor_mass - melt_mass)
    _set_scalar!(fields.melt, idx, _get_scalar(fields.melt, idx) + melt_mass)
    _set_scalar!(fields.runoff, idx, _get_scalar(fields.runoff, idx) + liquid_mass)
    _set_scalar!(fields.vapor_mass, idx, _get_scalar(fields.vapor_mass, idx) + vapor_fluxes.vapor_mass)
    _set_scalar!(fields.sublimation, idx, _get_scalar(fields.sublimation, idx) + vapor_fluxes.sublimation_mass)
    _set_scalar!(fields.latent_heat_flux_sum, idx, _get_scalar(fields.latent_heat_flux_sum, idx) + latent_heat_flux * forcing.dt_days)
    return nothing
end

Base.@propagate_inbounds function column_step_core!(
    fields,
    parameters::BESSIParameters,
    idx::Int,
    forcing::SnowpackStepForcing,
    workspace,
)
    N_storage = fields.N
    mass = fields.mass
    mass_w = fields.mass_w
    density = fields.density
    temperature = fields.temperature
    mass_base = fields.mass_base
    smb_ice = fields.smb_ice
    runoff = fields.runoff
    melt = fields.melt
    refreezing = fields.refreezing
    vapor_mass = fields.vapor_mass
    sublimation = fields.sublimation
    latent_heat_flux_sum = fields.latent_heat_flux_sum
    net_longwave_energy = workspace.net_longwave_energy
    absorbed_shortwave_energy = workspace.absorbed_shortwave_energy
    sensible_heat_energy = workspace.sensible_heat_energy
    rain_heat_energy = workspace.rain_heat_energy
    Tsrf = fields.Tsrf
    albedo_dynamic = fields.albedo
    c = parameters.c
    snow_age_days = fields.snow_age_days
    w_snow_max = fields.w_snow_max

    dt_seconds = forcing.dt_days * c.seconds_per_day
    # The geometric grid owns surface rebalancing. Applying the legacy
    # 100 kg m⁻² minimum during accumulation would merge away its gradient.
    surface_mass_min = ifelse(
        isfinite(parameters.near_surface_layer_max_thicknesses_m[1]),
        zero(parameters.mass_min),
        parameters.mass_min,
    )
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
        parameters.Ntot,
        parameters.mass_max,
        parameters.mass_split,
        surface_mass_min,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        dt_seconds,
        forcing.air_temperature,
        forcing.wind_speed,
    )
    # Snow on a previously bare surface takes the air temperature. Set it on
    # every new layer before remeshing splits the fresh snow.
    if forcing.snowfall_rate > zero(dt_seconds) && started_without_surface_snow
        @inbounds for layer_index in 1:_n_active(N_storage, idx)
            _set_layer!(temperature, layer_index, idx, forcing.air_temperature)
        end
    end
    _remesh_near_surface_layers!(
        N_storage, mass, mass_w, density, temperature,
        idx, parameters.Ntot, parameters.near_surface_layer_max_thicknesses_m, c,
    )

    has_surface_snow = _surface_has_snow(N_storage, mass, idx)
    if !has_surface_snow
        _uses_aging_albedo(c) && _set_scalar!(snow_age_days, idx, zero(dt_seconds))
        use_prescribed_albedo || _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        if parameters.ice_substrate_layers > 0
            return _step_bare_ice_substrate!(fields, parameters, idx, forcing, workspace, dt_seconds)
        end
        bare_ice_fluxes = _bare_ice_ablation_mass(c, forcing, dt_seconds)
        _set_scalar!(
            net_longwave_energy,
            idx,
            _get_scalar(net_longwave_energy, idx) + bare_ice_fluxes.longwave_flux * dt_seconds,
        )
        _set_scalar!(absorbed_shortwave_energy, idx,
            _get_scalar(absorbed_shortwave_energy, idx) + bare_ice_fluxes.absorbed_shortwave * dt_seconds)
        _set_scalar!(sensible_heat_energy, idx,
            _get_scalar(sensible_heat_energy, idx) + bare_ice_fluxes.sensible_heat_flux * dt_seconds)
        _set_scalar!(rain_heat_energy, idx,
            _get_scalar(rain_heat_energy, idx) + bare_ice_fluxes.rain_heat_flux * dt_seconds)
        rainfall_mass = max(forcing.rainfall_rate, zero(forcing.rainfall_rate)) * dt_seconds
        melt_mass = ifelse(forcing.has_prescribed_melt, max(forcing.prescribed_melt, zero(forcing.prescribed_melt)) * dt_seconds, bare_ice_fluxes.melt_mass)
        net_mass_change = bare_ice_fluxes.net_mass_change + bare_ice_fluxes.melt_mass - melt_mass
        _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) + net_mass_change)
        _set_scalar!(melt, idx, _get_scalar(melt, idx) + melt_mass)
        _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + rainfall_mass + melt_mass)
        _set_scalar!(vapor_mass, idx, _get_scalar(vapor_mass, idx) + bare_ice_fluxes.vapor_mass)
        _set_scalar!(sublimation, idx, _get_scalar(sublimation, idx) + bare_ice_fluxes.sublimation_mass)
        _set_scalar!(latent_heat_flux_sum, idx, _get_scalar(latent_heat_flux_sum, idx) + bare_ice_fluxes.latent_heat_flux * forcing.dt_days)
        return nothing
    end

    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    elseif _uses_semix_albedo(c)
        _update_semix_surface_albedo!(
            N_storage, mass, mass_w, temperature, albedo_dynamic, w_snow_max,
            idx, c, forcing.snowfall_rate, forcing,
        )
    elseif _uses_aging_albedo(c)
        _update_aging_surface_albedo_arrays!(
            N_storage, mass, temperature, albedo_dynamic, snow_age_days,
            idx, c, forcing.snowfall_rate, forcing.dt_days,
        )
    else
        _update_surface_albedo_arrays!(
            N_storage, mass, mass_w, density, temperature, albedo_dynamic,
            idx, c, forcing.dt_days,
        )
    end

    densification_tag = _densification_tag(c)
    n_liquid_water_before_energy = _copy_liquid_water_for_compaction!(
        densification_tag,
        N_storage,
        mass_w,
        workspace.liquid_water_before_energy,
        idx,
    )

    accumulation_rate = max(forcing.snowfall_rate, zero(dt_seconds)) + forcing.rainfall_rate
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
        forcing.wind_speed,
        fields.ice_temperature,
        parameters.ice_substrate_layers,
        parameters.ice_substrate_top_thickness_m,
    )
    # The implicit solve linearises Stefan--Boltzmann emission about the
    # previous surface temperature.  Record the physical, nonlinear flux at
    # the resolved surface temperature instead: the tangent approximation is
    # not an energy diagnostic when Ts changes appreciably within a substep.
    resolved_surface_temperature = _get_scalar(Tsrf, idx)
    longwave_down = ifelse(
        forcing.has_q_lw_down,
        forcing.q_lw_down,
        c.σ * c.ϵ_air * forcing.air_temperature^4,
    )
    exact_net_longwave = ifelse(
        _uses_semix_seb(c),
        c.ϵ_snow * (longwave_down - c.σ * resolved_surface_temperature^4),
        longwave_down - c.σ * c.ϵ_snow * resolved_surface_temperature^4,
    )
    _set_scalar!(
        net_longwave_energy,
        idx,
        _get_scalar(net_longwave_energy, idx) + exact_net_longwave * dt_seconds,
    )
    _set_scalar!(
        absorbed_shortwave_energy,
        idx,
        _get_scalar(absorbed_shortwave_energy, idx) + energy.absorbed_shortwave * dt_seconds,
    )
    _set_scalar!(
        sensible_heat_energy,
        idx,
        _get_scalar(sensible_heat_energy, idx) +
        (energy.sensible_heat_flux_constant - energy.sensible_heat_flux_linear * _get_scalar(Tsrf, idx)) * dt_seconds,
    )
    # The snow energy solver carries precipitation enthalpy in its effective
    # latent-heat coefficients. Record its rain component separately for the
    # energy diagnostic. Snowfall takes precedence when both phases occur,
    # matching `_diagnose_latent_heat_flux_coefficients`.
    rain_heat_flux = ifelse(
        forcing.snowfall_rate > zero(forcing.snowfall_rate),
        zero(dt_seconds),
        ifelse(
            forcing.rainfall_rate > zero(forcing.rainfall_rate),
            forcing.rainfall_rate * c.cw * (forcing.air_temperature - c.T0),
            zero(dt_seconds),
        ),
    )
    _set_scalar!(
        rain_heat_energy,
        idx,
        _get_scalar(rain_heat_energy, idx) + rain_heat_flux * dt_seconds,
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
        parameters.mass_split,
        surface_mass_min,
    )
    _set_scalar!(vapor_mass, idx, _get_scalar(vapor_mass, idx) + snow_vapor_fluxes.vapor_mass)
    _set_scalar!(sublimation, idx, _get_scalar(sublimation, idx) + snow_vapor_fluxes.sublimation_mass)
    _set_scalar!(latent_heat_flux_sum, idx, _get_scalar(latent_heat_flux_sum, idx) + snow_vapor_fluxes.latent_heat_flux * forcing.dt_days)

    endogenous_melt_mass = ifelse(energy.needs_melt, energy.melt_energy_available / c.Lm, zero(energy.melt_energy_available))
    melt_mass = ifelse(
        forcing.has_prescribed_melt,
        max(forcing.prescribed_melt, zero(forcing.prescribed_melt)) * dt_seconds,
        endogenous_melt_mass,
    )
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
        parameters.mass_split,
        surface_mass_min,
        melt_mass,
        c,
    )
    melt_reaches_ice = (melted_snow < melt_mass) & (_n_active(N_storage, idx) == 0)
    ice_melt = ifelse(melt_reaches_ice, melt_mass - melted_snow, zero(melt_mass))
    _set_scalar!(smb_ice, idx, _get_scalar(smb_ice, idx) - ice_melt)
    _set_scalar!(runoff, idx, _get_scalar(runoff, idx) + ice_melt)
    _set_scalar!(melt, idx, _get_scalar(melt, idx) + melt_mass)

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

    _apply_liquid_water_compaction!(
        densification_tag,
        N_storage,
        mass,
        mass_w,
        density,
        idx,
        workspace.liquid_water_before_energy,
        c.rho_i,
    )

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
        parameters.refreezing_correction,
    )

    _set_scalar!(refreezing, idx, _get_scalar(refreezing, idx) + refrozen_mass)

    # Melt, sublimation, compaction, and refreezing can leave a vanishingly
    # thin surface cell. Restore the prescribed near-surface geometry before
    # the next energy solve without changing any column-integrated reservoir.
    _remesh_near_surface_layers!(
        N_storage, mass, mass_w, density, temperature,
        idx, parameters.Ntot, parameters.near_surface_layer_max_thicknesses_m, c,
    )

    final_has_snow = _surface_has_snow(N_storage, mass, idx)
    if use_prescribed_albedo
        _set_prescribed_surface_albedo!(albedo_dynamic, idx, forcing)
    elseif !final_has_snow
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        _uses_aging_albedo(c) && _set_scalar!(snow_age_days, idx, zero(dt_seconds))
    end

    return nothing
end

"""Batch stepping over forcing fields with one KernelAbstractions path for CPU and GPU."""

"""
    _step_columns_single_kernel!(...)

KernelAbstractions kernel that advances each column for one forcing time
index. The absence of a dynamic time loop isolates Reactant compatibility
failures to the BESSI column-process call graph.
"""
@kernel function _step_columns_single_kernel!(
    fields,
    parameters::BESSIParameters,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing_fields,
    config,
    time_index::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        @inbounds begin
            idx = active_indices[active_idx]
            step_forcing = _step_forcing_at(forcing_fields, idx, time_index)
            column_step!(fields, parameters, idx, step_forcing, config, workspace)
        end
    end
end

"""
    _step_columns_kernel!(...)

KernelAbstractions kernel that extracts one time slice of the full forcing
and advances each column independently in-place.
"""
@kernel function _step_columns_kernel!(
    fields,
    parameters::BESSIParameters,
    workspace::ColumnarStepWorkspace,
    active_indices,
    forcing_fields,
    config,
    time_start::Int,
    time_stop::Int,
)
    active_idx = @index(Global)
    if active_idx <= length(active_indices)
        @inbounds begin
            idx = active_indices[active_idx]
            for time_index in time_start:time_stop
                step_forcing = _step_forcing_at(forcing_fields, idx, time_index)
                column_step!(fields, parameters, idx, step_forcing, config, workspace)
            end
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
    forcing,
    time_start::Int,
    time_stop::Int,
    workspace::ColumnarStepWorkspace,
    active_indices,
    config,
)
    time_stop >= time_start || return nothing
    backend = _ka_backend(state.mass)
    if time_start == time_stop
        kernel! = _step_columns_single_kernel!(backend, _step_kernel_workgroupsize(backend))
        return kernel!(
            get_fields(state),
            state.parameters,
            workspace,
            active_indices,
            _device_forcing_fields(forcing),
            config,
            time_start,
            ndrange=length(active_indices),
        )
    end
    kernel! = _step_columns_kernel!(backend, _step_kernel_workgroupsize(backend))
    return kernel!(
        get_fields(state),
        state.parameters,
        workspace,
        active_indices,
        _device_forcing_fields(forcing),
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
    forcing,
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
