"""
Simple dEBM-style diurnal shortwave partitioning helpers.
"""

@inline _debm_declination(day_of_year) = oftype(day_of_year, 23.44) * sind(oftype(day_of_year, 360) * (day_of_year - oftype(day_of_year, 79)) / oftype(day_of_year, 365))

function _debm_sunny_hours_q(latitude_deg, day_of_year, orbital_phase_deg=0.0)
    declination = _debm_declination(day_of_year + orbital_phase_deg / 360)
    cos_omega = clamp(-tand(latitude_deg) * tand(declination), -1.0, 1.0)
    omega = acos(cos_omega)
    hours = 24.0 * omega / π
    q = max(cosd(latitude_deg - declination), 0.0)
    fluxfac = q / π
    return (hours=hours, q=q, fluxfac=fluxfac)
end

function _debm_melt_window_fluxes(shortwave_down, baseline_nonshortwave_flux, latitude_deg, day_of_year)
    geometry = _debm_sunny_hours_q(latitude_deg, day_of_year, 0.0)
    sunny_fraction = clamp(geometry.hours / 24.0, 1.0 / 24.0, 1.0)
    shortwave_peak_flux = shortwave_down / sunny_fraction
    baseline_daily_flux = shortwave_down + baseline_nonshortwave_flux
    baseline_positive_daily_flux = max(baseline_daily_flux, 0.0)
    melt_window_daily_flux = max(shortwave_peak_flux + baseline_nonshortwave_flux, 0.0) * sunny_fraction
    nighttime_fraction = 1.0 - sunny_fraction
    nighttime_flux = baseline_nonshortwave_flux
    melt_period_fraction = sunny_fraction + (nighttime_flux > 0.0 ? nighttime_fraction : 0.0)
    refreezing_daily_flux = nighttime_flux < 0.0 ? -nighttime_flux * nighttime_fraction : 0.0
    return (
        baseline_daily_flux=baseline_daily_flux,
        baseline_positive_daily_flux=baseline_positive_daily_flux,
        melt_window_daily_flux=melt_window_daily_flux,
        refreezing_daily_flux=refreezing_daily_flux,
        melt_period_fraction=melt_period_fraction,
        melt_period_hours=24.0 * melt_period_fraction,
        sunny_hours=geometry.hours,
        q=geometry.q,
        fluxfac=geometry.fluxfac,
    )
end

function _diagnose_debm_diurnal_adjustment(
    c::SnowpackPhysicalConstants,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    dt_seconds,
    surface_temperature,
    surface_mass;
    q_sw_net=nothing,
    q_lw_down=nothing,
    q_sh=nothing,
    q_lh=nothing,
    latitude,
    day_of_year,
)
    shortwave_component = isnothing(q_sw_net) ? 0.0 : q_sw_net
    longwave_component = isnothing(q_lw_down) ? c.σ * (c.ϵ_air * air_temperature^4 - c.ϵ_snow * surface_temperature^4) :
        q_lw_down - c.σ * c.ϵ_snow * surface_temperature^4
    sensible_component = isnothing(q_sh) ? c.D_sh * (air_temperature - surface_temperature) : q_sh
    latent_component = isnothing(q_lh) ? 0.0 : q_lh
    rain_component = rainfall_rate * c.cw * (air_temperature - c.T0)
    baseline_nonshortwave_flux = longwave_component + sensible_component + latent_component + rain_component
    partition = _debm_melt_window_fluxes(shortwave_component, baseline_nonshortwave_flux, latitude, day_of_year)
    baseline_positive_flux = max(shortwave_component + baseline_nonshortwave_flux, 0.0)
    corrected_positive_flux = partition.melt_window_daily_flux
    extra_melt_energy = max(corrected_positive_flux - baseline_positive_flux, 0.0) * dt_seconds
    return (
        baseline_positive_flux=baseline_positive_flux,
        corrected_positive_flux=corrected_positive_flux,
        extra_melt_energy=extra_melt_energy,
        refreezing_recharge_energy=partition.refreezing_daily_flux * dt_seconds,
        refreezing_period_seconds=(1.0 - partition.melt_period_fraction) * dt_seconds,
        melt_period_seconds=partition.melt_period_fraction * dt_seconds,
        partition=partition,
    )
end

function _apply_diurnal_refreezing_recharge!(
    N_storage,
    mass,
    mass_w,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
    recharge_energy,
    refreezing_period_seconds,
)
    remaining_energy = max(recharge_energy, 0.0)
    applied_energy = 0.0

    if remaining_energy <= 0.0 || _n_active(N_storage, idx) <= 0 || !_column_has_liquid_water(N_storage, mass_w, idx)
        return (
            applied_recharge_energy=0.0,
            remaining_recharge_energy=remaining_energy,
            refreezing_period_seconds=refreezing_period_seconds,
        )
    end

    @inbounds for layer_index in 1:_n_active(N_storage, idx)
        liquid_water = _get_layer(mass_w, layer_index, idx)
        solid_mass = _get_layer(mass, layer_index, idx)
        if liquid_water <= 0.0 || solid_mass <= 0.0
            continue
        end

        layer_capacity = liquid_water * c.Lm
        applied_here = min(remaining_energy, layer_capacity)
        if applied_here > 0.0
            _set_layer!(temperature, layer_index, idx, _get_layer(temperature, layer_index, idx) - applied_here / (c.ci * solid_mass))
            remaining_energy -= applied_here
            applied_energy += applied_here
        end
        remaining_energy <= 0.0 && break
    end

    return (
        applied_recharge_energy=applied_energy,
        remaining_recharge_energy=remaining_energy,
        refreezing_period_seconds=refreezing_period_seconds,
    )
end
