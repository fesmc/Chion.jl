"""
Positive-degree-day mass-balance process.

This is a bulk, capped one-layer snow model synchronized with the Fortran
`fesmc/chion` implementation. Snow retained in the reservoir is distinct from
the ice-facing SMB flux.
"""

@inline function _pdd_normal_cdf(value)
    half = oftype(value, 0.5)
    sqrt_two = sqrt(oftype(value, 2))
    return half * erfc(-value / sqrt_two)
end

@inline function _pdd_expected_positive_temperature(mean_temperature_c, sigma)
    z = mean_temperature_c / sigma
    inv_sqrt_two_pi = oftype(z, 0.398942280401432678)
    pdf = inv_sqrt_two_pi * exp(-oftype(z, 0.5) * z * z)
    return sigma * pdf + mean_temperature_c * _pdd_normal_cdf(z)
end

@inline function _pdd_degree_days(
    air_temperature,
    dt_days,
    freezing_temperature,
    temperature_sigma,
    ::SimplePDD,
)
    mean_temperature_c = air_temperature - freezing_temperature
    return max(mean_temperature_c, zero(mean_temperature_c)) * dt_days
end

@inline function _pdd_degree_days(
    air_temperature,
    dt_days,
    freezing_temperature,
    temperature_sigma,
    ::PISMPDD,
)
    mean_temperature_c = air_temperature - freezing_temperature
    return dt_days * _pdd_expected_positive_temperature(mean_temperature_c, temperature_sigma)
end

@inline function _pdd_step_mass(rate, dt_days, seconds_per_day)
    return max(rate, zero(rate)) * dt_days * seconds_per_day
end

"""
    _pdd_apply_column!(..., snowfall, rainfall, pdd, model parameters...)

Apply one PDD step to one column. The budget follows the Fortran reference:

  * snow melt is limited by the available snow;
  * refreezing is limited by the snow remaining after melt;
  * refrozen water becomes superimposed ice instead of re-entering the snow;
  * snow above `H_snow_max` is converted to ice; and
  * `smb_ice` records only snow-to-ice conversion minus ice melt.

Consequently, every step closes as
`snowfall + rainfall == Δsnowpack_swe + Δsmb_ice + Δrunoff`.
"""
Base.@propagate_inbounds function _pdd_apply_column!(
    snowpack_swe,
    smb_ice,
    runoff,
    pdd_sum,
    idx::Int,
    snowfall,
    rainfall,
    pdd,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
    H_snow_max,
)
    pdd_sum[idx] += pdd

    available_snow = snowpack_swe[idx] + snowfall
    snow_melt = min(available_snow, ddf_snow * pdd)
    positive_degree_factor = ddf_snow > zero(ddf_snow)
    safe_degree_factor = ifelse(positive_degree_factor, ddf_snow, one(ddf_snow))
    remaining_pdd = ifelse(
        positive_degree_factor,
        max(pdd - snow_melt / safe_degree_factor, zero(pdd)),
        zero(pdd),
    )
    ice_melt = ddf_ice * remaining_pdd

    remaining_snow = available_snow - snow_melt
    refrozen = min(refreezing_fraction * remaining_snow, snow_melt)
    excess_snow = max(remaining_snow - H_snow_max, zero(remaining_snow))
    snow_to_ice = refrozen + excess_snow

    snowpack_swe[idx] = min(remaining_snow, H_snow_max)
    smb_ice[idx] += snow_to_ice - ice_melt
    runoff[idx] += rainfall + snow_melt - refrozen + ice_melt
    return nothing
end

@kernel function _pdd_step_kernel!(
    fields,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    time_index::Int,
    dt_days,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
    temperature_sigma,
    H_snow_max,
    pdd_method,
    freezing_temperature,
    seconds_per_day,
    active_indices,
)
    active_idx = @index(Global)
    @inbounds begin
        idx = active_indices[active_idx]
        snowfall = _pdd_step_mass(
            snowfall_rate[idx, time_index],
            dt_days,
            seconds_per_day,
        )
        rainfall = _pdd_step_mass(
            rainfall_rate[idx, time_index],
            dt_days,
            seconds_per_day,
        )
        pdd = _pdd_degree_days(
            air_temperature[idx, time_index],
            dt_days,
            freezing_temperature,
            temperature_sigma,
            pdd_method,
        )
        _pdd_apply_column!(
            fields.snowpack_swe,
            fields.smb_ice,
            fields.runoff,
            fields.pdd_sum,
            idx,
            snowfall,
            rainfall,
            pdd,
            ddf_snow,
            ddf_ice,
            refreezing_fraction,
            H_snow_max,
        )
    end
end

function _pdd_step_arrays!(
    state::PDDState,
    forcing,
    time_index::Int,
    params::PDDParams,
    active_indices,
)
    launch!(_pdd_step_kernel!, state.snowpack_swe, length(active_indices),
        get_fields(state),
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        time_index,
        _step_dt(forcing.dt_days, time_index),
        params.ddf_snow,
        params.ddf_ice,
        params.refreezing_fraction,
        params.temperature_sigma,
        params.H_snow_max,
        params.method,
        params.T0,
        params.seconds_per_day,
        active_indices,
    )
    return nothing
end
