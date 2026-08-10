"""Fortran-compatible insolation-temperature-melt (ITM) process."""

@inline _itm_transmissivity(z, a, b) = a + b * sqrt(max(z, zero(z)))
@inline _itm_latitude_offset(c, b, lat0, lat) = c + b * (lat - lat0)
@inline _itm_planetary_albedo(albedo, a, b) = a + b * albedo
@inline _itm_snow_albedo(model, ::Nothing) = model.alb_snow_wet
@inline _itm_snow_albedo(model, melt) = ifelse(melt > model.melt_crit, model.alb_snow_wet, model.alb_snow_dry)

@inline _itm_parameters(model::ITMModel) = (
    c=model.c, trans_a=model.trans_a, trans_b=model.trans_b, itm_c=model.itm_c,
    itm_t=model.itm_t, itm_b=model.itm_b, itm_lat0=model.itm_lat0,
    H_snow_max=model.H_snow_max, Pmaxfrac=model.Pmaxfrac,
    H_snow_crit_desert=model.H_snow_crit_desert, H_snow_crit_forest=model.H_snow_crit_forest,
    melt_crit=model.melt_crit, alb_ocean=model.alb_ocean, alb_land=model.alb_land,
    alb_forest=model.alb_forest, alb_ice=model.alb_ice, alb_snow_dry=model.alb_snow_dry,
    alb_snow_wet=model.alb_snow_wet, firn_fac=model.firn_fac,
)

@inline function _itm_surface_albedo(model, z, H_ice, H_snow, PDDs, melt=nothing)
    transitional_critical_depth = model.H_snow_crit_desert +
        (model.H_snow_crit_forest - model.H_snow_crit_desert) *
        (PDDs - oftype(PDDs, 100)) / oftype(PDDs, 900)
    H_crit = ifelse(
        PDDs <= oftype(PDDs, 100),
        model.H_snow_crit_desert,
        ifelse(
            PDDs <= oftype(PDDs, 1000),
            transitional_critical_depth,
            model.H_snow_crit_forest,
        ),
    )
    depth = min(H_snow / H_crit, one(H_snow))
    capped_pdd = min(PDDs, oftype(PDDs, 1000))
    land_background = model.alb_land * (oftype(PDDs, 1000) - capped_pdd) / oftype(PDDs, 1000) +
        model.alb_forest * capped_pdd / oftype(PDDs, 1000)
    background = ifelse(
        z <= zero(z),
        model.alb_ocean,
        ifelse(H_ice == zero(H_ice), land_background, model.alb_ice),
    )
    snow_albedo = _itm_snow_albedo(model, melt)
    return background + depth * (snow_albedo - background)
end

@inline function _itm_potential_melt(model, insolation, Tair, albedo, z, latitude)
    offset = ifelse(
        abs(model.itm_lat0) < oftype(model.itm_lat0, 90),
        _itm_latitude_offset(model.itm_c, model.itm_b, model.itm_lat0, latitude),
        model.itm_c,
    )
    energy = _itm_transmissivity(z, model.trans_a, model.trans_b) * (1 - albedo) * insolation +
             offset + model.itm_t * (Tair - model.c.T0)
    return max(energy / (model.c.rho_w * model.c.Lm), zero(energy)) * model.c.seconds_per_day * 1000
end

Base.@propagate_inbounds function _itm_apply_column!(
    H_snow, alb_s, smb, smbi, melt, runoff, refreezing, Tsrf, melt_net,
    smb_cum, smbi_cum, melt_cum, runoff_cum, refreezing_cum, idx,
    model, dt, latitude, Tair, snowfall_rate, rainfall_rate, shortwave_down,
    q_sw_net, has_q_sw_net, z, H_ice, PDDs,
)
    snow = H_snow[idx]
    sf = snowfall_rate * model.c.seconds_per_day
    rf = rainfall_rate * model.c.seconds_per_day
    precipitation = sf + rf
    albedo_pre = _itm_surface_albedo(model, z, H_ice, snow, PDDs)
    snow += sf * dt
    insolation = ifelse(has_q_sw_net, q_sw_net, shortwave_down)
    potential_melt = _itm_potential_melt(model, insolation, Tair, albedo_pre, z, latitude)
    snow_melt = min(potential_melt * dt, snow)
    ice_melt = max(potential_melt * dt - snow, zero(snow))
    total_melt = snow_melt + ice_melt
    snow = max(snow - snow_melt, zero(snow))
    albedo = _itm_surface_albedo(model, z, H_ice, snow, PDDs, total_melt / dt)
    rfac = model.Pmaxfrac * sf / max(precipitation, 1e-3)
    rfac += min(one(snow), snow / 1000) * (1 - rfac)
    refrz_rain = min(rf * dt * rfac, snow)
    refrz_snow = min(snow_melt * rfac, snow - refrz_rain)
    refrz = refrz_snow + refrz_rain
    snow_to_ice = refrz
    snow -= refrz
    excess = max(snow - model.H_snow_max, zero(snow))
    snow_to_ice += excess
    snow -= excess
    runoff_step = (snow_melt - refrz_snow) + (rf * dt - refrz_rain) + ice_melt
    smb_step = precipitation * dt - runoff_step
    smbi_step = snow_to_ice + refrz - ice_melt
    melt_net_step = ifelse(H_ice > zero(H_ice), refrz - total_melt, refrz - snow_melt)
    melt_rate = total_melt / dt
    runoff_rate = runoff_step / dt
    refrz_rate = refrz / dt
    smb_rate = smb_step / dt
    smbi_rate = smbi_step / dt
    melt_net_rate = melt_net_step / dt
    H_snow[idx] = snow
    alb_s[idx] = albedo
    smb[idx] = smb_rate
    smbi[idx] = smbi_rate
    melt[idx] = melt_rate
    runoff[idx] = runoff_rate
    refreezing[idx] = refrz_rate
    melt_net[idx] = melt_net_rate
    Tsrf[idx] = ifelse(
        H_ice > zero(H_ice),
        min(model.c.T0, Tair + model.firn_fac * max(melt_net_rate, zero(melt_net_rate))),
        Tair,
    )
    smb_cum[idx] += smb_rate * dt
    smbi_cum[idx] += smbi_rate * dt
    melt_cum[idx] += melt_rate * dt
    runoff_cum[idx] += runoff_rate * dt
    refreezing_cum[idx] += refrz_rate * dt
    return nothing
end

@kernel function _itm_step_kernel!(
    H_snow, alb_s, smb, smbi, melt, runoff, refreezing, Tsrf, melt_net,
    smb_cum, smb_ice, melt_cum, runoff_cum, refreezing_cum,
    air_temperature, snowfall_rate, rainfall_rate, shortwave_down, q_sw_net, has_q_sw_net,
    latitude_deg, surface_height, ice_thickness, annual_pdd, time_index, dt, parameters, active_indices,
)
    active_idx = @index(Global)
    @inbounds begin
        idx = active_indices[active_idx]
        _itm_apply_column!(H_snow, alb_s, smb, smbi, melt, runoff, refreezing, Tsrf,
            melt_net, smb_cum, smb_ice, melt_cum, runoff_cum, refreezing_cum,
            idx, parameters, dt,
            latitude_deg[idx, time_index], air_temperature[idx, time_index], snowfall_rate[idx, time_index],
            rainfall_rate[idx, time_index], shortwave_down[idx, time_index], q_sw_net[idx, time_index],
            has_q_sw_net[idx, time_index], surface_height[idx, time_index], ice_thickness[idx, time_index],
            annual_pdd[idx, time_index])
    end
end

function _itm_step_arrays!(state::ITMState, forcing, time_index::Int, model::ITMModel, active_indices)
    isempty(active_indices) && return nothing
    kernel! = _itm_step_kernel!(_ka_backend(state.H_snow))
    event = kernel!(state.H_snow, state.alb_s, state.smb, state.smbi, state.melt,
        state.runoff, state.refreezing, state.Tsrf, state.melt_net, state.smb_cum,
        state.smb_ice, state.melt_cum, state.runoff_cum, state.refreezing_cum,
        forcing.air_temperature, forcing.snowfall_rate, forcing.rainfall_rate, forcing.shortwave_down,
        forcing.q_sw_net, forcing.has_q_sw_net, forcing.latitude_deg, forcing.surface_height,
        forcing.ice_thickness, forcing.annual_pdd, time_index, _step_dt(forcing.dt_days, time_index),
        _itm_parameters(model), active_indices; ndrange=length(active_indices))
    _wait_kernel(event)
    return nothing
end
