"""
Surface albedo state and update rules for array-backed snow states.
"""

@inline function _semix_grain_size(temperature, snowfall_rate, c)
    cold = min(zero(temperature), temperature - (c.T0 - c.semix_dT_age))
    f_age_temperature = exp(c.semix_f_age_t * cold) + exp(cold)
    snow_per_day = max(snowfall_rate, zero(snowfall_rate)) * c.seconds_per_day
    f_p = f_age_temperature * (c.semix_snow_0 / max(oftype(snow_per_day, 1e-20), snow_per_day))^c.semix_snow_1
    f_age = one(f_p) - log(one(f_p) + f_p) / f_p
    return c.semix_snow_grain_fresh + (c.semix_snow_grain_old - c.semix_snow_grain_fresh) * f_age
end

@inline function _semix_dust_concentration(dust_deposition, snowfall_rate, snow_swe, snow_swe_max, c)
    dust = dust_deposition / max(oftype(snowfall_rate, 1e-7), snowfall_rate)
    melt_factor = ifelse(snow_swe > oftype(snow_swe, 1),
        clamp(one(snow_swe) + (snow_swe_max - snow_swe) / c.semix_w_snow_dust, one(snow_swe), oftype(snow_swe, 5)), one(snow_swe))
    dust = min(melt_factor * c.semix_dust_con_scale * dust, oftype(dust, 1000e-6))
    return ifelse(dust < oftype(dust, 1e-15), zero(dust), dust)
end

@inline function _semix_broadband_albedo(avd, and, avf, anf, cloud, c)
    cl = clamp(cloud, zero(cloud), one(cloud)); fv = clamp(c.semix_frac_vu, zero(c.semix_frac_vu), one(c.semix_frac_vu))
    return (one(cl) - cl) * (fv * avd + (one(fv) - fv) * and) + cl * (fv * avf + (one(fv) - fv) * anf)
end

@inline function _semix_ww_bands(grain, dust, coszm, c)
    fcos = max(zero(coszm), oftype(coszm, 0.5) * (oftype(coszm, 3) / (one(coszm) + oftype(coszm, 2) * coszm) - one(coszm)))
    d = min(dust * oftype(dust, 1e6), oftype(dust, 999))
    new_vis, aged_vis = ifelse(d <= oftype(d, 1.0001), (zero(d), zero(d)),
        ifelse(d < oftype(d, 10), (oftype(d, 0.02) * log(d)/log(oftype(d, 10)), oftype(d, 0.05) * log(d)/log(oftype(d, 10))),
        ifelse(d < oftype(d, 100), (oftype(d, 0.02) + oftype(d, 0.08)*(log(d)-log(oftype(d,10)))/(log(oftype(d,100))-log(oftype(d,10))), oftype(d,0.05)+oftype(d,0.10)*(log(d)-log(oftype(d,10)))/(log(oftype(d,100))-log(oftype(d,10)))),
        (oftype(d,0.10)+oftype(d,0.20)*(log(d)-log(oftype(d,100)))/(log(oftype(d,1000))-log(oftype(d,100))), oftype(d,0.15)+oftype(d,0.15)*(log(d)-log(oftype(d,100)))/(log(oftype(d,1000))-log(oftype(d,100)))))))
    fage = log10(one(grain) + (grain - c.semix_snow_grain_fresh) / oftype(grain, 200)) / log10(one(grain) + (c.semix_snow_grain_old-c.semix_snow_grain_fresh)/oftype(grain,200))
    avf = c.semix_alb_snow_vis_new - fage * (c.semix_d_alb_age_vis + aged_vis) - new_vis
    anf = c.semix_alb_snow_nir_new - fage * (c.semix_d_alb_age_nir + aged_vis / oftype(aged_vis, 2)) - new_vis / oftype(new_vis, 2)
    return avf + oftype(avf, .4)*fcos*(one(avf)-avf), anf + oftype(anf,.4)*fcos*(one(anf)-anf), avf, anf
end

@inline function _semix_dang_bands(grain, dust, coszm, z_sur_std, c)
    r0 = oftype(grain, 100); rn = log10(grain / r0)
    x = ifelse(dust > oftype(dust, 1e-8), log10(dust * oftype(dust, 1e6)), zero(dust))
    roughness = c.semix_k_sigma_orog * tanh(z_sur_std / c.semix_sigma_orog_crit)
    function visible(r, direct)
        rn_local = log10(r / r0)
        base = ifelse(direct, oftype(r, .9849) - oftype(r, .0215)*rn_local - oftype(r, .0132)*rn_local^2,
                      oftype(r, .9856) - oftype(r, .0202)*rn_local - oftype(r, .0125)*rn_local^2)
        f = ifelse(direct, oftype(r,155)+oftype(r,17.15)*x+oftype(r,.27)*x^2, oftype(r,152)+oftype(r,15.92)*x-oftype(r,.39)*x^2)
        h = dust / f / oftype(r,1e-6) * (r/r0)^oftype(r,.73)
        p = log10(max(h, oftype(h, 1e-30)))
        dark = ifelse(dust > oftype(dust,1e-8), oftype(r,10)^(-oftype(r,.05)*p^2 + ifelse(direct,oftype(r,.525),oftype(r,.514))*p - ifelse(direct,oftype(r,.893),oftype(r,.890))), zero(r))
        return min(one(r), base - dark) - roughness
    end
    avf = visible(grain, false)
    anf = min(one(grain), oftype(grain,.7493) - oftype(grain,.182)*rn - oftype(grain,.0388)*rn^2) - roughness
    rv = grain * (one(grain)+oftype(grain,.781)*(coszm-oftype(grain,.65))^2)
    rn_direct = log10((grain*(one(grain)+oftype(grain,.791)*(coszm-oftype(grain,.65))^2))/r0)
    return visible(rv, true), min(one(grain),oftype(grain,.6596)-oftype(grain,.1927)*rn_direct-oftype(grain,.0229)*rn_direct^2)-roughness, avf, anf
end

@inline function _semix_daily_coszm(latitude_deg, solar_longitude_deg)
    (!isfinite(latitude_deg) || !isfinite(solar_longitude_deg)) && return zero(latitude_deg)
    declination = asin(sin(oftype(latitude_deg, π) / oftype(latitude_deg,180) * oftype(latitude_deg,23.44)) * sin(solar_longitude_deg * oftype(latitude_deg,π)/oftype(latitude_deg,180)))
    latitude = latitude_deg * oftype(latitude_deg, π) / oftype(latitude_deg,180)
    h0 = acos(clamp(-tan(latitude)*tan(declination), -one(latitude), one(latitude)))
    return ifelse(h0 <= zero(h0), zero(h0), max(zero(h0), (h0*sin(latitude)*sin(declination)+cos(latitude)*cos(declination)*sin(h0))/h0))
end

@inline function _semix_column_swe(N_storage, mass, mass_w, idx)
    total = zero(_get_layer(mass, 1, idx))
    @inbounds for layer_index in 1:_n_active(N_storage, idx)
        total += max(_get_layer(mass, layer_index, idx), zero(total)) + max(_get_layer(mass_w, layer_index, idx), zero(total))
    end
    return total
end

@inline function _update_semix_surface_albedo!(N_storage, mass, mass_w, temperature, albedo, w_snow_max, idx, c, snowfall_rate, forcing)
    snow_swe = _semix_column_swe(N_storage, mass, mass_w, idx)
    previous_max = _get_scalar(w_snow_max, idx)
    seasonal_max = forcing.day_of_year <= forcing.dt_days ? snow_swe : max(previous_max, snow_swe)
    _set_scalar!(w_snow_max, idx, seasonal_max)
    grain = _semix_grain_size(_get_layer(temperature, 1, idx), snowfall_rate, c)
    coszm = ifelse(forcing.has_coszm, clamp(forcing.coszm, zero(snow_swe), one(snow_swe)),
                   _semix_daily_coszm(forcing.latitude_deg, forcing.solar_longitude_deg))
    cloud = ifelse(forcing.has_cloud, forcing.cloud, zero(snow_swe))
    z_sur_std = ifelse(forcing.has_z_sur_std, max(forcing.z_sur_std, zero(snow_swe)), zero(snow_swe))
    dust_deposition = ifelse(forcing.has_dust_deposition, max(forcing.dust_deposition, zero(snow_swe)), zero(snow_swe))
    dust = _semix_dust_concentration(dust_deposition, snowfall_rate, snow_swe, seasonal_max, c)
    bands = c.semix_snow_albedo == SEMIX_ALBEDO_DANG ?
            _semix_dang_bands(grain, dust, coszm, z_sur_std, c) :
            _semix_ww_bands(grain, dust, coszm, c)
    updated = clamp(_semix_broadband_albedo(bands..., cloud, c), zero(snow_swe), one(snow_swe))
    _set_scalar!(albedo, idx, updated)
    return updated
end

"""
    _constant_surface_albedo(N_storage, mass, temperature, idx, c)

Return the constant-scheme surface albedo for column `idx`, switching between
ice, dry snow, and wet snow according to surface state.
"""
@inline function _constant_surface_albedo(
    N_storage,
    mass,
    temperature,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
        return c.alpha_ice
    end
    return _get_layer(temperature, 1, idx) >= c.T0 ? c.alpha_wet : c.alpha_dry
end

"""
    _surface_liquid_water_content(N_storage, mass, mass_w, density, idx, c)

Estimate volumetric liquid-water content in the surface layer of column `idx`.
Returns zero when the layer is empty, ice-dense, or has no pore volume.
"""
@inline function _surface_liquid_water_content(
    N_storage,
    mass,
    mass_w,
    density,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_TINY
        return zero(eltype(mass))
    end

    surface_density = _get_layer(density, 1, idx)
    if surface_density <= EPS_TINY || surface_density >= c.rho_i - EPS_TINY
        return zero(surface_density)
    end

    pore_volume = _get_layer(mass, 1, idx) / surface_density - _get_layer(mass, 1, idx) / c.rho_i
    if pore_volume <= EPS_TINY
        return zero(pore_volume)
    end

    return max(_get_layer(mass_w, 1, idx), zero(eltype(mass))) / c.rho_w / pore_volume
end

"""
    _refresh_dynamic_albedo_from_snowfall!(albedo_dynamic, idx, c, snowfall_mass)

Refresh the dynamic surface albedo after snowfall on column `idx`. Mutates
`albedo_dynamic[idx]` and returns the updated albedo.
"""
function _refresh_dynamic_albedo_from_snowfall!(
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    snowfall_mass,
)
    snowfall_mass <= EPS_TINY && return _get_scalar(albedo_dynamic, idx)

    if _uses_constant_albedo(c)
        _set_scalar!(albedo_dynamic, idx, c.alpha_dry)
        return _get_scalar(albedo_dynamic, idx)
    elseif _uses_aging_albedo(c)
        _set_scalar!(albedo_dynamic, idx, c.alpha_dry)
        return _get_scalar(albedo_dynamic, idx)
    end

    updated = min(
        c.alpha_dry,
        _get_scalar(albedo_dynamic, idx) +
        (c.alpha_dry - c.alpha_wet) * (one(snowfall_mass) - exp(-snowfall_mass / oftype(snowfall_mass, 3))),
    )
    _set_scalar!(albedo_dynamic, idx, updated)
    return updated
end

"""
    _update_aging_surface_albedo_arrays!(N_storage, mass, temperature,
        albedo_dynamic, snow_age_days, idx, c, snowfall_rate, dt_days)

Update snow albedo from the time elapsed since the latest snowfall event.
Snowfall resets the age to zero and the albedo to `alpha_dry`.
Otherwise snow age advances by `dt_days` and albedo decays exponentially
toward `alpha_wet`, with separate cold- and melting-surface timescales.
Bare columns use `alpha_ice` and carry zero snow age.
"""
function _update_aging_surface_albedo_arrays!(
    N_storage,
    mass,
    temperature,
    albedo_dynamic,
    snow_age_days,
    idx::Int,
    c::SnowpackPhysicalConstants,
    snowfall_rate,
    dt_days,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
        _set_scalar!(snow_age_days, idx, zero(dt_days))
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        return c.alpha_ice
    end

    if snowfall_rate > zero(snowfall_rate)
        _set_scalar!(snow_age_days, idx, zero(dt_days))
        _set_scalar!(albedo_dynamic, idx, c.alpha_dry)
        return c.alpha_dry
    end
    age_days = max(_get_scalar(snow_age_days, idx), zero(dt_days)) + dt_days
    _set_scalar!(snow_age_days, idx, age_days)

    surface_temperature = _get_layer(temperature, 1, idx)
    timescale_days = surface_temperature >= c.T0 ?
        c.aging_melting_timescale_days :
        c.aging_cold_timescale_days
    previous_albedo = clamp(
        _get_scalar(albedo_dynamic, idx),
        c.alpha_wet,
        c.alpha_dry,
    )
    updated_albedo = c.alpha_wet +
                     (previous_albedo - c.alpha_wet) *
                     exp(-dt_days / timescale_days)
    updated_albedo = clamp(updated_albedo, c.alpha_wet, c.alpha_dry)
    _set_scalar!(albedo_dynamic, idx, updated_albedo)
    return updated_albedo
end

"""
    _update_surface_albedo_arrays!(N_storage, mass, mass_w, density, temperature, albedo_dynamic, idx, c, dt_days)

Update the diagnosed surface albedo for column `idx` in-place using either the
constant or dynamic albedo scheme. Dynamic snow aging is scaled by `dt_days`.
"""
function _update_surface_albedo_arrays!(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
    dt_days,
)
    if _n_active(N_storage, idx) <= 0 || _get_layer(mass, 1, idx) <= EPS_EMPTY_LAYER
        _set_scalar!(albedo_dynamic, idx, c.alpha_ice)
        return c.alpha_ice
    end

    if _uses_constant_albedo(c)
        updated = _constant_surface_albedo(N_storage, mass, temperature, idx, c)
        _set_scalar!(albedo_dynamic, idx, updated)
        return updated
    end

    previous_albedo = clamp(_get_scalar(albedo_dynamic, idx), c.alpha_wet, c.alpha_dry)
    surface_temperature = _get_layer(temperature, 1, idx)
    aging_decrement = (
        (surface_temperature - c.T0) * oftype(surface_temperature, 1.35e-3) +
        oftype(surface_temperature, 0.0278)
    ) * oftype(surface_temperature, dt_days)
    updated_albedo = min(
        previous_albedo,
        previous_albedo - aging_decrement,
    )
    updated_albedo = max(updated_albedo, c.alpha_wet)

    liquid_water_content = _surface_liquid_water_content(N_storage, mass, mass_w, density, idx, c)
    if liquid_water_content > zero(liquid_water_content) && c.max_lwc_albedo > EPS_TINY
        wet_adjusted_albedo = updated_albedo - (
            updated_albedo - c.alpha_wet
        ) * (liquid_water_content / c.max_lwc_albedo)
        updated_albedo = max(c.alpha_wet, min(updated_albedo, wet_adjusted_albedo))
    end

    updated_albedo = clamp(updated_albedo, c.alpha_wet, c.alpha_dry)
    _set_scalar!(albedo_dynamic, idx, updated_albedo)
    return updated_albedo
end

"""
    update_surface_albedo!(state, idx, dt_days=1)

Update the surface albedo of column `idx` in `state` and return the new
albedo. `dt_days` controls the elapsed time applied to dynamic snow aging.
"""
function update_surface_albedo!(state, idx::Int, dt_days=one(eltype(state.mass)))
    if _uses_aging_albedo(state.c)
        return _update_aging_surface_albedo_arrays!(
            state.N,
            state.mass,
            state.temperature,
            state.albedo,
            state.snow_age_days,
            idx,
            state.c,
            zero(dt_days),
            dt_days,
        )
    end
    return _update_surface_albedo_arrays!(
        state.N,
        state.mass,
        state.mass_w,
        state.density,
        state.temperature,
        state.albedo,
        idx,
        state.c,
        dt_days,
    )
end
