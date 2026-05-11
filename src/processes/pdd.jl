"""
Positive-degree-day mass-balance process.

This is a bulk snow/ice SMB model: snowfall builds a snow reservoir, positive
degree days melt snow first, and any remaining melt potential ablates ice.
"""

@inline _pdd_positive_degree_days(air_temperature, dt_days) =
    max(air_temperature - 273.15, zero(air_temperature)) * dt_days

@inline _pdd_step_mass(rate, dt_days) = max(rate, zero(rate)) * dt_days * 86_400.0
@inline _normal_pdf(x) = 0.3989422804014327 * exp(-0.5 * x * x)

@inline function _normal_cdf(x)
    t = inv(1.0 + 0.2316419 * abs(x))
    poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    cdf = 1.0 - _normal_pdf(x) * poly
    return x < 0 ? 1.0 - cdf : cdf
end

@inline function _expected_positive_temperature(mean_temperature_c, sigma)
    z = mean_temperature_c / sigma
    return sigma * _normal_pdf(z) + mean_temperature_c * _normal_cdf(z)
end

@inline function _pism_expected_positive_degree_days(air_temperature, dt_days, temperature_sigma)
    return dt_days * _expected_positive_temperature(air_temperature - 273.15, temperature_sigma)
end

function _pdd_apply_column_pdd!(
    snowpack_swe::AbstractVector,
    smb_ice::AbstractVector,
    runoff::AbstractVector,
    pdd_sum::AbstractVector,
    idx::Int,
    snowfall,
    rainfall,
    pdd,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
)
    pdd_sum[idx] += pdd

    available_snow = snowpack_swe[idx] + snowfall
    snow_melt_potential = ddf_snow * pdd
    snow_melt = min(available_snow, snow_melt_potential)
    remaining_pdd = ddf_snow > zero(ddf_snow) ? max(pdd - snow_melt / ddf_snow, zero(pdd)) : zero(pdd)
    ice_melt = ddf_ice * remaining_pdd
    refrozen = refreezing_fraction * snow_melt

    snowpack_swe[idx] = available_snow - snow_melt + refrozen
    smb_ice[idx] += snowfall - snow_melt + refrozen - ice_melt
    runoff[idx] += rainfall + snow_melt - refrozen + ice_melt
    return nothing
end

function _pdd_step_column!(
    snowpack_swe::AbstractVector,
    smb_ice::AbstractVector,
    runoff::AbstractVector,
    pdd_sum::AbstractVector,
    idx::Int,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    dt_days,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
)
    snowfall = _pdd_step_mass(snowfall_rate, dt_days)
    rainfall = _pdd_step_mass(rainfall_rate, dt_days)
    pdd = _pdd_positive_degree_days(air_temperature, dt_days)
    return _pdd_apply_column_pdd!(
        snowpack_swe,
        smb_ice,
        runoff,
        pdd_sum,
        idx,
        snowfall,
        rainfall,
        pdd,
        ddf_snow,
        ddf_ice,
        refreezing_fraction,
    )
end

function pdd_step!(model::PDDModel, state::PDDState, forcing::SnowpackForcing, time_index::Int)
    ncol = ncols(model.grid)
    size(forcing.air_temperature, 1) == ncol || error("Forcing column count must match the PDD model column count.")
    if forcing_step_kind(forcing, time_index) === :monthly
        scratch = (
            a=Vector{Float64}(undef, ncol),
            b=Vector{Float64}(undef, ncol),
            c=Vector{Float64}(undef, ncol),
            d=Vector{Float64}(undef, ncol),
            e=Vector{Float64}(undef, ncol),
            f=Vector{Float64}(undef, ncol),
        )
        pdd_monthly_step!(
            state.snowpack_swe,
            state.smb_ice,
            state.runoff,
            state.pdd_sum,
            forcing,
            time_index,
            model,
            scratch,
        )
        return nothing
    end
    pdd_step!(
        state.snowpack_swe,
        state.smb_ice,
        state.runoff,
        state.pdd_sum,
        forcing,
        time_index,
        model.ddf_snow,
        model.ddf_ice,
        model.refreezing_fraction,
    )
    return nothing
end

@kernel function _pdd_step_kernel!(
    snowpack_swe,
    smb_ice,
    runoff,
    pdd_sum,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    time_index::Int,
    dt_days,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
)
    idx = @index(Global)
    if idx <= length(snowpack_swe)
        _pdd_step_column!(
            snowpack_swe,
            smb_ice,
            runoff,
            pdd_sum,
            idx,
            air_temperature[idx, time_index],
            snowfall_rate[idx, time_index],
            rainfall_rate[idx, time_index],
            dt_days,
            ddf_snow,
            ddf_ice,
            refreezing_fraction,
        )
    end
end

@kernel function _pdd_monthly_step_kernel!(
    snowpack_swe,
    smb_ice,
    runoff,
    pdd_sum,
    air_temperature,
    snowfall_rate,
    rainfall_rate,
    time_index::Int,
    dt_days,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
    temperature_sigma,
)
    idx = @index(Global)
    if idx <= length(snowpack_swe)
        _pdd_apply_column_pdd!(
            snowpack_swe,
            smb_ice,
            runoff,
            pdd_sum,
            idx,
            _pdd_step_mass(snowfall_rate[idx, time_index], dt_days),
            _pdd_step_mass(rainfall_rate[idx, time_index], dt_days),
            _pism_expected_positive_degree_days(air_temperature[idx, time_index], dt_days, temperature_sigma),
            ddf_snow,
            ddf_ice,
            refreezing_fraction,
        )
    end
end

# CPU vectorised path — no per-column loop, operates on whole arrays at once.
# scratch is a NamedTuple of 6 pre-allocated Vector{Float64} buffers of length ncol.
function pdd_step!(
    snowpack_swe::Vector{Float64},
    smb_ice::Vector{Float64},
    runoff::Vector{Float64},
    pdd_sum::Vector{Float64},
    forcing::SnowpackForcing,
    time_index::Int,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
    scratch,
)
    dt = _step_dt(forcing.dt_days, time_index)
    T  = @view forcing.air_temperature[:, time_index]
    sf = @view forcing.snowfall_rate[:, time_index]
    rf = @view forcing.rainfall_rate[:, time_index]

    snowfall       = scratch.a
    rainfall       = scratch.b
    pdd            = scratch.c
    available_snow = scratch.d
    snow_melt      = scratch.e
    remaining_pdd  = scratch.f

    scale = dt * 86_400.0
    @. snowfall       = max(sf, 0.0) * scale
    @. rainfall       = max(rf, 0.0) * scale
    @. pdd            = max(T - 273.15, 0.0) * dt
    @. pdd_sum       += pdd
    @. available_snow = snowpack_swe + snowfall
    @. snow_melt      = min(available_snow, ddf_snow * pdd)
    @. remaining_pdd  = max(pdd - snow_melt / ddf_snow, 0.0)

    @inbounds for i in eachindex(snowpack_swe)
        refrozen        = refreezing_fraction * snow_melt[i]
        ice_melt        = ddf_ice * remaining_pdd[i]
        snowpack_swe[i] = available_snow[i] - snow_melt[i] + refrozen
        smb_ice[i]     += snowfall[i] - snow_melt[i] + refrozen - ice_melt
        runoff[i]      += rainfall[i] + snow_melt[i] - refrozen + ice_melt
    end
    return nothing
end

function pdd_monthly_step!(
    snowpack_swe::Vector{Float64},
    smb_ice::Vector{Float64},
    runoff::Vector{Float64},
    pdd_sum::Vector{Float64},
    forcing::SnowpackForcing,
    time_index::Int,
    model::PDDModel,
    scratch,
)
    dt = _step_dt(forcing.dt_days, time_index)
    T  = @view forcing.air_temperature[:, time_index]
    sf = @view forcing.snowfall_rate[:, time_index]
    rf = @view forcing.rainfall_rate[:, time_index]

    snowfall       = scratch.a
    rainfall       = scratch.b
    pdd            = scratch.c
    available_snow = scratch.d
    snow_melt      = scratch.e
    remaining_pdd  = scratch.f

    scale = dt * 86_400.0
    sigma = model.monthly_method.temperature_sigma
    @. snowfall       = max(sf, 0.0) * scale
    @. rainfall       = max(rf, 0.0) * scale
    @. pdd            = _pism_expected_positive_degree_days(T, dt, sigma)
    @. pdd_sum       += pdd
    @. available_snow = snowpack_swe + snowfall
    @. snow_melt      = min(available_snow, model.ddf_snow * pdd)
    @. remaining_pdd  = max(pdd - snow_melt / model.ddf_snow, 0.0)

    @inbounds for i in eachindex(snowpack_swe)
        refrozen        = model.refreezing_fraction * snow_melt[i]
        ice_melt        = model.ddf_ice * remaining_pdd[i]
        snowpack_swe[i] = available_snow[i] - snow_melt[i] + refrozen
        smb_ice[i]     += snowfall[i] - snow_melt[i] + refrozen - ice_melt
        runoff[i]      += rainfall[i] + snow_melt[i] - refrozen + ice_melt
    end
    return nothing
end

function pdd_step!(
    snowpack_swe::Vector{Float64},
    smb_ice::Vector{Float64},
    runoff::Vector{Float64},
    pdd_sum::Vector{Float64},
    forcing::SnowpackForcing,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
    scratch,
)
    for time_index in 1:_step_time_count(forcing)
        pdd_step!(snowpack_swe, smb_ice, runoff, pdd_sum, forcing, time_index, ddf_snow, ddf_ice, refreezing_fraction, scratch)
    end
    return nothing
end

# GPU / fallback path — KernelAbstractions kernel, one launch per time step.
function pdd_step!(
    snowpack_swe::AbstractVector,
    smb_ice::AbstractVector,
    runoff::AbstractVector,
    pdd_sum::AbstractVector,
    forcing::SnowpackForcing,
    time_index::Int,
    ddf_snow,
    ddf_ice,
    refreezing_fraction,
)
    size(forcing.air_temperature, 1) == length(snowpack_swe) || error("Forcing column count must match the PDD state column count.")
    length(pdd_sum) == length(snowpack_swe) || error("PDD diagnostic length must match the PDD state column count.")
    kernel! = _pdd_step_kernel!(_ka_backend(snowpack_swe))
    event = kernel!(
        snowpack_swe,
        smb_ice,
        runoff,
        pdd_sum,
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        time_index,
        _step_dt(forcing.dt_days, time_index),
        ddf_snow,
        ddf_ice,
        refreezing_fraction;
        ndrange=length(snowpack_swe),
    )
    _wait_kernel(event)
    return nothing
end

function pdd_monthly_step!(
    snowpack_swe::AbstractVector,
    smb_ice::AbstractVector,
    runoff::AbstractVector,
    pdd_sum::AbstractVector,
    forcing::SnowpackForcing,
    time_index::Int,
    model::PDDModel,
)
    size(forcing.air_temperature, 1) == length(snowpack_swe) || error("Forcing column count must match the PDD state column count.")
    length(pdd_sum) == length(snowpack_swe) || error("PDD diagnostic length must match the PDD state column count.")
    kernel! = _pdd_monthly_step_kernel!(_ka_backend(snowpack_swe))
    event = kernel!(
        snowpack_swe,
        smb_ice,
        runoff,
        pdd_sum,
        forcing.air_temperature,
        forcing.snowfall_rate,
        forcing.rainfall_rate,
        time_index,
        _step_dt(forcing.dt_days, time_index),
        model.ddf_snow,
        model.ddf_ice,
        model.refreezing_fraction,
        model.monthly_method.temperature_sigma;
        ndrange=length(snowpack_swe),
    )
    _wait_kernel(event)
    return nothing
end

function pdd_step!(model::PDDModel, state::PDDState, forcing::SnowpackForcing)
    for time_index in 1:_step_time_count(forcing)
        pdd_step!(model, state, forcing, time_index)
    end
    return nothing
end
