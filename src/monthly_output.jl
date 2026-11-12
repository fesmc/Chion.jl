"""BESSI monthly aggregation and buffered NetCDF output."""

mutable struct MonthlyState{VT <: AbstractVector{<:AbstractFloat}}
    surface_smb::VT
    smb_ice::VT
    runoff::VT
    melt::VT
    refreezing::VT
    sublimation::VT
    latent_heat_flux::VT
    shortwave_down::VT
    shortwave_net::VT
    longwave_down::VT
    longwave_net::VT
    sensible_heat_flux::VT
    rain_heat_flux::VT
    albedo::VT
    count::Int
    days::Float64
    prev_smb_ice::VT
    prev_vapor_mass::VT
    prev_runoff::VT
    prev_melt::VT
    prev_refreezing::VT
    prev_sublimation::VT
    prev_latent_heat_flux_sum::VT
    prev_shortwave_down_sum::VT
    prev_shortwave_net_sum::VT
    prev_longwave_down_sum::VT
    prev_longwave_net_sum::VT
    prev_sensible_heat_flux_sum::VT
    prev_rain_heat_flux_sum::VT
end

mutable struct MonthlyOutputBuffer{MT <: AbstractMatrix{<:AbstractFloat}}
    surface_smb::MT
    smb_ice::MT
    runoff::MT
    melt::MT
    refreezing::MT
    sublimation::MT
    latent_heat_flux::MT
    shortwave_down::MT
    shortwave_net::MT
    longwave_down::MT
    longwave_net::MT
    sensible_heat_flux::MT
    rain_heat_flux::MT
    albedo::MT
    count::Int
end

_zero_like(v) = fill!(similar(v), zero(eltype(v)))
_prev_name(flux_sum::Symbol) = Symbol(:prev_, flux_sum)

function MonthlyState(state::BESSIState)
    return MonthlyState(
        _zero_like(state.smb_ice),
        _zero_like(state.smb_ice),
        _zero_like(state.runoff),
        _zero_like(state.melt),
        _zero_like(state.refreezing),
        _zero_like(state.sublimation),
        (_zero_like(getfield(state, flux_sum)) for (_, flux_sum) in MONTHLY_FLUX_SUMS)...,
        _zero_like(state.albedo),
        0,
        0.0,
        copy(state.smb_ice),
        copy(state.vapor_mass),
        copy(state.runoff),
        copy(state.melt),
        copy(state.refreezing),
        copy(state.sublimation),
        (copy(getfield(state, flux_sum)) for (_, flux_sum) in MONTHLY_FLUX_SUMS)...,
    )
end

function MonthlyOutputBuffer(state::BESSIState; nmonth::Integer=12)
    return MonthlyOutputBuffer(
        (similar(state.runoff, Float32, Int(nmonth), state.ncol) for _ in MONTHLY_OUTPUT_VARS)...,
        0,
    )
end

@kernel function _accumulate_monthly_precipitation!(surface_smb, snowfall_rate, rainfall_rate, time_index::Int, dt_seconds)
    idx = @index(Global)
    @inbounds surface_smb[idx] += (snowfall_rate[idx, time_index] + rainfall_rate[idx, time_index]) * dt_seconds
end

function accumulate_monthly!(monthly::MonthlyState, state::BESSIState, forcing_fields, time_index::Int, dt_days::Real=1.0)
    kernel! = _accumulate_monthly_precipitation!(_ka_backend(monthly.surface_smb))
    kernel!(
        monthly.surface_smb,
        forcing_fields.snowfall_rate,
        forcing_fields.rainfall_rate,
        time_index,
        Float64(dt_days) * state.c.seconds_per_day;
        ndrange=length(monthly.surface_smb),
    )
    monthly.albedo .+= state.albedo
    monthly.count += 1
    monthly.days += Float64(dt_days)
    return monthly
end

function finalize_monthly!(monthly::MonthlyState, state::BESSIState)
    monthly.surface_smb .-= state.runoff .- monthly.prev_runoff
    monthly.surface_smb .+= state.vapor_mass .- monthly.prev_vapor_mass
    monthly.smb_ice .= state.smb_ice .- monthly.prev_smb_ice
    monthly.runoff .= state.runoff .- monthly.prev_runoff
    monthly.melt .= state.melt .- monthly.prev_melt
    monthly.refreezing .= state.refreezing .- monthly.prev_refreezing
    monthly.sublimation .= state.sublimation .- monthly.prev_sublimation
    for (flux, flux_sum) in MONTHLY_FLUX_SUMS
        mean = getfield(monthly, flux)
        if monthly.days > 0.0
            mean .= (getfield(state, flux_sum) .- getfield(monthly, _prev_name(flux_sum))) ./ monthly.days
        else
            fill!(mean, zero(eltype(mean)))
        end
    end
    if monthly.count > 0
        monthly.albedo ./= monthly.count
    else
        monthly.albedo .= state.albedo
    end
    monthly.prev_smb_ice .= state.smb_ice
    monthly.prev_vapor_mass .= state.vapor_mass
    monthly.prev_runoff .= state.runoff
    monthly.prev_melt .= state.melt
    monthly.prev_refreezing .= state.refreezing
    monthly.prev_sublimation .= state.sublimation
    for (_, flux_sum) in MONTHLY_FLUX_SUMS
        getfield(monthly, _prev_name(flux_sum)) .= getfield(state, flux_sum)
    end
    return monthly
end

function reset_monthly!(monthly::MonthlyState)
    for field in MONTHLY_OUTPUT_VARS
        values = getfield(monthly, field)
        fill!(values, zero(eltype(values)))
    end
    monthly.count = 0
    monthly.days = 0.0
    return monthly
end

@inline function _named_fields(source, names::NTuple{N,Symbol}) where {N}
    return NamedTuple{names}(ntuple(index -> getfield(source, names[index]), Val(N)))
end

@generated function _store_named_fields_at!(
    output::NamedTuple{Names},
    source::NamedTuple{Names},
    row::Int,
    idx::Int,
) where {Names}
    assignments = [
        :(getfield(output, $(QuoteNode(name)))[row, idx] =
          getfield(source, $(QuoteNode(name)))[idx])
        for name in Names
    ]
    return :(Base.@inbounds begin
        $(assignments...)
        nothing
    end)
end

@kernel function _store_named_monthly_fields_kernel!(output, source, row::Int)
    idx = @index(Global)
    @inbounds begin
        _store_named_fields_at!(output, source, row, idx)
    end
end
function _store_named_monthly_fields!(output, source, names, row::Int, reference)
    output_fields = _named_fields(output, names)
    source_fields = _named_fields(source, names)
    kernel! = _store_named_monthly_fields_kernel!(_ka_backend(reference))
    kernel!(output_fields, source_fields, row; ndrange=length(reference))
    return output
end

function store_monthly!(output::MonthlyOutputBuffer, monthly::MonthlyState)
    row = output.count + 1
    row <= size(output.runoff, 1) || error("Monthly output buffer is full.")
    _store_named_monthly_fields!(output, monthly, MONTHLY_OUTPUT_VARS, row, monthly.runoff)
    output.count = row
    return output
end

function pdd_monthly_output_buffer(state::PDDState; nmonth::Integer)
    dims = (Int(nmonth), length(state.snowpack_swe))
    allocate(field) = similar(field, Float32, dims)
    return (
        snowpack_swe=allocate(state.snowpack_swe),
        smb_ice=allocate(state.smb_ice),
        runoff=allocate(state.runoff),
        pdd_sum=allocate(state.pdd_sum),
    )
end

function store_pdd_monthly!(output, monthly, row::Int)
    row <= size(output.runoff, 1) || error("PDD monthly output buffer is full.")
    return _store_named_monthly_fields!(output, monthly, PDD_OUTPUT_VARS, row, monthly.runoff)
end

function itm_monthly_output_buffer(state::ITMState; nmonth::Integer)
    dims = (Int(nmonth), length(state.H_snow))
    allocate(field) = similar(field, Float32, dims)
    return (; (name => allocate(getfield(state, name)) for name in ITM_OUTPUT_VARS)...)
end

function store_itm_monthly!(output, state::ITMState, row::Int)
    row <= size(output.H_snow, 1) || error("ITM monthly output buffer is full.")
    return _store_named_monthly_fields!(output, state, ITM_OUTPUT_VARS, row, state.H_snow)
end
