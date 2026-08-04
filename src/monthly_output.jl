"""BESSI monthly aggregation and buffered NetCDF output."""

mutable struct MonthlyState{VT <: AbstractVector{<:AbstractFloat}}
    smb_ice::VT
    runoff::VT
    melt::VT
    refreezing::VT
    sublimation::VT
    latent_heat_flux::VT
    albedo::VT
    count::Int
    days::Float64
    prev_smb_ice::VT
    prev_runoff::VT
    prev_melt::VT
    prev_refreezing::VT
    prev_sublimation::VT
    prev_latent_heat_flux_sum::VT
end

mutable struct MonthlyOutputBuffer{MT <: AbstractMatrix{<:AbstractFloat}}
    smb_ice::MT
    runoff::MT
    melt::MT
    refreezing::MT
    sublimation::MT
    latent_heat_flux::MT
    albedo::MT
    count::Int
end

_zero_like(v) = fill!(similar(v), zero(eltype(v)))

function MonthlyState(state::BESSIState)
    return MonthlyState(
        _zero_like(state.smb_ice),
        _zero_like(state.runoff),
        _zero_like(state.melt),
        _zero_like(state.refreezing),
        _zero_like(state.sublimation),
        _zero_like(state.latent_heat_flux_sum),
        _zero_like(state.albedo),
        0,
        0.0,
        copy(state.smb_ice),
        copy(state.runoff),
        copy(state.melt),
        copy(state.refreezing),
        copy(state.sublimation),
        copy(state.latent_heat_flux_sum),
    )
end

function MonthlyOutputBuffer(state::BESSIState; nmonth::Integer=12)
    return MonthlyOutputBuffer(
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        similar(state.runoff, Float32, Int(nmonth), state.ncol),
        0,
    )
end

function accumulate_monthly!(monthly::MonthlyState, state::BESSIState, dt_days::Real=1.0)
    monthly.albedo .+= state.albedo
    monthly.count += 1
    monthly.days += Float64(dt_days)
    return monthly
end

function finalize_monthly!(monthly::MonthlyState, state::BESSIState)
    monthly.smb_ice .= state.smb_ice .- monthly.prev_smb_ice
    monthly.runoff .= state.runoff .- monthly.prev_runoff
    monthly.melt .= state.melt .- monthly.prev_melt
    monthly.refreezing .= state.refreezing .- monthly.prev_refreezing
    monthly.sublimation .= state.sublimation .- monthly.prev_sublimation
    if monthly.days > 0.0
        monthly.latent_heat_flux .=
            (state.latent_heat_flux_sum .- monthly.prev_latent_heat_flux_sum) ./ monthly.days
    else
        fill!(monthly.latent_heat_flux, zero(eltype(monthly.latent_heat_flux)))
    end
    if monthly.count > 0
        monthly.albedo ./= monthly.count
    else
        monthly.albedo .= state.albedo
    end
    monthly.prev_smb_ice .= state.smb_ice
    monthly.prev_runoff .= state.runoff
    monthly.prev_melt .= state.melt
    monthly.prev_refreezing .= state.refreezing
    monthly.prev_sublimation .= state.sublimation
    monthly.prev_latent_heat_flux_sum .= state.latent_heat_flux_sum
    return monthly
end

function reset_monthly!(monthly::MonthlyState)
    for field in (:smb_ice, :runoff, :melt, :refreezing, :sublimation, :latent_heat_flux, :albedo)
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
    return Expr(:block, assignments..., :(nothing))
end

@kernel function _store_named_monthly_fields_kernel!(output, source, row::Int, ncol::Int)
    idx = @index(Global)
    if idx <= ncol
        _store_named_fields_at!(output, source, row, idx)
    end
end
function _store_named_monthly_fields!(output, source, names, row::Int, reference)
    output_fields = _named_fields(output, names)
    source_fields = _named_fields(source, names)
    kernel! = _store_named_monthly_fields_kernel!(_ka_backend(reference))
    event = kernel!(output_fields, source_fields, row, length(reference); ndrange=length(reference))
    _wait_kernel(event)
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
