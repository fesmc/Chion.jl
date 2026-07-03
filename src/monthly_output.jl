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

@inline _wait_monthly_event(event, array) = array isa Array ? _wait_kernel(event) : nothing

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

@kernel function _store_monthly_fields_kernel!(
    output_smb_ice,
    output_runoff,
    output_melt,
    output_refreezing,
    output_sublimation,
    output_latent_heat_flux,
    output_albedo,
    monthly_smb_ice,
    monthly_runoff,
    monthly_melt,
    monthly_refreezing,
    monthly_sublimation,
    monthly_latent_heat_flux,
    monthly_albedo,
    row::Int,
)
    idx = @index(Global)
    if idx <= length(monthly_runoff)
        output_smb_ice[row, idx] = monthly_smb_ice[idx]
        output_runoff[row, idx] = monthly_runoff[idx]
        output_melt[row, idx] = monthly_melt[idx]
        output_refreezing[row, idx] = monthly_refreezing[idx]
        output_sublimation[row, idx] = monthly_sublimation[idx]
        output_latent_heat_flux[row, idx] = monthly_latent_heat_flux[idx]
        output_albedo[row, idx] = monthly_albedo[idx]
    end
end

function store_monthly!(output::MonthlyOutputBuffer, monthly::MonthlyState)
    row = output.count + 1
    row <= size(output.runoff, 1) || error("Monthly output buffer is full.")
    kernel! = _store_monthly_fields_kernel!(_ka_backend(monthly.runoff))
    event = kernel!(
        output.smb_ice,
        output.runoff,
        output.melt,
        output.refreezing,
        output.sublimation,
        output.latent_heat_flux,
        output.albedo,
        monthly.smb_ice,
        monthly.runoff,
        monthly.melt,
        monthly.refreezing,
        monthly.sublimation,
        monthly.latent_heat_flux,
        monthly.albedo,
        row;
        ndrange=length(monthly.runoff),
    )
    _wait_monthly_event(event, monthly.runoff)
    output.count = row
    return output
end

function pdd_monthly_output_buffer(runtime; nmonth::Integer)
    dims = (Int(nmonth), length(runtime.snowpack_swe))
    allocate(field) = similar(field, Float32, dims)
    return (
        snowpack_swe=allocate(runtime.snowpack_swe),
        smb_ice=allocate(runtime.smb_ice),
        runoff=allocate(runtime.runoff),
        pdd_sum=allocate(runtime.pdd_sum),
    )
end

@kernel function _store_pdd_monthly_fields_kernel!(
    output_snowpack_swe,
    output_smb_ice,
    output_runoff,
    output_pdd_sum,
    snowpack_swe,
    smb_ice,
    runoff,
    pdd_sum,
    row::Int,
)
    idx = @index(Global)
    if idx <= length(snowpack_swe)
        output_snowpack_swe[row, idx] = snowpack_swe[idx]
        output_smb_ice[row, idx] = smb_ice[idx]
        output_runoff[row, idx] = runoff[idx]
        output_pdd_sum[row, idx] = pdd_sum[idx]
    end
end

function store_pdd_monthly!(output, monthly, row::Int)
    row <= size(output.runoff, 1) || error("PDD monthly output buffer is full.")
    kernel! = _store_pdd_monthly_fields_kernel!(_ka_backend(monthly.runoff))
    event = kernel!(
        output.snowpack_swe,
        output.smb_ice,
        output.runoff,
        output.pdd_sum,
        monthly.snowpack_swe,
        monthly.smb_ice,
        monthly.runoff,
        monthly.pdd_sum,
        row;
        ndrange=length(monthly.runoff),
    )
    _wait_monthly_event(event, monthly.runoff)
    return output
end
