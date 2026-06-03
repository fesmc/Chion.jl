"""State containers owned by `Simulation.ref` and `Simulation.now`."""

"""Abstract supertype for model state owned by `Simulation`."""
abstract type AbstractSnowModelState <: AbstractState end

"""
    CurrentState

Flat BESSI state container. Evolving arrays live directly on the state so user
code can inspect fields as `sim.now.mass`, `sim.now.runoff`, and
`sim.now.thickness` without going through an intermediate domain wrapper.
"""
mutable struct CurrentState{
        NF <: AbstractFloat,
        NI <: AbstractVector{<:Integer},
        MT <: AbstractMatrix{NF},
        VT <: AbstractVector{NF},
    } <: AbstractSnowpackDomain{NF}
    c::SnowpackPhysicalConstants{NF}
    Ntot::Int
    ncol::Int
    mass_max::NF
    mass_split::NF
    mass_min::NF
    rho_max::NF
    N::NI
    mass::MT
    mass_w::MT
    density::MT
    temperature::MT
    mass_base::VT
    smb_ice::VT
    runoff::VT
    melt::VT
    refreezing::VT
    vapor_mass::VT
    sublimation::VT
    latent_heat_flux_sum::VT
    Tsrf::VT
    albedo::VT
    thickness::VT
    wet_mass::VT
    bulk_density::VT
    liquid_water::VT
end

"""
    ReferenceState

Snapshot of a BESSI state at simulation initialization. Arrays are copied from
the initial `CurrentState` and are not advanced by the runtime.
"""
struct ReferenceState{
        NF <: AbstractFloat,
        NI <: AbstractVector{<:Integer},
        MT <: AbstractMatrix{NF},
        VT <: AbstractVector{NF},
    } <: AbstractState
    c::SnowpackPhysicalConstants{NF}
    Ntot::Int
    ncol::Int
    mass_max::NF
    mass_split::NF
    mass_min::NF
    rho_max::NF
    N::NI
    mass::MT
    mass_w::MT
    density::MT
    temperature::MT
    mass_base::VT
    smb_ice::VT
    runoff::VT
    melt::VT
    refreezing::VT
    vapor_mass::VT
    sublimation::VT
    latent_heat_flux_sum::VT
    Tsrf::VT
    albedo::VT
    thickness::VT
    wet_mass::VT
    bulk_density::VT
    liquid_water::VT
end

function CurrentState(domain::SnowpackDomain; density_init::Real=DEFAULT_DENSITY_INIT, temperature_init::Real=DEFAULT_TEMPERATURE_INIT)
    NF = number_type(domain.c)
    return CurrentState(
        domain.c,
        domain.Ntot,
        domain.ncol,
        domain.mass_max,
        domain.mass_split,
        domain.mass_min,
        domain.rho_max,
        zeros(Int, domain.ncol),
        zeros(NF, domain.Ntot, domain.ncol),
        zeros(NF, domain.Ntot, domain.ncol),
        fill(convert(NF, density_init), domain.Ntot, domain.ncol),
        fill(convert(NF, temperature_init), domain.Ntot, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        fill(domain.c.T0, domain.ncol),
        fill(domain.c.alpha_dry, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
        zeros(NF, domain.ncol),
    )
end

function CurrentState(model::BESSIModel)
    return CurrentState(SnowpackDomain(model); density_init=model.density_init, temperature_init=model.temperature_init)
end

ReferenceState(state::CurrentState) = ReferenceState(
    state.c,
    state.Ntot,
    state.ncol,
    state.mass_max,
    state.mass_split,
    state.mass_min,
    state.rho_max,
    copy(state.N),
    copy(state.mass),
    copy(state.mass_w),
    copy(state.density),
    copy(state.temperature),
    copy(state.mass_base),
    copy(state.smb_ice),
    copy(state.runoff),
    copy(state.melt),
    copy(state.refreezing),
    copy(state.vapor_mass),
    copy(state.sublimation),
    copy(state.latent_heat_flux_sum),
    copy(state.Tsrf),
    copy(state.albedo),
    copy(state.thickness),
    copy(state.wet_mass),
    copy(state.bulk_density),
    copy(state.liquid_water),
)

@inline Base.getproperty(state::Union{CurrentState,ReferenceState}, name::Symbol) =
    name === :albedo_dynamic ? getfield(state, :albedo) : getfield(state, name)

function Base.propertynames(state::Union{CurrentState,ReferenceState}, private::Bool=false)
    names = fieldnames(typeof(state))
    return private ? names : (names..., :albedo_dynamic)
end

cpu_state(state::CurrentState) = adapt(Array, state)
cpu_domain(state::CurrentState) = cpu_state(state)

function gpu_state(state::CurrentState, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return adapt(storage_type, state)
end
gpu_domain(state::CurrentState, storage_type=gpu_storage_type()) = gpu_state(state, storage_type)

@adapt_structure CurrentState

mutable struct MonthlyState{VT <: AbstractVector{<:AbstractFloat}, MT <: AbstractMatrix{<:AbstractFloat}} <: AbstractState
    smb_ice::VT
    runoff::VT
    melt::VT
    refreezing::VT
    sublimation::VT
    latent_heat_flux::VT
    albedo::VT
    packed::MT
    count::Int
    days::Float64
    prev_smb_ice::VT
    prev_runoff::VT
    prev_melt::VT
    prev_refreezing::VT
    prev_sublimation::VT
    prev_latent_heat_flux_sum::VT
end

mutable struct MonthlyYearState{MT <: AbstractMatrix{<:AbstractFloat}} <: AbstractState
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

function MonthlyState(state::CurrentState)
    return MonthlyState(
        _zero_like(state.smb_ice),
        _zero_like(state.runoff),
        _zero_like(state.melt),
        _zero_like(state.refreezing),
        _zero_like(state.sublimation),
        _zero_like(state.latent_heat_flux_sum),
        _zero_like(state.albedo),
        similar(state.runoff, Float32, 7, length(state.runoff)),
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

function MonthlyYearState(state::CurrentState; nmonth::Integer=12)
    return MonthlyYearState(
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

cpu_state(state::MonthlyState) = adapt(Array, state)

@adapt_structure MonthlyState

@inline _wait_monthly_event(event, array) = array isa Array ? _wait_kernel(event) : nothing

@kernel function _snapshot_monthly_kernel!(
    packed,
    prev_smb_ice,
    prev_runoff,
    prev_melt,
    prev_refreezing,
    prev_sublimation,
    prev_latent_heat_flux_sum,
    smb_ice,
    runoff,
    melt,
    refreezing,
    sublimation,
    latent_heat_flux_sum,
    albedo,
    days,
)
    idx = @index(Global)
    if idx <= length(runoff)
        smb_now = smb_ice[idx]
        runoff_now = runoff[idx]
        melt_now = melt[idx]
        refreezing_now = refreezing[idx]
        sublimation_now = sublimation[idx]
        latent_heat_flux_sum_now = latent_heat_flux_sum[idx]
        packed[1, idx] = smb_now - prev_smb_ice[idx]
        packed[2, idx] = runoff_now - prev_runoff[idx]
        packed[3, idx] = melt_now - prev_melt[idx]
        packed[4, idx] = refreezing_now - prev_refreezing[idx]
        packed[5, idx] = sublimation_now - prev_sublimation[idx]
        packed[6, idx] = days > 0 ? (latent_heat_flux_sum_now - prev_latent_heat_flux_sum[idx]) / days : zero(days)
        packed[7, idx] = albedo[idx]
        prev_smb_ice[idx] = smb_now
        prev_runoff[idx] = runoff_now
        prev_melt[idx] = melt_now
        prev_refreezing[idx] = refreezing_now
        prev_sublimation[idx] = sublimation_now
        prev_latent_heat_flux_sum[idx] = latent_heat_flux_sum_now
    end
end

function snapshot_monthly!(monthly::MonthlyState, state::CurrentState; days::Real=monthly.days)
    kernel! = _snapshot_monthly_kernel!(_ka_backend(monthly.packed))
    event = kernel!(
        monthly.packed,
        monthly.prev_smb_ice,
        monthly.prev_runoff,
        monthly.prev_melt,
        monthly.prev_refreezing,
        monthly.prev_sublimation,
        monthly.prev_latent_heat_flux_sum,
        state.smb_ice,
        state.runoff,
        state.melt,
        state.refreezing,
        state.sublimation,
        state.latent_heat_flux_sum,
        state.albedo,
        Float64(days);
        ndrange=length(state.runoff),
    )
    _wait_monthly_event(event, monthly.packed)
    return monthly
end

function accumulate_monthly!(monthly::MonthlyState, state::CurrentState, dt_days::Real=1.0)
    monthly.albedo .+= state.albedo
    monthly.count += 1
    monthly.days += Float64(dt_days)
    return monthly
end

function finalize_monthly!(monthly::MonthlyState, state::CurrentState)
    monthly.smb_ice .= state.smb_ice .- monthly.prev_smb_ice
    monthly.runoff .= state.runoff .- monthly.prev_runoff
    monthly.melt .= state.melt .- monthly.prev_melt
    monthly.refreezing .= state.refreezing .- monthly.prev_refreezing
    monthly.sublimation .= state.sublimation .- monthly.prev_sublimation
    if monthly.days > 0.0
        monthly.latent_heat_flux .= (state.latent_heat_flux_sum .- monthly.prev_latent_heat_flux_sum) ./ monthly.days
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
    fill!(monthly.smb_ice, zero(eltype(monthly.smb_ice)))
    fill!(monthly.runoff, zero(eltype(monthly.runoff)))
    fill!(monthly.melt, zero(eltype(monthly.melt)))
    fill!(monthly.refreezing, zero(eltype(monthly.refreezing)))
    fill!(monthly.sublimation, zero(eltype(monthly.sublimation)))
    fill!(monthly.latent_heat_flux, zero(eltype(monthly.latent_heat_flux)))
    fill!(monthly.albedo, zero(eltype(monthly.albedo)))
    monthly.count = 0
    monthly.days = 0.0
    return monthly
end

@kernel function _store_monthly_fields_kernel!(
    year_smb_ice,
    year_runoff,
    year_melt,
    year_refreezing,
    year_sublimation,
    year_latent_heat_flux,
    year_albedo,
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
        year_smb_ice[row, idx] = monthly_smb_ice[idx]
        year_runoff[row, idx] = monthly_runoff[idx]
        year_melt[row, idx] = monthly_melt[idx]
        year_refreezing[row, idx] = monthly_refreezing[idx]
        year_sublimation[row, idx] = monthly_sublimation[idx]
        year_latent_heat_flux[row, idx] = monthly_latent_heat_flux[idx]
        year_albedo[row, idx] = monthly_albedo[idx]
    end
end

function store_monthly!(year_state::MonthlyYearState, monthly::MonthlyState)
    row = year_state.count + 1
    row <= size(year_state.runoff, 1) || error("Monthly output year buffer is full.")
    kernel! = _store_monthly_fields_kernel!(_ka_backend(monthly.runoff))
    event = kernel!(
        year_state.smb_ice,
        year_state.runoff,
        year_state.melt,
        year_state.refreezing,
        year_state.sublimation,
        year_state.latent_heat_flux,
        year_state.albedo,
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
    year_state.count = row
    return year_state
end

function reset_monthly_year!(year_state::MonthlyYearState)
    year_state.count = 0
    return year_state
end

"""Mutable state for `PDDModel`."""
struct PDDState{ST <: AbstractVector{Float64}} <: AbstractSnowModelState
    snowpack_swe::ST
    smb_ice::ST
    runoff::ST
    pdd_sum::ST
end

function PDDState(model::PDDModel)
    ncol = ncols(model.grid)
    return PDDState(
        zeros(Float64, ncol),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
    )
end

initial_state(model::PDDModel) = PDDState(model)

initial_state(model::BESSIModel) = CurrentState(model)
reference_state(state::CurrentState) = ReferenceState(state)
reference_state(state::AbstractSnowModelState) = deepcopy(state)

function _copy_current_state!(dest::CurrentState, src::CurrentState)
    dest.Ntot == src.Ntot || error("Cannot copy state with different `Ntot`.")
    dest.ncol == src.ncol || error("Cannot copy state with different column count.")
    dest.N .= src.N
    dest.mass .= src.mass
    dest.mass_w .= src.mass_w
    dest.density .= src.density
    dest.temperature .= src.temperature
    dest.mass_base .= src.mass_base
    dest.smb_ice .= src.smb_ice
    dest.runoff .= src.runoff
    dest.melt .= src.melt
    dest.refreezing .= src.refreezing
    dest.vapor_mass .= src.vapor_mass
    dest.sublimation .= src.sublimation
    dest.latent_heat_flux_sum .= src.latent_heat_flux_sum
    dest.Tsrf .= src.Tsrf
    dest.albedo .= src.albedo
    dest.thickness .= src.thickness
    dest.wet_mass .= src.wet_mass
    dest.bulk_density .= src.bulk_density
    dest.liquid_water .= src.liquid_water
    return dest
end

function update_diagnostics!(state::CurrentState)
    summarize_domain_state!(
        state.thickness,
        state.wet_mass,
        state.bulk_density,
        state.mass_base,
        state.smb_ice,
        state.liquid_water,
        state.runoff,
        state.melt,
        state.refreezing,
        state.vapor_mass,
        state.sublimation,
        state.latent_heat_flux_sum,
        state.albedo,
        state,
    )
    return state
end
