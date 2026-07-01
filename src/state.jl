"""State containers owned by `Simulation.ref` and `Simulation.now`."""

"""
    CurrentState

Flat BESSI state container. Evolving arrays live directly on the state so user
code can inspect fields as `sim.now.mass`, `sim.now.runoff`, and
`sim.now.thickness` without going through an intermediate domain wrapper.
"""
struct CurrentState{
        NF <: AbstractFloat,
        NI <: AbstractVector{<:Integer},
        MT <: AbstractMatrix{NF},
        VT <: AbstractVector{NF},
    }
    c::SnowpackPhysicalConstants{NF}
    Ntot::Int
    ncol::Int
    mass_max::NF
    mass_split::NF
    mass_min::NF
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

@kernel function _initialize_current_state_kernel!(
    N,
    mass,
    mass_w,
    density,
    temperature,
    mass_base,
    smb_ice,
    runoff,
    melt,
    refreezing,
    vapor_mass,
    sublimation,
    latent_heat_flux_sum,
    Tsrf,
    albedo,
    thickness,
    wet_mass,
    bulk_density,
    liquid_water,
    Ntot::Int,
    density_init,
    temperature_init,
    surface_temperature_init,
    albedo_init,
)
    idx = @index(Global)
    if idx <= length(N)
        N[idx] = 0
        mass_base[idx] = zero(density_init)
        smb_ice[idx] = zero(density_init)
        runoff[idx] = zero(density_init)
        melt[idx] = zero(density_init)
        refreezing[idx] = zero(density_init)
        vapor_mass[idx] = zero(density_init)
        sublimation[idx] = zero(density_init)
        latent_heat_flux_sum[idx] = zero(density_init)
        Tsrf[idx] = surface_temperature_init
        albedo[idx] = albedo_init
        thickness[idx] = zero(density_init)
        wet_mass[idx] = zero(density_init)
        bulk_density[idx] = zero(density_init)
        liquid_water[idx] = zero(density_init)
        for layer_index in 1:Ntot
            mass[layer_index, idx] = zero(density_init)
            mass_w[layer_index, idx] = zero(density_init)
            density[layer_index, idx] = density_init
            temperature[layer_index, idx] = temperature_init
        end
    end
end

@inline _state_init_workgroupsize(backend) =
    backend isa KernelAbstractions.CPU ? 1024 : 256

function _initialize_current_state_arrays!(
    state::CurrentState,
    density_init,
    temperature_init,
)
    backend = _ka_backend(state.mass)
    kernel! = _initialize_current_state_kernel!(backend, _state_init_workgroupsize(backend))
    event = kernel!(
        state.N,
        state.mass,
        state.mass_w,
        state.density,
        state.temperature,
        state.mass_base,
        state.smb_ice,
        state.runoff,
        state.melt,
        state.refreezing,
        state.vapor_mass,
        state.sublimation,
        state.latent_heat_flux_sum,
        state.Tsrf,
        state.albedo,
        state.thickness,
        state.wet_mass,
        state.bulk_density,
        state.liquid_water,
        state.Ntot,
        convert(eltype(state.mass), density_init),
        convert(eltype(state.mass), temperature_init),
        state.c.T0,
        state.c.alpha_dry;
        ndrange=state.ncol,
    )
    _wait_kernel(event)
    return state
end

function CurrentState(model::BESSIModel)
    NF = number_type(model.c)
    ncol = ncols(model.grid)
    state = CurrentState(
        model.c,
        model.Ntot,
        ncol,
        model.mass_max,
        model.mass_split,
        model.mass_min,
        Vector{Int}(undef, ncol),
        Matrix{NF}(undef, model.Ntot, ncol),
        Matrix{NF}(undef, model.Ntot, ncol),
        Matrix{NF}(undef, model.Ntot, ncol),
        Matrix{NF}(undef, model.Ntot, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
    )
    return _initialize_current_state_arrays!(state, model.density_init, model.temperature_init)
end

cpu_state(state::CurrentState) = adapt(Array, state)

function gpu_state(state::CurrentState, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return adapt(storage_type, state)
end

@adapt_structure CurrentState

mutable struct MonthlyState{VT <: AbstractVector{<:AbstractFloat}, MT <: AbstractMatrix{<:AbstractFloat}}
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

function MonthlyOutputBuffer(state::CurrentState; nmonth::Integer=12)
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

function reset_monthly_output!(output::MonthlyOutputBuffer)
    output.count = 0
    return output
end

"""Mutable state for `PDDModel`."""
struct PDDState{ST <: AbstractVector{Float64}}
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
reference_state(state::CurrentState) = deepcopy(state)
reference_state(state::PDDState) = deepcopy(state)

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
