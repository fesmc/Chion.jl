"""State containers owned by `Simulation.ref` and `Simulation.now`."""

"""
    BESSIState

Flat BESSI state container. Evolving arrays live directly on the state so user
code can inspect fields as `sim.now.mass`, `sim.now.runoff`, and
`sim.now.thickness` without going through an intermediate domain wrapper.
"""
struct BESSIState{
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

@kernel function _initialize_bessi_state_kernel!(
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

function _initialize_bessi_state_arrays!(
    state::BESSIState,
    density_init,
    temperature_init,
)
    backend = _ka_backend(state.mass)
    kernel! = _initialize_bessi_state_kernel!(backend, _state_init_workgroupsize(backend))
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

function BESSIState(model::BESSIModel)
    NF = number_type(model.c)
    ncol = ncols(model.grid)
    state = BESSIState(
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
    return _initialize_bessi_state_arrays!(state, model.density_init, model.temperature_init)
end

cpu_state(state::BESSIState) = adapt(Array, state)

function gpu_state(state::BESSIState, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return adapt(storage_type, state)
end

@adapt_structure BESSIState
ncols(state::BESSIState) = state.ncol

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

initial_state(model::BESSIModel) = BESSIState(model)
reference_state(state::BESSIState) = deepcopy(state)
reference_state(state::PDDState) = deepcopy(state)

function _copy_bessi_state!(dest::BESSIState, src::BESSIState)
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

function update_diagnostics!(state::BESSIState)
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
