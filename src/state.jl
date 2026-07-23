"""State containers owned by `Simulation.ref` and `Simulation.now`."""

"""
Logical `(layer, column)` matrix backed by `(column, layer)` storage.

The transposed physical layout gives consecutive GPU threads coalesced access
while preserving the existing layer-first indexing used by the physics code.
"""
struct TransposedLayerMatrix{T, M <: AbstractMatrix{T}} <: AbstractMatrix{T}
    parent::M
end

Base.size(matrix::TransposedLayerMatrix) = reverse(size(matrix.parent))
@inline function Base.getindex(matrix::TransposedLayerMatrix, layer::Int, column::Int)
    @boundscheck checkbounds(matrix, layer, column)
    return @inbounds matrix.parent[column, layer]
end
@inline function Base.setindex!(matrix::TransposedLayerMatrix, value, layer::Int, column::Int)
    @boundscheck checkbounds(matrix, layer, column)
    @inbounds matrix.parent[column, layer] = value
    return value
end
Base.IndexStyle(::Type{<:TransposedLayerMatrix}) = IndexCartesian()
KernelAbstractions.get_backend(matrix::TransposedLayerMatrix) =
    KernelAbstractions.get_backend(matrix.parent)
Base.similar(matrix::TransposedLayerMatrix, ::Type{T}, nlayer::Int, ncol::Int) where {T} =
    TransposedLayerMatrix(similar(matrix.parent, T, ncol, nlayer))
Base.similar(matrix::TransposedLayerMatrix, ::Type{T}, n::Int) where {T} =
    similar(matrix.parent, T, n)
Base.Array(matrix::TransposedLayerMatrix) = permutedims(Array(matrix.parent), (2, 1))
Adapt.@adapt_structure TransposedLayerMatrix

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

@inline _cpu_layer_matrix(matrix::TransposedLayerMatrix) = Array(matrix)
@inline _cpu_layer_matrix(matrix) = Array(matrix)

function cpu_state(state::BESSIState)
    return BESSIState(
        state.c,
        state.Ntot,
        state.ncol,
        state.mass_max,
        state.mass_split,
        state.mass_min,
        Array(state.N),
        _cpu_layer_matrix(state.mass),
        _cpu_layer_matrix(state.mass_w),
        _cpu_layer_matrix(state.density),
        _cpu_layer_matrix(state.temperature),
        Array(state.mass_base),
        Array(state.smb_ice),
        Array(state.runoff),
        Array(state.melt),
        Array(state.refreezing),
        Array(state.vapor_mass),
        Array(state.sublimation),
        Array(state.latent_heat_flux_sum),
        Array(state.Tsrf),
        Array(state.albedo),
        Array(state.thickness),
        Array(state.wet_mass),
        Array(state.bulk_density),
        Array(state.liquid_water),
    )
end

@inline _gpu_layer_matrix(matrix, storage_type) =
    TransposedLayerMatrix(adapt(storage_type, permutedims(Array(matrix), (2, 1))))

function gpu_state(state::BESSIState, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return BESSIState(
        state.c,
        state.Ntot,
        state.ncol,
        state.mass_max,
        state.mass_split,
        state.mass_min,
        adapt(storage_type, state.N),
        _gpu_layer_matrix(state.mass, storage_type),
        _gpu_layer_matrix(state.mass_w, storage_type),
        _gpu_layer_matrix(state.density, storage_type),
        _gpu_layer_matrix(state.temperature, storage_type),
        adapt(storage_type, state.mass_base),
        adapt(storage_type, state.smb_ice),
        adapt(storage_type, state.runoff),
        adapt(storage_type, state.melt),
        adapt(storage_type, state.refreezing),
        adapt(storage_type, state.vapor_mass),
        adapt(storage_type, state.sublimation),
        adapt(storage_type, state.latent_heat_flux_sum),
        adapt(storage_type, state.Tsrf),
        adapt(storage_type, state.albedo),
        adapt(storage_type, state.thickness),
        adapt(storage_type, state.wet_mass),
        adapt(storage_type, state.bulk_density),
        adapt(storage_type, state.liquid_water),
    )
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

@inline _copy_state_array!(dest, src) = copyto!(dest, src)
function _copy_state_array!(dest::AbstractMatrix, src::TransposedLayerMatrix)
    permutedims!(dest, Array(src.parent), (2, 1))
    return dest
end

function _copy_bessi_state!(dest::BESSIState, src::BESSIState)
    dest.Ntot == src.Ntot || error("Cannot copy state with different `Ntot`.")
    dest.ncol == src.ncol || error("Cannot copy state with different column count.")
    _copy_state_array!(dest.N, src.N)
    _copy_state_array!(dest.mass, src.mass)
    _copy_state_array!(dest.mass_w, src.mass_w)
    _copy_state_array!(dest.density, src.density)
    _copy_state_array!(dest.temperature, src.temperature)
    _copy_state_array!(dest.mass_base, src.mass_base)
    _copy_state_array!(dest.smb_ice, src.smb_ice)
    _copy_state_array!(dest.runoff, src.runoff)
    _copy_state_array!(dest.melt, src.melt)
    _copy_state_array!(dest.refreezing, src.refreezing)
    _copy_state_array!(dest.vapor_mass, src.vapor_mass)
    _copy_state_array!(dest.sublimation, src.sublimation)
    _copy_state_array!(dest.latent_heat_flux_sum, src.latent_heat_flux_sum)
    _copy_state_array!(dest.Tsrf, src.Tsrf)
    _copy_state_array!(dest.albedo, src.albedo)
    _copy_state_array!(dest.thickness, src.thickness)
    _copy_state_array!(dest.wet_mass, src.wet_mass)
    _copy_state_array!(dest.bulk_density, src.bulk_density)
    _copy_state_array!(dest.liquid_water, src.liquid_water)
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
