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

const _BESSI_LAYER_FIELD_NAMES = (:mass, :mass_w, :density, :temperature)
const _BESSI_ARRAY_FIELD_NAMES = (
    :N,
    _BESSI_LAYER_FIELD_NAMES...,
    :mass_base,
    :smb_ice,
    :runoff,
    :melt,
    :refreezing,
    :vapor_mass,
    :sublimation,
    :latent_heat_flux_sum,
    :Tsrf,
    :albedo,
    :thickness,
    :wet_mass,
    :bulk_density,
    :liquid_water,
)

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
@inline _cpu_state_array(state::BESSIState, name::Symbol) =
    name in _BESSI_LAYER_FIELD_NAMES ?
    _cpu_layer_matrix(getfield(state, name)) :
    Array(getfield(state, name))

function cpu_state(state::BESSIState)
    return BESSIState(
        state.c,
        state.Ntot,
        state.ncol,
        state.mass_max,
        state.mass_split,
        state.mass_min,
        (_cpu_state_array(state, name) for name in _BESSI_ARRAY_FIELD_NAMES)...,
    )
end

@inline _gpu_layer_matrix(matrix, storage_type) =
    TransposedLayerMatrix(adapt(storage_type, permutedims(Array(matrix), (2, 1))))
@inline _gpu_state_array(state::BESSIState, name::Symbol, storage_type) =
    name in _BESSI_LAYER_FIELD_NAMES ?
    _gpu_layer_matrix(getfield(state, name), storage_type) :
    adapt(storage_type, getfield(state, name))

function gpu_state(state::BESSIState, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return BESSIState(
        state.c,
        state.Ntot,
        state.ncol,
        state.mass_max,
        state.mass_split,
        state.mass_min,
        (_gpu_state_array(state, name, storage_type) for name in _BESSI_ARRAY_FIELD_NAMES)...,
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

"""Mutable state for the bulk Fortran-compatible `ITMModel`.

Instantaneous budget fields use ITM's native mm w.e. day⁻¹ units; cumulative
fields use mm w.e. and are the quantities exposed in regular output.
"""
struct ITMState{ST <: AbstractVector{Float64}}
    H_snow::ST
    alb_s::ST
    smb::ST
    smbi::ST
    melt::ST
    runoff::ST
    refreezing::ST
    Tsrf::ST
    melt_net::ST
    smb_cum::ST
    smb_ice::ST
    melt_cum::ST
    runoff_cum::ST
    refreezing_cum::ST
end

function ITMState(model::ITMModel)
    ncol = ncols(model.grid)
    return ITMState(
        fill(model.H_snow_max, ncol), fill(model.alb_snow_dry, ncol), zeros(ncol), zeros(ncol),
        zeros(ncol), zeros(ncol), zeros(ncol), fill(model.c.T0, ncol), zeros(ncol),
        zeros(ncol), zeros(ncol), zeros(ncol), zeros(ncol), zeros(ncol),
    )
end

initial_state(model::ITMModel) = ITMState(model)

initial_state(model::BESSIModel) = BESSIState(model)
reference_state(state::BESSIState) = deepcopy(state)
reference_state(state::PDDState) = deepcopy(state)
reference_state(state::ITMState) = deepcopy(state)

@inline _copy_state_array!(dest, src) = copyto!(dest, src)
function _copy_state_array!(dest::AbstractMatrix, src::TransposedLayerMatrix)
    permutedims!(dest, Array(src.parent), (2, 1))
    return dest
end

function _copy_bessi_state!(dest::BESSIState, src::BESSIState)
    dest.Ntot == src.Ntot || error("Cannot copy state with different `Ntot`.")
    dest.ncol == src.ncol || error("Cannot copy state with different column count.")
    for name in _BESSI_ARRAY_FIELD_NAMES
        _copy_state_array!(getfield(dest, name), getfield(src, name))
    end
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
