#=
Prognostic state of the models: the arrays that a time step reads and updates. Diagnostics and
scratch arrays are in the model's cache (`model.jl`).
=#

"""
Logical `(layer, column)` matrix backed by `(column, layer)` storage.

The transposed physical layout gives consecutive GPU threads coalesced access
while preserving the existing layer-first indexing used by the physics code.
"""
struct TransposedLayerMatrix{T, M <: AbstractMatrix{T}} <: AbstractMatrix{T}
    parent::M
end

Base.size(matrix::TransposedLayerMatrix) = reverse(size(matrix.parent))
Base.@propagate_inbounds function Base.getindex(matrix::TransposedLayerMatrix, layer::Int, column::Int)
    @boundscheck checkbounds(matrix, layer, column)
    return @inbounds matrix.parent[column, layer]
end
Base.@propagate_inbounds function Base.setindex!(matrix::TransposedLayerMatrix, value, layer::Int, column::Int)
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
$(TYPEDEF)

Prognostic state of a [`BESSIModel`](@ref): the snow and firn layers `(layer, column)` (on a GPU
stored `(column, layer)` for coalesced access, see `TransposedLayerMatrix`), the ice-substrate
temperatures, and per column the number of active layers, the cumulative mass fluxes (kg m⁻²),
the time-integrated surface energy balance (W m⁻² day, positive into the surface) and the surface
properties.

# Fields
$(TYPEDFIELDS)
"""
struct BESSIState{NI <: AbstractVector{<:Integer}, MT <: AbstractMatrix, VT <: AbstractVector}
    "number of active snow and firn layers"
    N::NI
    "solid mass of the layers (kg m⁻²)"
    mass::MT
    "liquid water of the layers (kg m⁻²)"
    mass_w::MT
    "density of the layers (kg m⁻³)"
    density::MT
    "temperature of the layers (K)"
    temperature::MT
    "temperature of the ice-substrate layers (K)"
    ice_temperature::MT
    "firn mass passed to the ice below (kg m⁻², cumulative)"
    mass_base::VT
    "net mass forcing of the ice sheet (kg m⁻², cumulative)"
    smb_ice::VT
    "runoff (kg m⁻², cumulative)"
    runoff::VT
    "melt (kg m⁻², cumulative)"
    melt::VT
    "refreezing (kg m⁻², cumulative)"
    refreezing::VT
    "net vapour exchange, deposition positive (kg m⁻², cumulative)"
    vapor_mass::VT
    "sublimation (kg m⁻², cumulative)"
    sublimation::VT
    "latent heat flux (W m⁻² day, cumulative)"
    latent_heat_flux_sum::VT
    "shortwave down (W m⁻² day, cumulative)"
    shortwave_down_sum::VT
    "absorbed shortwave (W m⁻² day, cumulative)"
    shortwave_net_sum::VT
    "longwave down (W m⁻² day, cumulative)"
    longwave_down_sum::VT
    "net longwave (W m⁻² day, cumulative)"
    longwave_net_sum::VT
    "sensible heat flux (W m⁻² day, cumulative)"
    sensible_heat_flux_sum::VT
    "heat of rain (W m⁻² day, cumulative)"
    rain_heat_flux_sum::VT
    "surface temperature (K)"
    Tsrf::VT
    "surface albedo (–)"
    albedo::VT
    "time since the latest snowfall, `AgingAlbedo` (day)"
    snow_age_days::VT
    "seasonal maximum of the snow water equivalent, `SemixAlbedo` (kg m⁻²)"
    w_snow_max::VT
end

# Layer-first matrices, transposed on a device.
const _BESSI_LAYER_FIELD_NAMES = (:mass, :mass_w, :density, :temperature, :ice_temperature)
# Surface energy balance accumulators, in the order of the state fields.
const _BESSI_SEB_SUM_NAMES = (
    :shortwave_down_sum,
    :shortwave_net_sum,
    :longwave_down_sum,
    :longwave_net_sum,
    :sensible_heat_flux_sum,
    :rain_heat_flux_sum,
)
const _BESSI_ARRAY_FIELD_NAMES = fieldnames(BESSIState)

"""The arrays of a state as a NamedTuple, the form the kernels receive."""
@inline get_fields(state) = _fields(state)

ncols(state::BESSIState) = length(state.N)
number_type(state::BESSIState) = eltype(state.mass)

@inline function _zero_seb_sums!(fields, idx)
    @inbounds begin
        fields.shortwave_down_sum[idx] = zero(eltype(fields.shortwave_down_sum))
        fields.shortwave_net_sum[idx] = zero(eltype(fields.shortwave_net_sum))
        fields.longwave_down_sum[idx] = zero(eltype(fields.longwave_down_sum))
        fields.longwave_net_sum[idx] = zero(eltype(fields.longwave_net_sum))
        fields.sensible_heat_flux_sum[idx] = zero(eltype(fields.sensible_heat_flux_sum))
        fields.rain_heat_flux_sum[idx] = zero(eltype(fields.rain_heat_flux_sum))
    end
    return nothing
end

"""Reset the columns `indices` (all with `1:ncol`) to the initial state: no snow layers, the
initial density and temperature in the empty layers and the substrate, zero counters."""
@kernel function _initialize_bessi_columns_kernel!(
    fields,
    indices,
    density_init,
    temperature_init,
    surface_temperature_init,
    albedo_init,
)
    item = @index(Global)
    @inbounds begin
        idx = indices[item]
        fields.N[idx] = 0
        fields.mass_base[idx] = zero(density_init)
        fields.smb_ice[idx] = zero(density_init)
        fields.runoff[idx] = zero(density_init)
        fields.melt[idx] = zero(density_init)
        fields.refreezing[idx] = zero(density_init)
        fields.vapor_mass[idx] = zero(density_init)
        fields.sublimation[idx] = zero(density_init)
        fields.latent_heat_flux_sum[idx] = zero(density_init)
        _zero_seb_sums!(fields, idx)
        fields.Tsrf[idx] = surface_temperature_init
        fields.albedo[idx] = albedo_init
        fields.snow_age_days[idx] = zero(density_init)
        fields.w_snow_max[idx] = zero(density_init)
        for layer_index in 1:size(fields.mass, 1)
            fields.mass[layer_index, idx] = zero(density_init)
            fields.mass_w[layer_index, idx] = zero(density_init)
            fields.density[layer_index, idx] = density_init
            fields.temperature[layer_index, idx] = temperature_init
        end
        for layer_index in 1:size(fields.ice_temperature, 1)
            fields.ice_temperature[layer_index, idx] = temperature_init
        end
    end
end

@inline _state_init_workgroupsize(::KernelAbstractions.CPU) = 1024
@inline _state_init_workgroupsize(::KernelAbstractions.Backend) = 256

function _initialize_columns!(state::BESSIState, params::Params, indices)
    launch!(_initialize_bessi_columns_kernel!, state.mass, length(indices),
        get_fields(state),
        indices,
        params.density_init,
        params.temperature_init,
        params.T0,
        params.alpha_dry;
        workgroupsize = _state_init_workgroupsize(_ka_backend(state.mass)),
    )
    return state
end

"""
    initial_state(grid, params) -> BESSIState | PDDState | ITMState

Initial state of the columns of `grid` on the CPU: no snow, the initial density and temperature
of `params` (a [`Params`](@ref)), or the initial bulk state for `PDDParams` and `ITMParams`.
"""
function initial_state(grid::Grid, params::Params)
    NF = number_type(params)
    ncol = ncols(grid)
    nlayer = grid.layering.Ntot
    layer() = Matrix{NF}(undef, nlayer, ncol)
    column() = Vector{NF}(undef, ncol)
    state = BESSIState(
        Vector{Int}(undef, ncol),
        layer(), layer(), layer(), layer(),
        Matrix{NF}(undef, grid.layering.ice_substrate_layers, ncol),
        ntuple(_ -> column(), fieldcount(BESSIState) - 6)...,
    )
    return _initialize_columns!(state, params, 1:ncol)
end

# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

"""A copy of the array `a` on the KernelAbstractions `backend`."""
_to_device(backend, a::AbstractArray) = (b = KernelAbstractions.allocate(backend, eltype(a), size(a)); copyto!(b, a); b)

# A layer matrix: `(layer, column)` on the CPU, transposed on a device.
_layer_to_backend(matrix, ::KernelAbstractions.CPU) = Array(matrix)
_layer_to_backend(matrix, backend) =
    TransposedLayerMatrix(_to_device(backend, permutedims(Array(matrix), (2, 1))))
_vector_to_backend(vector, ::KernelAbstractions.CPU) = Array(vector)
_vector_to_backend(vector, backend) = _to_device(backend, Array(vector))

"""
    to_backend(state, backend)

A copy of `state` with every array on the KernelAbstractions `backend` (`CPU()`, or e.g.
`CUDABackend()` after `using CUDA`).
"""
function to_backend(state::BESSIState, backend)
    arrays = map(fieldnames(BESSIState)) do name
        value = getfield(state, name)
        name in _BESSI_LAYER_FIELD_NAMES ? _layer_to_backend(value, backend) : _vector_to_backend(value, backend)
    end
    return BESSIState(arrays...)
end

"""
$(TYPEDEF)

State of a [`PDDModel`](@ref): the snow reservoir and the cumulative mass fluxes (mm w.e.).

# Fields
$(TYPEDFIELDS)
"""
struct PDDState{ST <: AbstractVector{<:Real}}
    "snow water equivalent (mm w.e.)"
    snowpack_swe::ST
    "net mass forcing of the ice sheet (mm w.e., cumulative)"
    smb_ice::ST
    "runoff (mm w.e., cumulative)"
    runoff::ST
    "positive degree days (K day, cumulative)"
    pdd_sum::ST
end

const _PDD_STATE_FIELD_NAMES = fieldnames(PDDState)

function initial_state(grid::Grid, params::PDDParams)
    NF = number_type(params)
    ncol = ncols(grid)
    return PDDState(zeros(NF, ncol), zeros(NF, ncol), zeros(NF, ncol), zeros(NF, ncol))
end

"""
$(TYPEDEF)

State of an [`ITMModel`](@ref). The instantaneous budget fields are in ITM's native
mm w.e. day⁻¹; the cumulative fields (`*_cum`, `smb_ice`) in mm w.e.

# Fields
$(TYPEDFIELDS)
"""
struct ITMState{ST <: AbstractVector{<:Real}}
    "snow water equivalent (mm w.e.)"
    H_snow::ST
    "surface albedo (–)"
    alb_s::ST
    "surface mass balance (mm w.e. day⁻¹)"
    smb::ST
    "ice-sheet mass balance (mm w.e. day⁻¹)"
    smbi::ST
    "melt (mm w.e. day⁻¹)"
    melt::ST
    "runoff (mm w.e. day⁻¹)"
    runoff::ST
    "refreezing (mm w.e. day⁻¹)"
    refreezing::ST
    "surface temperature (K)"
    Tsrf::ST
    "net melt (mm w.e. day⁻¹)"
    melt_net::ST
    "surface mass balance (mm w.e., cumulative)"
    smb_cum::ST
    "ice-sheet mass balance (mm w.e., cumulative)"
    smb_ice::ST
    "melt (mm w.e., cumulative)"
    melt_cum::ST
    "runoff (mm w.e., cumulative)"
    runoff_cum::ST
    "refreezing (mm w.e., cumulative)"
    refreezing_cum::ST
end

const _ITM_STATE_FIELD_NAMES = fieldnames(ITMState)

function initial_state(grid::Grid, params::ITMParams)
    NF = number_type(params)
    ncol = ncols(grid)
    return ITMState(
        fill(params.H_snow_max, ncol), fill(params.alb_snow_dry, ncol), zeros(NF, ncol), zeros(NF, ncol),
        zeros(NF, ncol), zeros(NF, ncol), zeros(NF, ncol), fill(params.T0, ncol), zeros(NF, ncol),
        zeros(NF, ncol), zeros(NF, ncol), zeros(NF, ncol), zeros(NF, ncol), zeros(NF, ncol),
    )
end

ncols(state::Union{PDDState, ITMState}) = length(state.runoff)
number_type(state::Union{PDDState, ITMState}) = eltype(state.runoff)
to_backend(state::Union{PDDState, ITMState}, backend) = _mapfields(v -> _vector_to_backend(v, backend), state)

Adapt.@adapt_structure BESSIState
Adapt.@adapt_structure PDDState
Adapt.@adapt_structure ITMState
