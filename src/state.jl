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
    BESSIState

Flat BESSI state container. Evolving arrays live directly on the state so user
code can inspect fields as `sim.now.mass`, `sim.now.runoff`, and
`sim.now.thickness` without going through an intermediate domain wrapper.
"""
struct BESSIState{
        NF <: AbstractFloat,
        P <: BESSIParameters,
        NI <: AbstractVector{<:Integer},
        MT <: AbstractMatrix{NF},
        VT <: AbstractVector{NF},
    }
    parameters::P
    ncol::Int
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
    snow_age_days::VT
    w_snow_max::VT
    thickness::VT
    wet_mass::VT
    bulk_density::VT
    liquid_water::VT
    # Temperatures of the optional fixed-geometry ice substrate below the
    # snow/firn column, `(ice_substrate_layers, column)`.
    ice_temperature::MT
end

const _BESSI_STATE_PARAMETER_NAMES = (
    :c, :Ntot, :mass_max, :mass_split, :mass_min, :refreezing_correction,
)

@inline function Base.getproperty(state::BESSIState, name::Symbol)
    hasfield(typeof(state), name) && return getfield(state, name)
    name in _BESSI_STATE_PARAMETER_NAMES &&
        return getproperty(getfield(state, :parameters), name)
    return getfield(state, name)
end

function Base.propertynames(state::BESSIState, private::Bool=false)
    public = (_BESSI_STATE_PARAMETER_NAMES..., :ncol, _BESSI_ARRAY_FIELD_NAMES...)
    return private ? (:parameters, public...) : public
end

const _BESSI_LAYER_FIELD_NAMES = (:mass, :mass_w, :density, :temperature)
# Layer-first matrices that share the layer storage conversions.
@inline _is_bessi_layer_field(name::Symbol) = name in _BESSI_LAYER_FIELD_NAMES || name === :ice_temperature
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
    :snow_age_days,
    :w_snow_max,
    :thickness,
    :wet_mass,
    :bulk_density,
    :liquid_water,
    :ice_temperature,
)

@inline get_fields(state::BESSIState) = NamedTuple{_BESSI_ARRAY_FIELD_NAMES}(
    ntuple(index -> getfield(state, _BESSI_ARRAY_FIELD_NAMES[index]), Val(length(_BESSI_ARRAY_FIELD_NAMES))),
)

@kernel function _initialize_bessi_state_kernel!(
    fields,
    Ntot::Int,
    density_init,
    temperature_init,
    surface_temperature_init,
    albedo_init,
)
    idx = @index(Global)
    @inbounds begin
        fields.N[idx] = 0
        fields.mass_base[idx] = zero(density_init)
        fields.smb_ice[idx] = zero(density_init)
        fields.runoff[idx] = zero(density_init)
        fields.melt[idx] = zero(density_init)
        fields.refreezing[idx] = zero(density_init)
        fields.vapor_mass[idx] = zero(density_init)
        fields.sublimation[idx] = zero(density_init)
        fields.latent_heat_flux_sum[idx] = zero(density_init)
        fields.Tsrf[idx] = surface_temperature_init
        fields.albedo[idx] = albedo_init
        fields.snow_age_days[idx] = zero(density_init)
        fields.w_snow_max[idx] = zero(density_init)
        fields.thickness[idx] = zero(density_init)
        fields.wet_mass[idx] = zero(density_init)
        fields.bulk_density[idx] = zero(density_init)
        fields.liquid_water[idx] = zero(density_init)
        for layer_index in 1:Ntot
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
        get_fields(state),
        state.Ntot,
        convert(eltype(state.mass), density_init),
        convert(eltype(state.mass), temperature_init),
        state.c.T0,
        _initial_snow_albedo(state.c);
        ndrange=state.ncol,
    )
    _wait_kernel(event)
    return state
end

@kernel function _initialize_bessi_snowpack_kernel!(
    fields,
    swe,
    T_init,
    parameters::BESSIParameters,
    density_firn,
    albedo_snow,
)
    idx = @index(Global)
    @inbounds begin
        c = parameters.c
        temperature = min(T_init[idx], c.T0)
        fields.Tsrf[idx] = temperature
        for layer_index in 1:size(fields.ice_temperature, 1)
            fields.ice_temperature[layer_index, idx] = temperature
        end
        if swe[idx] > zero(swe[idx])
            _initial_firn_layers!(
                fields.N, fields.mass, fields.mass_w, fields.density, fields.temperature,
                fields.mass_base, fields.smb_ice, idx, swe[idx], density_firn, temperature,
                parameters.Ntot, oftype(density_firn, parameters.mass_split),
                oftype(density_firn, parameters.mass_min),
                parameters.near_surface_layer_max_thicknesses_m, c,
            )
            fields.albedo[idx] = albedo_snow
            fields.snow_age_days[idx] = zero(density_firn)
            thickness, wet_mass, bulk_density, liquid_water =
                _column_summary(fields.N, fields.mass, fields.mass_w, fields.density, idx, fields.thickness)
            fields.thickness[idx] = thickness
            fields.wet_mass[idx] = wet_mass
            fields.bulk_density[idx] = bulk_density
            fields.liquid_water[idx] = liquid_water
        end
    end
end

_initial_column_values(::Type{NF}, value::Real, ncol::Int, ::Symbol) where {NF} = fill(NF(value), ncol)
function _initial_column_values(::Type{NF}, values::AbstractVector{<:Real}, ncol::Int, name::Symbol) where {NF}
    length(values) == ncol || error("`$name` must have one value per column ($ncol), got $(length(values)).")
    return Vector{NF}(values)
end

function _initialize_bessi_snowpack!(state::BESSIState, swe, T_init)
    NF = eltype(state.mass)
    ncol = state.ncol
    swe_columns = isnothing(swe) ? zeros(NF, ncol) : _initial_column_values(NF, swe, ncol, :swe)
    T_columns = _initial_column_values(NF, T_init, ncol, :T_init)
    all(value -> isfinite(value) && value >= 0, swe_columns) ||
        error("`swe` must be finite and non-negative.")
    all(value -> isfinite(value) && value > 0, T_columns) ||
        error("`T_init` must be finite and positive [K].")
    density_firn = state.parameters.density_init_firn
    max_swe = density_firn * BESSI_REFERENCE_SNOW_DEPTH_M
    all(<=(max_swe), swe_columns) || error(
        "`swe` exceeds the snow-depth cap: at most $(max_swe) kg m-2 at " *
        "`density_init_firn` = $(density_firn) kg m-3 ($(BESSI_REFERENCE_SNOW_DEPTH_M) m).",
    )
    backend = _ka_backend(state.mass)
    kernel! = _initialize_bessi_snowpack_kernel!(backend, _state_init_workgroupsize(backend))
    event = kernel!(
        get_fields(state),
        copyto!(similar(state.Tsrf, NF, ncol), swe_columns),
        copyto!(similar(state.Tsrf, NF, ncol), T_columns),
        state.parameters,
        convert(NF, density_firn),
        _initial_snow_albedo(state.c);
        ndrange=ncol,
    )
    _wait_kernel(event)
    return state
end

"""
    BESSIState(model; swe=nothing, T_init=nothing)

Initial BESSI state for `model`. Without keywords every column starts
snow-free, with its ice substrate at `model.temperature_init` and `Tsrf = T0`.

`swe` [kg m⁻²] and `T_init` [K], each a number or one value per column, start
the columns with a dry firn column (a cold start for ice-sheet hosts). A column
with `swe > 0` holds `swe` of firn at `model.density_init_firn`, layered as an
accumulated column (near-surface layers at their maximum thicknesses, then
`mass_split` layers); its layers, ice substrate and `Tsrf` are at
`min(T_init, T0)`, its albedo is the dry-snow albedo and its snow age is zero. A
column with `swe = 0` stays bare, with its substrate and `Tsrf` at
`min(T_init, T0)`. `swe` requires `T_init`, and its depth `swe /
density_init_firn` must not exceed the snow-depth cap. The initial firn is
storage: the cumulative mass and flux diagnostics start at zero.
"""
function BESSIState(model::BESSIModel; swe=nothing, T_init=nothing)
    isnothing(swe) || !isnothing(T_init) ||
        error("`swe` requires `T_init`: give the temperature of the initial snowpack.")
    NF = number_type(model.c)
    ncol = ncols(model.grid)
    state = BESSIState(
        model.parameters,
        ncol,
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
        Vector{NF}(undef, ncol),
        Vector{NF}(undef, ncol),
        Matrix{NF}(undef, model.ice_substrate_layers, ncol),
    )
    _initialize_bessi_state_arrays!(state, model.density_init, model.temperature_init)
    isnothing(T_init) && return state
    return _initialize_bessi_snowpack!(state, swe, T_init)
end

@inline _cpu_layer_matrix(matrix::TransposedLayerMatrix) = Array(matrix)
@inline _cpu_layer_matrix(matrix) = Array(matrix)
@inline _cpu_state_array(state::BESSIState, name::Symbol) =
    _is_bessi_layer_field(name) ?
    _cpu_layer_matrix(getfield(state, name)) :
    Array(getfield(state, name))

function cpu_state(state::BESSIState)
    return BESSIState(
        state.parameters,
        state.ncol,
        (_cpu_state_array(state, name) for name in _BESSI_ARRAY_FIELD_NAMES)...,
    )
end

@inline _gpu_layer_matrix(matrix, storage_type) =
    TransposedLayerMatrix(adapt(storage_type, permutedims(Array(matrix), (2, 1))))
@inline _gpu_state_array(state::BESSIState, name::Symbol, storage_type) =
    _is_bessi_layer_field(name) ?
    _gpu_layer_matrix(getfield(state, name), storage_type) :
    adapt(storage_type, getfield(state, name))

function gpu_state(state::BESSIState, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return BESSIState(
        state.parameters,
        state.ncol,
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

const _PDD_STATE_FIELD_NAMES = (:snowpack_swe, :smb_ice, :runoff, :pdd_sum)
@inline get_fields(state::PDDState) = NamedTuple{_PDD_STATE_FIELD_NAMES}(
    ntuple(index -> getfield(state, _PDD_STATE_FIELD_NAMES[index]), Val(length(_PDD_STATE_FIELD_NAMES))),
)

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
Adapt.@adapt_structure PDDState

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

const _ITM_STATE_FIELD_NAMES = (
    :H_snow, :alb_s, :smb, :smbi, :melt, :runoff, :refreezing, :Tsrf,
    :melt_net, :smb_cum, :smb_ice, :melt_cum, :runoff_cum, :refreezing_cum,
)
@inline get_fields(state::ITMState) = NamedTuple{_ITM_STATE_FIELD_NAMES}(
    ntuple(index -> getfield(state, _ITM_STATE_FIELD_NAMES[index]), Val(length(_ITM_STATE_FIELD_NAMES))),
)

function ITMState(model::ITMModel)
    ncol = ncols(model.grid)
    return ITMState(
        fill(model.H_snow_max, ncol), fill(model.alb_snow_dry, ncol), zeros(ncol), zeros(ncol),
        zeros(ncol), zeros(ncol), zeros(ncol), fill(model.c.T0, ncol), zeros(ncol),
        zeros(ncol), zeros(ncol), zeros(ncol), zeros(ncol), zeros(ncol),
    )
end

initial_state(model::ITMModel) = ITMState(model)
Adapt.@adapt_structure ITMState

"""
    initial_state(model; kwargs...)

Initial state of `model`. For `BESSIModel` the keywords `swe` and `T_init`
start the columns with a firn column; see [`BESSIState`](@ref).
"""
initial_state(model::BESSIModel; kwargs...) = BESSIState(model; kwargs...)
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
