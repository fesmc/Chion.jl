#=
The models: grid, prognostic state, cache (diagnostics, scratch, active columns), parameters and
the atmospheric forcing, as Laddie's `Model`. A model knows nothing about time beyond the loaded
forcing year; a `Simulation` advances it.
=#

abstract type AbstractSnowModel end

# ---------------------------------------------------------------------------
# Atmospheric forcing in the model
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

The forcing of a model: its [`Forcing`](@ref) `source` and the loaded calendar year, as the host
[`SnowpackForcing`](@ref) `data` and the kernel view `fields` on the model's backend in the model's
number type. [`update_forcing!`](@ref) loads another year.

# Fields
$(TYPEDFIELDS)
"""
mutable struct AtmosphericForcing{S <: Forcing, D <: SnowpackForcing, F}
    "where the forcing comes from, one calendar year at a time"
    source::S
    "the loaded calendar year"
    year::Int
    "the loaded year on the host"
    data::D
    "the loaded year as the kernels read it, on the model's backend"
    fields::F
end

_as_forcing(forcing::Forcing) = forcing
_as_forcing(forcing::SnowpackForcing) = Forcing(forcing)

# The forcing fields a model's kernel reads, in its number type, on `backend`.
_kernel_forcing(::Type{<:Params}, data) = get_fields(data)
_kernel_forcing(::Type{<:PDDParams}, data) = PDDForcing(data)
_kernel_forcing(::Type{<:ITMParams}, data) = ITMForcing(data)

_forcing_on(::KernelAbstractions.CPU, fields) = fields
_forcing_on(backend, fields) = _forcing_to_device(backend, fields)
_forcing_to_device(backend, fields::NamedTuple) = adapt(backend, fields)
# PDD and ITM keep `dt_days` on the host: their kernels receive the step length as a scalar.
_forcing_to_device(backend, f::PDDForcing) =
    PDDForcing(f.dt_days, adapt(backend, f.air_temperature), adapt(backend, f.snowfall_rate), adapt(backend, f.rainfall_rate))
_forcing_to_device(backend, f::ITMForcing) = ITMForcing(f.dt_days,
    (adapt(backend, getfield(f, name)) for name in fieldnames(ITMForcing)[2:end])...)

function _kernel_fields(params, data, backend)
    fields = forcing_fields_as(number_type(params), _kernel_forcing(typeof(params), data))
    return _forcing_on(backend, fields)
end

function AtmosphericForcing(source::Forcing, year::Integer, params, backend)
    data = load_forcing_year(source, year)
    return AtmosphericForcing(source, Int(year), data, _kernel_fields(params, data, backend))
end

"""
    update_forcing!(model, year)

Load the calendar year `year` of the model's forcing source, unless it is loaded already.
"""
function update_forcing!(model::AbstractSnowModel, year::Integer)
    forcing = model.forcing
    forcing.year == year && return model
    data = load_forcing_year(forcing.source, year)
    _validate_forcing(model.params, data)
    fields = _kernel_fields(model.params, data, _backend(model))
    typeof(data) === typeof(forcing.data) && typeof(fields) === typeof(forcing.fields) ||
        error("The forcing of $year has other field types than that of $(forcing.year) " *
              "(e.g. a field constant in one year and varying in the other).")
    forcing.year = Int(year)
    forcing.data = data
    forcing.fields = fields
    return model
end

"""
    sync_forcing!(model)

Rebuild the kernel view of the loaded forcing year from its host data, after the fields of
`model.forcing.data` (a [`SnowpackForcing`](@ref)) were changed in place, e.g. by a coupled
ice-sheet model.
"""
function sync_forcing!(model::AbstractSnowModel)
    forcing = model.forcing
    forcing.fields = _kernel_fields(model.params, forcing.data, _backend(model))
    return model
end

"""
    set_input!(model, name, values)

Replace the forcing field `name` (e.g. `:surface_height` from an ice-sheet model) of the loaded
year by `values`: a column vector (constant in time) or a `(column, time)` matrix. The air pressure
derived from the surface height is not recomputed.
"""
function set_input!(model::AbstractSnowModel, name::Symbol, values)
    forcing = model.forcing
    field = getproperty(forcing.data, name)
    _set_forcing_values!(field, values)
    forcing.fields = _kernel_fields(model.params, forcing.data, _backend(model))
    return model
end

_set_forcing_values!(field::ColumnForcingMatrix, values::AbstractVector) = (_set_column_values!(field, values); field)
_set_forcing_values!(field::AbstractMatrix, values) = (field .= values; field)
_set_forcing_values!(field, values) =
    error("Cannot set a forcing field stored as $(nameof(typeof(field))); give the field per column in the forcing.")

# ---------------------------------------------------------------------------
# Validation of the forcing against the model
# ---------------------------------------------------------------------------

function _validate_forcing(params::Params, forcing::SnowpackForcing)
    _validate_diurnal_forcing(params.diurnal_cycle, forcing)
    _validate_albedo_forcing(params.albedo_scheme, forcing)
    return nothing
end

_validate_diurnal_forcing(::NoDiurnalCycle, forcing) = nothing
_validate_diurnal_forcing(::DiurnalSubsteps, forcing) = all(isfinite, forcing.latitude_deg) ||
    throw(ArgumentError("The forcing needs `latitude_deg` for the diurnal substeps; provide it or use " *
        "`Params(diurnal_cycle = NoDiurnalCycle())` (e.g. for sub-daily forcing)."))

_validate_albedo_forcing(::AbstractAlbedoScheme, forcing) = nothing
_validate_albedo_forcing(::PrescribedAlbedo, forcing) = all(forcing.has_prescribed_albedo) ||
    throw(ArgumentError("`PrescribedAlbedo` needs `prescribed_albedo` for every column and step."))

_validate_forcing(::PDDParams, forcing::SnowpackForcing) = nothing

function _validate_forcing(::ITMParams, forcing::SnowpackForcing)
    for name in (:surface_height, :ice_thickness, :annual_pdd, :latitude_deg)
        all(isfinite, getproperty(forcing, name)) || throw(ArgumentError("ITMModel needs `$name` in the forcing."))
    end
    all(>=(0.0), forcing.ice_thickness) || throw(ArgumentError("`ice_thickness` must be non-negative for ITMModel."))
    all(>=(0.0), forcing.annual_pdd) || throw(ArgumentError("`annual_pdd` must be non-negative for ITMModel."))
    return nothing
end

# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

Diagnostic and scratch arrays of a [`BESSIModel`](@ref) and the columns that are stepped.

# Fields
$(TYPEDFIELDS)
"""
mutable struct BESSICache{W <: ColumnarStepWorkspace, V <: AbstractVector, A <: AbstractVector{Int}}
    "scratch of the column step (heat solve) and the cumulative `boundary_enthalpy` (J m⁻²)"
    workspace::W
    "snow and firn depth (m), from [`update_diagnostics!`](@ref)"
    thickness::V
    "solid plus liquid mass of the snow and firn (kg m⁻²), from `update_diagnostics!`"
    wet_mass::V
    "bulk density of the snow and firn (kg m⁻³), from `update_diagnostics!`"
    bulk_density::V
    "liquid water in the snow and firn (kg m⁻²), from `update_diagnostics!`"
    liquid_water::V
    "stepped columns (host)"
    active::Vector{Bool}
    "indices of the stepped columns, on the model's backend"
    active_indices::A
end

function BESSICache(state::BESSIState, layering::Layering)
    zeros_like() = fill!(similar(state.runoff), zero(eltype(state.runoff)))
    ncol = ncols(state)
    return BESSICache(ColumnarStepWorkspace(state, layering), zeros_like(), zeros_like(), zeros_like(),
        zeros_like(), fill(true, ncol), _active_indices(state.runoff, 1:ncol))
end

"""
$(TYPEDEF)

The columns that a bulk model ([`PDDModel`](@ref), [`ITMModel`](@ref)) steps.

# Fields
$(TYPEDFIELDS)
"""
mutable struct BulkCache{A <: AbstractVector{Int}}
    "stepped columns (host)"
    active::Vector{Bool}
    "indices of the stepped columns, on the model's backend"
    active_indices::A
end

BulkCache(state) = BulkCache(fill(true, ncols(state)), _active_indices(state.runoff, 1:ncols(state)))

_cache(state::BESSIState, grid::Grid) = BESSICache(state, grid.layering)
_cache(state, grid::Grid) = BulkCache(state)

# Column indices on the backend of `reference`.
_active_indices(reference, indices) = _to_device(KernelAbstractions.get_backend(reference), collect(Int, indices))

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

const _MODEL_FIELDS = (:grid, :state, :cache, :params, :boundary, :forcing)
const _RESERVED_MODEL_PROPERTIES = (_MODEL_FIELDS..., :NF, :layering, :year)

# Flat forwarding resolves `model.x` through state, cache, params, grid and forcing data, so a name
# in two of them would be shadowed silently: reject it at construction.
function _check_property_collisions(grid, state, cache, params, forcing)
    seen = Dict{Symbol, String}()
    for (label, names) in (("Grid", propertynames(grid)), ("State", fieldnames(typeof(state))),
            ("Cache", fieldnames(typeof(cache))), ("Params", fieldnames(typeof(params))))
        for name in names
            name in _RESERVED_MODEL_PROPERTIES && name !== :layering &&
                error("Field `$name` of $label collides with the model property `$name`.")
            haskey(seen, name) && error("Model property forwarding is ambiguous: `$name` exists in $(seen[name]) and $label.")
            seen[name] = label
        end
    end
    return nothing
end

"""
$(TYPEDEF)

The layered BESSI snowpack and firn model: a [`Grid`](@ref), the [`BESSIState`](@ref), the
[`BESSICache`](@ref), the [`Params`](@ref) and the [`AtmosphericForcing`](@ref). Wrap it in a
[`Simulation`](@ref) to run it.

    BESSIModel(grid; forcing, params = Params(), boundary = BoundaryConditions(), backend = CPU(), state = nothing)

`forcing` is a [`Forcing`](@ref) or an in-memory [`SnowpackForcing`](@ref); its first year is
loaded. `backend` is a KernelAbstractions backend (`CPU()`, `CUDABackend()`) or its name
(`:threads`, `:gpu`). `state` (a host `BESSIState`, e.g. a restart) replaces the initial state of `params`.
Fields resolve through the model: `model.mass` (state), `model.thickness` (cache), `model.alpha_dry`
(params), `model.ncol`, `model.Ntot` (grid), `model.year` (loaded forcing year).

```jldoctest
julia> grid = Grid(2);

julia> forcing = SnowpackForcing(dt_days = ones(3), ncol = 2, air_temperature_c = -5.0, snowfall_mm_day = 2.0,
           rainfall_mm_day = 0.0, shortwave_down = 100.0, latitude_deg = 70.0);

julia> model = BESSIModel(grid; forcing, params = Params(albedo_scheme = ConstantAlbedo()));

julia> model.albedo_scheme, size(model.mass)
(ConstantAlbedo(), (15, 2))
```

# Fields
$(TYPEDFIELDS)
"""
mutable struct BESSIModel{G <: Grid, S <: BESSIState, C <: BESSICache, P <: Params, B <: BoundaryConditions,
        F <: AtmosphericForcing} <: AbstractSnowModel
    "columns, layout and layering"
    grid::G
    "prognostic arrays"
    state::S
    "diagnostics, scratch and active columns"
    cache::C
    "physical constants and schemes"
    params::P
    "surface and bottom boundary conditions of the heat solve"
    boundary::B
    "the forcing source and the loaded year"
    forcing::F

    function BESSIModel{G, S, C, P, B, F}(grid, state, cache, params, boundary, forcing) where {G, S, C, P, B, F}
        _check_property_collisions(grid, state, cache, params, forcing)
        return new{G, S, C, P, B, F}(grid, state, cache, params, boundary, forcing)
    end
end

BESSIModel(grid::G, state::S, cache::C, params::P, boundary::B, forcing::F) where {G, S, C, P, B, F} =
    BESSIModel{G, S, C, P, B, F}(grid, state, cache, params, boundary, forcing)

"""
$(TYPEDEF)

Bulk positive-degree-day surface mass balance: a one-layer snow reservoir capped at `H_snow_max`,
melt from the [`PDDParams`](@ref).

    PDDModel(grid; forcing, params = PDDParams(), backend = CPU(), state = nothing)

# Fields
$(TYPEDFIELDS)
"""
mutable struct PDDModel{G <: Grid, S <: PDDState, C <: BulkCache, P <: PDDParams, F <: AtmosphericForcing} <:
               AbstractSnowModel
    "columns and layout"
    grid::G
    "prognostic arrays"
    state::S
    "active columns"
    cache::C
    "degree-day factors and constants"
    params::P
    "the forcing source and the loaded year"
    forcing::F

    function PDDModel{G, S, C, P, F}(grid, state, cache, params, forcing) where {G, S, C, P, F}
        _check_property_collisions(grid, state, cache, params, forcing)
        return new{G, S, C, P, F}(grid, state, cache, params, forcing)
    end
end

"""
$(TYPEDEF)

Bulk insolation–temperature–melt surface mass balance, ported from the Fortran `fesmc/chion`/smbpal
code. Needs `surface_height`, `ice_thickness`, `annual_pdd` and `latitude_deg` in the forcing;
`q_sw_net` is the insolation driver where it is given, otherwise `shortwave_down`.

    ITMModel(grid; forcing, params = ITMParams(), backend = CPU(), state = nothing)

# Fields
$(TYPEDFIELDS)
"""
mutable struct ITMModel{G <: Grid, S <: ITMState, C <: BulkCache, P <: ITMParams, F <: AtmosphericForcing} <:
               AbstractSnowModel
    "columns and layout"
    grid::G
    "prognostic arrays"
    state::S
    "active columns"
    cache::C
    "ITM coefficients and constants"
    params::P
    "the forcing source and the loaded year"
    forcing::F

    function ITMModel{G, S, C, P, F}(grid, state, cache, params, forcing) where {G, S, C, P, F}
        _check_property_collisions(grid, state, cache, params, forcing)
        return new{G, S, C, P, F}(grid, state, cache, params, forcing)
    end
end

for M in (:PDDModel, :ITMModel)
    @eval $M(grid::G, state::S, cache::C, params::P, forcing::F) where {G, S, C, P, F} =
        $M{G, S, C, P, F}(grid, state, cache, params, forcing)
end

const _DEFAULT_PARAMS = (BESSIModel = Params, PDDModel = PDDParams, ITMModel = ITMParams)

for M in (:BESSIModel, :PDDModel, :ITMModel)
    @eval function $M(grid::Grid; forcing, params = $(_DEFAULT_PARAMS[M])(; NF = number_type(grid)),
            backend = KernelAbstractions.CPU(), state = nothing, year::Union{Nothing, Integer} = nothing,
            kwargs...)
        return _build_model($M, grid, _as_forcing(forcing), params, backend, state, year; kwargs...)
    end
end

function _build_model(M, grid::Grid, source::Forcing, params, backend_choice, state, year; kwargs...)
    backend = _ka_backend_from(backend_choice)
    number_type(grid) === number_type(params) || throw(ArgumentError(
        "The grid is in $(number_type(grid)) and the parameters in $(number_type(params)); " *
        "build both with the same `NF`."))
    ncols(forcing_grid(source)) == ncols(grid) || throw(ArgumentError(
        "The grid has $(ncols(grid)) columns and the forcing $(ncols(forcing_grid(source)))."))
    first_year = something(year, first(forcing_years(source)))
    host_state = state === nothing ? initial_state(grid, params) : state
    _check_state(host_state, grid, params)
    device_state = to_backend(host_state, backend)
    forcing = AtmosphericForcing(source, first_year, params, backend)
    _validate_forcing(params, forcing.data)
    return _assemble(M, grid, device_state, _cache(device_state, grid), params, forcing; kwargs...)
end

# The model from its parts; only BESSI has boundary conditions (in the model's number type).
function _assemble(::Type{BESSIModel}, grid, state, cache, params, forcing;
        boundary::BoundaryConditions = BoundaryConditions())
    validate(boundary)
    return BESSIModel(grid, state, cache, params, _promote_param(boundary, number_type(params)), forcing)
end
_assemble(M, grid, state, cache, params, forcing) = M(grid, state, cache, params, forcing)

"""`model` with some fields replaced (the others shared), e.g. another forcing or state."""
function _replace(model::M; kwargs...) where {M <: AbstractSnowModel}
    values = map(name -> get(kwargs, name, getfield(model, name)), fieldnames(M))
    return Base.typename(M).wrapper(values...)
end

function _check_state(state::BESSIState, grid::Grid, params::Params)
    ncols(state) == ncols(grid) || throw(ArgumentError("The state has $(ncols(state)) columns, the grid $(ncols(grid))."))
    size(state.mass, 1) == grid.layering.Ntot && size(state.ice_temperature, 1) == grid.layering.ice_substrate_layers ||
        throw(ArgumentError("The state does not have the layers of the grid's layering."))
    eltype(state.mass) === number_type(params) ||
        throw(ArgumentError("The state is in $(eltype(state.mass)), the parameters in $(number_type(params))."))
    return nothing
end
_check_state(state::PDDState, grid, params) = _check_vector_state(state, grid, params)
_check_state(state::ITMState, grid, params) = _check_vector_state(state, grid, params)
function _check_vector_state(state, grid, params)
    ncols(state) == ncols(grid) || throw(ArgumentError("The state has $(ncols(state)) columns, the grid $(ncols(grid))."))
    eltype(state.runoff) === number_type(params) ||
        throw(ArgumentError("The state is in $(eltype(state.runoff)), the parameters in $(number_type(params))."))
    return nothing
end

"""
    BESSIModel(forcing::Forcing; NF = Float64, layering = Layering(), kwargs...)

The model on the columns of `forcing`: `BESSIModel(Grid(forcing; NF, layering); forcing, kwargs...)`
with `params = Params(; NF)` unless given. Likewise `PDDModel(forcing)` and `ITMModel(forcing)`.
"""
function BESSIModel(forcing::Forcing; NF::Type{<:Real} = Float64, layering::Layering = Layering(),
        params = Params(; NF), kwargs...)
    return BESSIModel(Grid(forcing; NF = number_type(params), layering); forcing, params, kwargs...)
end
PDDModel(forcing::Forcing; NF::Type{<:Real} = Float64, params = PDDParams(; NF), kwargs...) =
    PDDModel(Grid(forcing; NF = number_type(params)); forcing, params, kwargs...)
ITMModel(forcing::Forcing; NF::Type{<:Real} = Float64, params = ITMParams(; NF), kwargs...) =
    ITMModel(Grid(forcing; NF = number_type(params)); forcing, params, kwargs...)

"""
    Grid(forcing::Forcing; NF = Float64, layering = Layering())

The grid of the columns (and their x–y layout) of a file-backed or in-memory forcing.
"""
function Grid(forcing::Forcing; NF::Type{<:Real} = Float64, layering::Layering = Layering())
    g = forcing_grid(forcing)
    return Grid(g.ncol, g.layout, (validate(layering); _promote_param(layering, NF)))
end

# ---------------------------------------------------------------------------
# Flat property forwarding
# ---------------------------------------------------------------------------

@inline function Base.getproperty(model::AbstractSnowModel, name::Symbol)
    hasfield(typeof(model), name) && return getfield(model, name)
    name === :NF && return number_type(getfield(model, :params))
    name === :year && return getfield(model, :forcing).year
    state = getfield(model, :state)
    hasfield(typeof(state), name) && return getfield(state, name)
    cache = getfield(model, :cache)
    hasfield(typeof(cache), name) && return getfield(cache, name)
    params = getfield(model, :params)
    hasfield(typeof(params), name) && return getfield(params, name)
    grid = getfield(model, :grid)
    name in propertynames(grid) && return getproperty(grid, name)
    error("$(nameof(typeof(model))) has no property `$name`.")
end

function Base.propertynames(model::AbstractSnowModel, private::Bool = false)
    names(x) = fieldnames(typeof(x))
    return (fieldnames(typeof(model))..., :NF, :year, names(model.state)..., names(model.cache)..., names(model.params)...,
        propertynames(model.grid)...)
end

number_type(model::AbstractSnowModel) = number_type(model.params)
ncols(model::AbstractSnowModel) = ncols(model.grid)
"""How a model's kernels run, from the type of its state arrays: `NativeExecution()` launches
KernelAbstractions kernels; the Reactant extension adds a compiled execution."""
struct NativeExecution end
_execution(array) = NativeExecution()

_backend(model::AbstractSnowModel) = _backend(_execution(model.state.runoff), model)
_backend(::NativeExecution, model) = KernelAbstractions.get_backend(model.state.runoff)

# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

"""
    to_backend(model, backend)

A new model with the state, cache and loaded forcing on the KernelAbstractions `backend`
(`CPU()`, or `CUDABackend()` after `using CUDA`); the original is not modified. The grid, the
parameters and the forcing source are shared.

```julia
using CUDA
model = to_backend(BESSIModel(grid; forcing), CUDABackend())
```
"""
to_backend(model::AbstractSnowModel, backend::Symbol) = to_backend(model, _ka_backend_from(backend))

function to_backend(model::M, backend) where {M <: AbstractSnowModel}
    state = to_backend(model.state, backend)
    cache = _cache(state, model.grid)
    _copy_cache!(cache, model.cache)
    cache.active = copy(model.cache.active)
    cache.active_indices = _active_indices(state.runoff, findall(cache.active))
    source = model.forcing
    forcing = AtmosphericForcing(source.source, source.year, source.data, _kernel_fields(model.params, source.data, backend))
    return _replace(model; state, cache, forcing)
end

_copy_cache!(dest::BulkCache, src::BulkCache) = dest
function _copy_cache!(dest::BESSICache, src::BESSICache)
    copyto!(dest.workspace.boundary_enthalpy, Array(src.workspace.boundary_enthalpy))
    for name in (:thickness, :wet_mass, :bulk_density, :liquid_water)
        copyto!(getfield(dest, name), Array(getfield(src, name)))
    end
    return dest
end

"""
    cpu_state(model) -> BESSIState | PDDState | ITMState

A copy of the model's state on the CPU (layers `(layer, column)`), e.g. to inspect or save the
state of a GPU run.
"""
cpu_state(model::AbstractSnowModel) = to_backend(model.state, KernelAbstractions.CPU())

"""
    set_state!(model, state)

Copy the host `state` (e.g. from [`cpu_state`](@ref) or a restart) into the model.
"""
function set_state!(model::AbstractSnowModel, state)
    _check_state(state, model.grid, model.params)
    device_state = to_backend(state, _backend(model))
    for name in fieldnames(typeof(device_state))
        _copy_array!(getfield(model.state, name), getfield(device_state, name))
    end
    return model
end

_copy_array!(dest::TransposedLayerMatrix, src::TransposedLayerMatrix) = (copyto!(dest.parent, src.parent); dest)
_copy_array!(dest, src) = copyto!(dest, src)

# ---------------------------------------------------------------------------
# Active columns
# ---------------------------------------------------------------------------

@inline _active_mask_value(value::Bool) = value
@inline _active_mask_value(value) = isfinite(Float64(value)) && Float64(value) != 0.0

function _active_mask_vector(mask, grid, ncol::Int)
    data = collect(mask)
    if ndims(data) == 1
        length(data) == ncol || throw(ArgumentError("The mask has $(length(data)) values for $ncol columns."))
        return [_active_mask_value(data[idx]) for idx in eachindex(data)]
    elseif ndims(data) == 2
        has_spatial_coords(grid) || throw(ArgumentError("A 2-D mask needs a grid with spatial coordinates."))
        size(data) == size(grid.mask) ||
            throw(ArgumentError("The mask must have size $(size(grid.mask)), got $(size(data))."))
        return [_active_mask_value(data[grid.js[col], grid.is[col]]) for col in 1:ncol]
    end
    throw(ArgumentError("The mask must be a column vector or a y–x matrix matching `grid.mask`."))
end

"""
    set_active_mask!(model, mask; reset_newly_inactive = true)

Step only the columns where `mask` is true (a column vector or a y–x matrix on the grid). Columns
that become inactive are reset to the initial state unless `reset_newly_inactive = false`.
"""
function set_active_mask!(model::AbstractSnowModel, mask; reset_newly_inactive::Bool = true)
    active = _active_mask_vector(mask, model.grid, ncols(model))
    any(active) || throw(ArgumentError("The mask keeps no column."))
    newly_inactive = findall(model.cache.active .& .!active)
    reset_newly_inactive && !isempty(newly_inactive) &&
        _reset_columns!(model, _active_indices(model.state.runoff, newly_inactive))
    model.cache.active = active
    model.cache.active_indices = _active_indices(model.state.runoff, findall(active))
    return model
end

_reset_columns!(model::BESSIModel, indices) = _reset_columns!(_execution(model.state.runoff), model, indices)
_reset_columns!(::NativeExecution, model::BESSIModel, indices) = _initialize_columns!(model.state, model.params, indices)

@kernel function _reset_pdd_columns_kernel!(fields, indices)
    item = @index(Global)
    @inbounds begin
        idx = indices[item]
        fields.snowpack_swe[idx] = zero(eltype(fields.snowpack_swe))
        fields.smb_ice[idx] = zero(eltype(fields.smb_ice))
        fields.runoff[idx] = zero(eltype(fields.runoff))
        fields.pdd_sum[idx] = zero(eltype(fields.pdd_sum))
    end
end

_reset_columns!(model::PDDModel, indices) =
    launch!(_reset_pdd_columns_kernel!, model.state.runoff, length(indices), get_fields(model.state), indices)

@kernel function _reset_itm_columns_kernel!(fields, indices, H_snow_init, albedo_init, T0)
    item = @index(Global)
    @inbounds begin
        idx = indices[item]
        zero_value = zero(H_snow_init)
        fields.H_snow[idx] = H_snow_init
        fields.alb_s[idx] = albedo_init
        fields.smb[idx] = zero_value
        fields.smbi[idx] = zero_value
        fields.melt[idx] = zero_value
        fields.runoff[idx] = zero_value
        fields.refreezing[idx] = zero_value
        fields.Tsrf[idx] = T0
        fields.melt_net[idx] = zero_value
        fields.smb_cum[idx] = zero_value
        fields.smb_ice[idx] = zero_value
        fields.melt_cum[idx] = zero_value
        fields.runoff_cum[idx] = zero_value
        fields.refreezing_cum[idx] = zero_value
    end
end

_reset_columns!(model::ITMModel, indices) = launch!(_reset_itm_columns_kernel!, model.state.runoff, length(indices),
    get_fields(model.state), indices, model.params.H_snow_max, model.params.alb_snow_dry, model.params.T0)

# ---------------------------------------------------------------------------
# Stepping the loaded forcing year
# ---------------------------------------------------------------------------

"""Number of forcing steps of the loaded year."""
forcing_steps(model::AbstractSnowModel) = length(model.forcing.data.time_values)

"""
    _step!(model, time_start, time_stop)

Advance the active columns over the steps `time_start:time_stop` of the loaded forcing year.
"""
_step!(model::BESSIModel, time_start::Int, time_stop::Int) =
    _step!(_execution(model.state.runoff), model, time_start, time_stop)
_step!(::NativeExecution, model::BESSIModel, time_start::Int, time_stop::Int) =
    _step_columns!(model, time_start, time_stop)

function _step!(model::PDDModel, time_start::Int, time_stop::Int)
    for time_index in time_start:time_stop
        _pdd_step_arrays!(model.state, model.forcing.fields, time_index, model.params, model.cache.active_indices)
    end
    return nothing
end

function _step!(model::ITMModel, time_start::Int, time_stop::Int)
    for time_index in time_start:time_stop
        _itm_step_arrays!(model.state, model.forcing.fields, time_index, model.params, model.cache.active_indices)
    end
    return nothing
end

"""
    integrate!(model, steps)

Advance `model` over the forcing steps `steps` (a range of step indices of the loaded forcing
year) without a [`Simulation`](@ref): no output and no change of the forcing year. The entry
point for automatic differentiation.
"""
function integrate!(model::AbstractSnowModel, steps::AbstractUnitRange{<:Integer})
    first(steps) >= 1 && last(steps) <= forcing_steps(model) ||
        throw(ArgumentError("Steps $steps outside the $(forcing_steps(model)) steps of the loaded year."))
    _step!(model, Int(first(steps)), Int(last(steps)))
    return model
end

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

"""
    update_diagnostics!(model::BESSIModel)

Recompute the column diagnostics `thickness`, `wet_mass`, `bulk_density` and `liquid_water` of the
cache from the layers. [`run!`](@ref) calls it at the end of every run.
"""
update_diagnostics!(model::BESSIModel) = _update_diagnostics!(_execution(model.state.runoff), model)

function _update_diagnostics!(::NativeExecution, model::BESSIModel)
    cache = model.cache
    launch!(_summarize_column_kernel!, model.state.runoff, ncols(model),
        (thickness = cache.thickness, wet_mass = cache.wet_mass, bulk_density = cache.bulk_density,
         liquid_water = cache.liquid_water), get_fields(model.state))
    return model
end
update_diagnostics!(model::AbstractSnowModel) = model

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", model::AbstractSnowModel)
    println(io, nameof(typeof(model)), "{", number_type(model), "} on ", nameof(typeof(_backend(model))))
    println(io, "  grid:    ", ncols(model), " columns", has_spatial_coords(model.grid) ? " (x–y layout)" : "",
        model isa BESSIModel ? ", up to $(model.grid.layering.Ntot) layers" : "")
    println(io, "  forcing: ", model.forcing.source.description, ", year ", model.forcing.year, " loaded")
    active = count(model.cache.active)
    active < ncols(model) && println(io, "  active:  ", active, " of ", ncols(model), " columns")
    print(io, "  params:  ", _params_summary(model.params))
    return nothing
end

_params_summary(params::Params) = join(("$(name) = $(getfield(params, name))" for name in
    (:albedo_scheme, :densification_scheme, :turbulent_flux_scheme, :longwave_scheme, :diurnal_cycle)), ", ")
_params_summary(params::PDDParams) = "method = $(params.method)"
_params_summary(params::ITMParams) = "ITM"
