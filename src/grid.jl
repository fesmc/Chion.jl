#=
The discretisation: independent columns, their optional 2D layout, and the vertical layering of
the snow, firn and ice substrate in every column.
=#

"""
$(TYPEDEF)

Vertical discretisation of a column: up to `Ntot` snow and firn layers whose masses are kept
between `mass_min` and `mass_max` (a layer above `mass_max` splits into `mass_split` and the rest),
finer layers near the surface (`near_surface_layer_max_thicknesses_m`, the maximum thicknesses of
the top four layers), and a thermal ice substrate of `ice_substrate_layers` layers below the firn
whose thicknesses double downward from `ice_substrate_top_thickness_m`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct Layering{NF <: AbstractFloat}
    "maximum number of snow and firn layers (default `15`)"
    Ntot::Int = DEFAULT_NTOT
    "mass above which a layer splits (kg m⁻², default `500`)"
    mass_max::NF = DEFAULT_MASS_MAX
    "mass of the upper part of a split layer (kg m⁻², default `300`)"
    mass_split::NF = DEFAULT_MASS_SPLIT
    "mass below which a layer merges with its neighbour (kg m⁻², default `100`)"
    mass_min::NF = DEFAULT_MASS_MIN
    "maximum thicknesses of the top four layers (m, default `(0.02, 0.05, 0.10, 0.30)`)"
    near_surface_layer_max_thicknesses_m::NTuple{4, NF} = (0.02, 0.05, 0.10, 0.30)
    "number of ice-substrate layers below the firn (default `5`; `0` for snow on land)"
    ice_substrate_layers::Int = 5
    "thickness of the top ice-substrate layer (m, default `0.05`)"
    ice_substrate_top_thickness_m::NF = 0.05
end

function Layering(Ntot::Integer, mass_max::Real, mass_split::Real, mass_min::Real,
        near_surface_layer_max_thicknesses_m::NTuple{4, Real}, ice_substrate_layers::Integer,
        ice_substrate_top_thickness_m::Real)
    values = promote(map(float, (mass_max, mass_split, mass_min, near_surface_layer_max_thicknesses_m...,
        ice_substrate_top_thickness_m))...)
    return Layering(Int(Ntot), values[1:3]..., values[4:7], Int(ice_substrate_layers), values[8])
end

function validate(layering::Layering)
    _require(layering.Ntot >= 1, "`Ntot` must be at least 1.")
    _require(layering.mass_split < layering.mass_max, "`mass_split` must be smaller than `mass_max`.")
    _require(layering.mass_min < layering.mass_split, "`mass_min` must be smaller than `mass_split`.")
    _require(layering.mass_split / layering.mass_max >= 0.5, "`mass_split / mass_max` must be at least 0.5.")
    _require(all(>(0), layering.near_surface_layer_max_thicknesses_m),
        "`near_surface_layer_max_thicknesses_m` must contain positive values.")
    _require(layering.ice_substrate_layers >= 0, "`ice_substrate_layers` must be non-negative.")
    _require(layering.ice_substrate_top_thickness_m > 0, "`ice_substrate_top_thickness_m` must be positive.")
    return nothing
end

"""
$(TYPEDEF)

Positions of the columns on a regular x–y grid: column `k` is the cell `(js[k], is[k])` of the
`(length(y), length(x))` matrix `mask`.

# Fields
$(TYPEDFIELDS)
"""
struct SpatialLayout
    "x coordinates of the grid (m or degrees)"
    x::Vector{Float64}
    "y coordinates of the grid"
    y::Vector{Float64}
    "row (y index) of every column"
    js::Vector{Int}
    "column (x index) of every column"
    is::Vector{Int}
    "mask of the grid cells (1 where a column is)"
    mask::Matrix{Float64}
end

"""
$(TYPEDEF)

The discretisation of a model: `ncol` independent columns, their optional x–y `layout`, and the
vertical [`Layering`](@ref) of every column, at the floating-point type `NF` of the model.

    Grid(ncol; NF = Float64, layering = Layering(), x, y, js, is, mask)
    Grid(forcing::Forcing; NF = Float64, layering = Layering())

`x`, `y`, `js`, `is` (and optionally `mask`) give the [`SpatialLayout`](@ref); without them the
columns have no coordinates. `Grid(forcing)` takes the columns and layout of a file-backed forcing.

```jldoctest
julia> grid = Grid(4; layering = Layering(Ntot = 20));

julia> grid.Ntot, grid.ncol
(20, 4)
```

# Fields
$(TYPEDFIELDS)
"""
struct Grid{NF, L <: Union{Nothing, SpatialLayout}}
    "number of columns"
    ncol::Int
    "x–y positions of the columns, or `nothing`"
    layout::L
    "vertical discretisation of every column"
    layering::Layering{NF}
end

function Grid(ncol::Integer; NF::Type{<:AbstractFloat} = Float64, layering::Layering = Layering(),
        x = nothing, y = nothing, js = nothing, is = nothing, mask = nothing)
    _require(ncol > 0, "`ncol` must be positive.")
    validate(layering)
    layering = _promote_param(layering, NF)
    coordinates = (x, y, js, is)
    all(isnothing, coordinates) && return Grid(Int(ncol), nothing, layering)
    any(isnothing, coordinates) && throw(ArgumentError("Provide all of x, y, js, is for spatial coordinates."))
    x_values = Float64.(collect(x))
    y_values = Float64.(collect(y))
    js_values = Int.(collect(js))
    is_values = Int.(collect(is))
    _require(length(js_values) == ncol, "`js` length must equal `ncol`.")
    _require(length(is_values) == ncol, "`is` length must equal `ncol`.")
    mask_values = isnothing(mask) ? ones(Float64, length(y_values), length(x_values)) : Matrix{Float64}(mask)
    _require(size(mask_values) == (length(y_values), length(x_values)), "`mask` must have size (length(y), length(x)).")
    layout = SpatialLayout(x_values, y_values, js_values, is_values, mask_values)
    return Grid(Int(ncol), layout, layering)
end

"""The grid of `grid`'s columns and layout with the layering `layering` and number type `NF`."""
Grid(grid::Grid; NF::Type{<:AbstractFloat} = number_type(grid), layering::Layering = grid.layering) =
    (validate(layering); Grid(grid.ncol, grid.layout, _promote_param(layering, NF)))

const _GRID_LAYOUT_PROPERTIES = (:x, :y, :js, :is, :mask)

@inline function Base.getproperty(grid::Grid, name::Symbol)
    hasfield(Grid, name) && return getfield(grid, name)
    if name in _GRID_LAYOUT_PROPERTIES
        layout = getfield(grid, :layout)
        return isnothing(layout) ? nothing : getfield(layout, name)
    end
    return getfield(getfield(grid, :layering), name)
end

Base.propertynames(grid::Grid, private::Bool = false) =
    (fieldnames(Grid)..., _GRID_LAYOUT_PROPERTIES..., fieldnames(Layering)...)

number_type(::Grid{NF}) where {NF} = NF
ncols(grid::Grid) = grid.ncol
has_spatial_coords(::Grid{NF, Nothing}) where {NF} = false
has_spatial_coords(::Grid{NF, SpatialLayout}) where {NF} = true

function Base.show(io::IO, ::MIME"text/plain", grid::Grid{NF}) where {NF}
    println(io, "Grid{", NF, "}: ", grid.ncol, " columns", has_spatial_coords(grid) ?
        " on a $(length(grid.y))×$(length(grid.x)) x–y grid" : "")
    l = grid.layering
    print(io, "  layering: up to ", l.Ntot, " layers of ", l.mass_min, "–", l.mass_max, " kg m⁻², ",
        l.ice_substrate_layers, " ice-substrate layers")
    return nothing
end
