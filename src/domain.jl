"""Grid and backend helpers shared across Chion."""

@inline cuda_available() = CUDA.functional()
@inline gpu_storage_type() = CUDA.CuArray

@inline function _ka_backend(array)
    backend = KernelAbstractions.get_backend(array)
    return backend isa KernelAbstractions.CPU ? KernelAbstractions.CPU(; static=true) : backend
end

@inline function _wait_kernel(event)
    isnothing(event) || KernelAbstractions.wait(event)
    return nothing
end

struct SpatialLayout
    x::Vector{Float64}
    y::Vector{Float64}
    js::Vector{Int}
    is::Vector{Int}
    mask::Matrix{Float64}
end

"""Spatial discretization for a set of independent snowpack columns."""
struct SnowpackGrid{L<:Union{Nothing,SpatialLayout}}
    ncol::Int
    layout::L
end


@inline function Base.getproperty(grid::SnowpackGrid, name::Symbol)
    name === :ncol && return getfield(grid, :ncol)
    name === :layout && return getfield(grid, :layout)
    if name in (:x, :y, :js, :is, :mask)
        layout = getfield(grid, :layout)
        return isnothing(layout) ? nothing : getfield(layout, name)
    end
    return getfield(grid, name)
end

Base.propertynames(::SnowpackGrid, private::Bool=false) = private ?
    (:ncol, :layout, :x, :y, :js, :is, :mask) :
    (:ncol, :x, :y, :js, :is, :mask)

function SnowpackGrid(
    ncol::Integer;
    x=nothing,
    y=nothing,
    js=nothing,
    is=nothing,
    mask=nothing,
)
    ncol > 0 || error("`ncol` must be positive.")
    has_spatial = !isnothing(x) || !isnothing(y) || !isnothing(js) || !isnothing(is)
    if has_spatial
        (isnothing(x) || isnothing(y) || isnothing(js) || isnothing(is)) &&
            error("Provide all of x, y, js, is when supplying spatial coordinates.")
        x_v = Float64.(collect(x))
        y_v = Float64.(collect(y))
        js_v = Int.(collect(js))
        is_v = Int.(collect(is))
        length(js_v) == ncol || error("`js` length must equal `ncol`.")
        length(is_v) == ncol || error("`is` length must equal `ncol`.")
        mask_m = isnothing(mask) ? ones(Float64, length(y_v), length(x_v)) : Matrix{Float64}(mask)
        size(mask_m, 1) == length(y_v) || error("`mask` y-dimension must match `y`.")
        size(mask_m, 2) == length(x_v) || error("`mask` x-dimension must match `x`.")
        return SnowpackGrid(Int(ncol), SpatialLayout(x_v, y_v, js_v, is_v, mask_m))
    end
    return SnowpackGrid(Int(ncol), nothing)
end

ncols(grid::SnowpackGrid) = grid.ncol

has_spatial_coords(::SnowpackGrid{Nothing}) = false
has_spatial_coords(::SnowpackGrid{SpatialLayout}) = true

@inline function _validate_mass_partition(mass_max, mass_split, mass_min)
    mass_split < mass_max || error("`mass_split` must be smaller than `mass_max`.")
    mass_min < mass_split || error("`mass_min` must be smaller than `mass_split`.")
    mass_split / mass_max >= 0.5 || error("`mass_split / mass_max` must be at least 0.5.")
    return nothing
end

@inline function _domain_thresholds(::Type{NF}, mass_max, mass_split, mass_min) where {NF}
    _validate_mass_partition(mass_max, mass_split, mass_min)
    return convert(NF, mass_max), convert(NF, mass_split), convert(NF, mass_min)
end
