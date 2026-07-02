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

"""Spatial discretization for a set of independent snowpack columns."""
struct SnowpackGrid
    ncol::Int
    x::Union{Nothing, Vector{Float64}}
    y::Union{Nothing, Vector{Float64}}
    js::Union{Nothing, Vector{Int}}
    is::Union{Nothing, Vector{Int}}
    mask::Union{Nothing, Matrix{Float64}}
end

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
        return SnowpackGrid(Int(ncol), x_v, y_v, js_v, is_v, mask_m)
    end
    return SnowpackGrid(Int(ncol), nothing, nothing, nothing, nothing, nothing)
end

ncols(grid::SnowpackGrid) = grid.ncol

has_spatial_coords(grid::SnowpackGrid) =
    !isnothing(grid.x) && !isnothing(grid.y) && !isnothing(grid.js) &&
    !isnothing(grid.is) && !isnothing(grid.mask)

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
