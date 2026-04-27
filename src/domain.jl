"""
Abstract types shared across Chion's public simulation-first API.
"""

abstract type AbstractSnowpackGrid end
abstract type AbstractSnowModel{G <: AbstractSnowpackGrid} end

ncols(g::AbstractSnowpackGrid) = g.ncol
grid(m::AbstractSnowModel) = m.grid

"""
    AbstractSnowpackDomain{NF}

Abstract supertype for array-backed snowpack state containers that expose the
core layer arrays and auxiliary diagnostics used by the process kernels.
"""
abstract type AbstractSnowpackDomain{NF} end

"""
    cuda_available()

Return `true` when CUDA is functional in the current Julia session. This does
not allocate any device arrays; it only checks backend availability.
"""
@inline cuda_available() = CUDA.functional()

"""
    gpu_storage_type()

Return the default array storage type used for GPU adaptation in the current
build.
"""
@inline gpu_storage_type() = CUDA.CuArray

"""
    _ka_backend(array)

Return the KernelAbstractions backend associated with `array`. The result is
used to launch backend-specific kernels for CPU or GPU storage.
"""
@inline function _ka_backend(array)
    return KernelAbstractions.get_backend(array)
end

"""
    _wait_kernel(event)

Synchronize a KernelAbstractions launch event when one was returned. The
function is a no-op for `nothing` and always returns `nothing`.
"""
@inline function _wait_kernel(event)
    if !isnothing(event)
        KernelAbstractions.wait(event)
    end
    return nothing
end

"""Spatial discretization for a set of independent snowpack columns."""
struct SnowpackGrid{DEV} <: AbstractSnowpackGrid
    device::DEV
    ncol::Int
    x::Union{Nothing, Vector{Float64}}
    y::Union{Nothing, Vector{Float64}}
    js::Union{Nothing, Vector{Int}}
    is::Union{Nothing, Vector{Int}}
    mask::Union{Nothing, Matrix{Float64}}
end

function SnowpackGrid(
    device,
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
        return SnowpackGrid(device, Int(ncol), x_v, y_v, js_v, is_v, mask_m)
    end
    return SnowpackGrid(device, Int(ncol), nothing, nothing, nothing, nothing, nothing)
end

has_spatial_coords(g::SnowpackGrid) =
    !isnothing(g.x) && !isnothing(g.y) && !isnothing(g.js) && !isnothing(g.is) && !isnothing(g.mask)
has_spatial_coords(::Nothing) = false

"""
Snowpack domain state container.
"""

Base.eltype(::AbstractSnowpackDomain{NF}) where {NF} = NF

"""
    column_count(domain)

Return the number of snow columns stored in `domain`.
"""
@inline column_count(domain::AbstractSnowpackDomain) = domain.ncol

"""
    _validate_mass_partition(mass_max, mass_split, mass_min)

Check that the layer split thresholds are ordered consistently for the layer
management routines.
"""
@inline function _validate_mass_partition(mass_max, mass_split, mass_min)
    mass_split < mass_max || error("`mass_split` must be smaller than `mass_max`.")
    mass_min < mass_split || error("`mass_min` must be smaller than `mass_split`.")
    mass_split / mass_max >= 0.5 || error("`mass_split / mass_max` must be at least 0.5.")
    return nothing
end

"""
    _validate_domain_vector(name, values, ncol)

Validate that a vector-valued state field has one entry per column.
"""
@inline function _validate_domain_vector(name::AbstractString, values, ncol::Int)
    length(values) == ncol || error("`$name` must match `N`.")
    return nothing
end

"""
    SnowpackDomain{NF,...}

Mutable array-backed snowpack state for one or more columns. Each column stores
layer-wise solid mass, liquid water, density, and temperature together with
auxiliary diagnostics such as runoff, snow cover, and surface albedo.
"""
mutable struct SnowpackDomain{
        NF <: AbstractFloat,
        NI <: AbstractVector{<:Integer},
        MT <: AbstractMatrix{NF},
        VT <: AbstractVector{NF},
    } <: AbstractSnowpackDomain{NF}
    c::SnowpackPhysicalConstants{NF}
    Ntot::Int
    ncol::Int
    mass_max::NF
    mass_split::NF
    mass_min::NF
    rho_max::NF
    N::NI
    mass::MT
    mass_w::MT
    density::MT
    temperature::MT
    mass_base::VT
    smb_ice::VT
    runoff::VT
    Tsrf::VT
    snow_cover::VT
    albedo_dynamic::VT
end

"""
    SnowpackDomain(; c=SnowpackPhysicalConstants(), Ntot=DEFAULT_NTOT, ncol=1, ...)

Allocate a new array-backed snowpack domain with `ncol` columns and `Ntot`
maximum layers per column. State arrays are initialized to simple defaults and
mutated in-place by the model.
"""
function SnowpackDomain(;
    c::SnowpackPhysicalConstants=SnowpackPhysicalConstants(),
    Ntot::Int=DEFAULT_NTOT,
    ncol::Int=1,
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    rho_max::Real=DEFAULT_RHO_MAX,
    density_init::Real=DEFAULT_DENSITY_INIT,
    temperature_init::Real=DEFAULT_TEMPERATURE_INIT,
)
    ncol > 0 || error("`ncol` must be positive.")
    _validate_mass_partition(mass_max, mass_split, mass_min)

    NF = number_type(c)
    return SnowpackDomain(
        c,
        Ntot,
        ncol,
        convert(NF, mass_max),
        convert(NF, mass_split),
        convert(NF, mass_min),
        convert(NF, rho_max),
        zeros(Int, ncol),
        zeros(NF, Ntot, ncol),
        zeros(NF, Ntot, ncol),
        fill(convert(NF, density_init), Ntot, ncol),
        fill(convert(NF, temperature_init), Ntot, ncol),
        zeros(NF, ncol),
        zeros(NF, ncol),
        zeros(NF, ncol),
        fill(c.T0, ncol),
        zeros(NF, ncol),
        fill(c.alpha_dry, ncol),
    )
end

"""
    SnowpackDomain(N, mass, mass_w, density, temperature, mass_base, smb_ice, runoff, Tsrf, snow_cover, albedo_dynamic; c=..., ...)

Wrap existing state arrays as a `SnowpackDomain`. Array shapes and per-column
vector lengths are validated but the input arrays are not copied.
"""
function SnowpackDomain(
    N::AbstractVector{<:Integer},
    mass::AbstractMatrix{NF},
    mass_w::AbstractMatrix{NF},
    density::AbstractMatrix{NF},
    temperature::AbstractMatrix{NF},
    mass_base::AbstractVector{NF},
    smb_ice::AbstractVector{NF},
    runoff::AbstractVector{NF},
    Tsrf::AbstractVector{NF},
    snow_cover::AbstractVector{NF},
    albedo_dynamic::AbstractVector{NF};
    c::SnowpackPhysicalConstants{NF}=SnowpackPhysicalConstants(NF),
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    rho_max::Real=DEFAULT_RHO_MAX,
) where {NF <: AbstractFloat}
    ncol = length(N)
    size(mass, 2) == ncol || error("`mass` must have one column per entry of `N`.")
    size(mass_w) == size(mass) || error("`mass_w` must match `mass`.")
    size(density) == size(mass) || error("`density` must match `mass`.")
    size(temperature) == size(mass) || error("`temperature` must match `mass`.")
    _validate_domain_vector("mass_base", mass_base, ncol)
    _validate_domain_vector("smb_ice", smb_ice, ncol)
    _validate_domain_vector("runoff", runoff, ncol)
    _validate_domain_vector("Tsrf", Tsrf, ncol)
    _validate_domain_vector("snow_cover", snow_cover, ncol)
    _validate_domain_vector("albedo_dynamic", albedo_dynamic, ncol)
    _validate_mass_partition(mass_max, mass_split, mass_min)

    return SnowpackDomain(
        c,
        size(mass, 1),
        ncol,
        convert(NF, mass_max),
        convert(NF, mass_split),
        convert(NF, mass_min),
        convert(NF, rho_max),
        N,
        mass,
        mass_w,
        density,
        temperature,
        mass_base,
        smb_ice,
        runoff,
        Tsrf,
        snow_cover,
        albedo_dynamic,
    )
end

"""
    cpu_domain(domain)

Return a copy of `domain` adapted to CPU `Array` storage.
"""
cpu_domain(domain::SnowpackDomain) = adapt(Array, domain)

"""
    gpu_domain(domain)

Return a copy of `domain` adapted to the default GPU storage type. Throws if
CUDA is not functional in the current session.
"""
function gpu_domain(domain::SnowpackDomain, storage_type=gpu_storage_type())
    cuda_available() || error("CUDA is not functional in the current environment.")
    return adapt(storage_type, domain)
end

@adapt_structure SnowpackDomain

"""Thin model-facing view of evolving snowpack state arrays."""
struct SnowpackState{VT, MT, NI}
    mass::MT
    mass_w::MT
    density::MT
    temperature::MT
    N::NI
    surface_temperature::VT
    albedo::VT
    snow_cover::VT
    runoff::VT
    smb::VT
end

function SnowpackState(domain::SnowpackDomain)
    return SnowpackState(
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.N,
        domain.Tsrf,
        domain.albedo_dynamic,
        domain.snow_cover,
        domain.runoff,
        domain.smb_ice,
    )
end

"""
State accessors and formatted state output.
"""

"""
    _state_dict(N_storage, mass, mass_w, density, temperature, smb_ice, albedo_dynamic, idx, c)

Build a dictionary snapshot for column `idx`. The result includes active-layer
profiles, bulk totals, and surface diagnostics, and allocates new Julia arrays
for the returned profile data.
"""
function _state_dict(
    N_storage,
    mass,
    mass_w,
    density,
    temperature,
    smb_ice,
    albedo_dynamic,
    idx::Int,
    c::SnowpackPhysicalConstants,
)
    snow_cover = _snow_cover_fraction(N_storage, mass, mass_w, density, idx)
    n = _n_active(N_storage, idx)
    if n == 0
        return Dict(
            "N" => 0,
            "n_active" => 0,
            "mass" => Float64[],
            "solid_mass" => Float64[],
            "mass_w" => Float64[],
            "liquid_water_mass" => Float64[],
            "density" => Float64[],
            "total_mass" => 0.0,
            "total_liquid_water" => 0.0,
            "total_wet_mass" => 0.0,
            "thickness" => Float64[],
            "total_thickness" => 0.0,
            "surface_temperature" => c.T0,
            "snow_cover" => snow_cover,
            "surface_albedo" => _get_scalar(albedo_dynamic, idx),
            "albedo_dynamic" => _get_scalar(albedo_dynamic, idx),
            "smb_ice" => _get_scalar(smb_ice, idx),
            "ice_sheet_smb" => _get_scalar(smb_ice, idx),
        )
    end

    active_solid_mass = [@inbounds _get_layer(mass, layer_index, idx) for layer_index in 1:n]
    active_liquid_water_mass = [@inbounds _get_layer(mass_w, layer_index, idx) for layer_index in 1:n]
    active_density = [@inbounds _get_layer(density, layer_index, idx) for layer_index in 1:n]
    thickness = active_solid_mass ./ active_density

    return Dict(
        "N" => n,
        "n_active" => n,
        "mass" => active_solid_mass,
        "solid_mass" => active_solid_mass,
        "mass_w" => active_liquid_water_mass,
        "liquid_water_mass" => active_liquid_water_mass,
        "density" => active_density,
        "total_mass" => sum(active_solid_mass),
        "total_liquid_water" => sum(active_liquid_water_mass),
        "total_wet_mass" => sum(active_solid_mass) + sum(active_liquid_water_mass),
        "thickness" => thickness,
        "total_thickness" => sum(thickness),
        "surface_temperature" => _get_layer(temperature, 1, idx),
        "snow_cover" => snow_cover,
        "surface_albedo" => _get_scalar(albedo_dynamic, idx),
        "albedo_dynamic" => _get_scalar(albedo_dynamic, idx),
        "smb_ice" => _get_scalar(smb_ice, idx),
        "ice_sheet_smb" => _get_scalar(smb_ice, idx),
    )
end

"""
    get_state(domain, idx=1)

Return a dictionary snapshot for column `idx` of `domain`. Values are copied
into plain Julia containers so callers can inspect state without mutating the
domain.
"""
function get_state(domain::AbstractSnowpackDomain, idx::Int=1)
    return _state_dict(
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.temperature,
        domain.smb_ice,
        domain.albedo_dynamic,
        idx,
        domain.c,
    )
end

"""
    print_state(domain, idx=1)

Print a short formatted summary of column `idx` to standard output. This is a
diagnostic convenience wrapper around `get_state`.
"""
function print_state(domain::AbstractSnowpackDomain, idx::Int=1)
    state = get_state(domain, idx)
    println("=" ^ 60)
    println("Snowpack Domain Column State")
    println("=" ^ 60)
    println("Column index: ", idx)
    println("Active layers: ", state["N"])
    println("Total mass: ", round(state["total_mass"], digits=2), " kg/m^2")
    println("Total thickness: ", round(state["total_thickness"], digits=3), " m")
    println("Snow cover: ", round(state["snow_cover"], digits=3))
    println("Surface albedo: ", round(state["surface_albedo"], digits=3))
    println()
end

"""
    compute_auxiliary!(domain, idx)

Recompute derived diagnostics for column `idx` in-place. Currently this
updates snow-cover fraction and surface albedo.
"""
function compute_auxiliary!(domain::AbstractSnowpackDomain, idx::Int)
    update_snow_cover!(domain, idx)
    update_surface_albedo!(domain, idx)
    return nothing
end

"""
    compute_auxiliary!(domain)

Recompute derived diagnostics for every column in `domain`. Mutates the
domain’s auxiliary fields in-place and returns `nothing`.
"""
function compute_auxiliary!(domain::AbstractSnowpackDomain)
    for idx in 1:column_count(domain)
        compute_auxiliary!(domain, idx)
    end
    return nothing
end

"""
Domain-wide summary helpers.
"""

"""
    _summarize_domain_state_kernel!(...)

KernelAbstractions kernel that summarizes each column into thickness, wet
mass, bulk density, basal mass, ice SMB, liquid water, and runoff arrays. All
output arrays are mutated in-place.
"""
@kernel function _summarize_domain_state_kernel!(
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    smb_ice,
    liquid_water,
    runoff,
    N,
    mass,
    mass_w,
    density,
    mass_base_state,
    smb_ice_state,
    runoff_state,
)
    idx = @index(Global)
    if idx <= length(N)
        n = N[idx]
        thickness_local = zero(eltype(thickness))
        wet_mass_local = zero(eltype(wet_mass))
        solid_mass_local = zero(eltype(wet_mass))
        liquid_water_local = zero(eltype(liquid_water))
        for layer_index in 1:n
            solid = mass[layer_index, idx]
            liquid = mass_w[layer_index, idx]
            rho = density[layer_index, idx]
            solid_mass_local += solid
            wet_mass_local += solid + liquid
            liquid_water_local += liquid
            if solid > zero(solid) && rho > EPS_TINY
                thickness_local += solid / rho
            end
        end
        thickness[idx] = thickness_local
        wet_mass[idx] = wet_mass_local
        bulk_density[idx] = thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(bulk_density))
        base_mass[idx] = mass_base_state[idx]
        smb_ice[idx] = smb_ice_state[idx]
        liquid_water[idx] = liquid_water_local
        runoff[idx] = runoff_state[idx]
    end
end

"""
    _summarize_cycle_state_kernel!(...)

KernelAbstractions kernel that summarizes the subset of column diagnostics
needed for equilibrium-cycle tracking. Mutates the supplied output arrays
in-place.
"""
@kernel function _summarize_cycle_state_kernel!(
    thickness,
    wet_mass,
    bulk_density,
    base_mass,
    N,
    mass,
    mass_w,
    density,
    mass_base_state,
)
    idx = @index(Global)
    if idx <= length(N)
        n = N[idx]
        thickness_local = zero(eltype(thickness))
        wet_mass_local = zero(eltype(wet_mass))
        solid_mass_local = zero(eltype(wet_mass))
        for layer_index in 1:n
            solid = mass[layer_index, idx]
            liquid = mass_w[layer_index, idx]
            rho = density[layer_index, idx]
            solid_mass_local += solid
            wet_mass_local += solid + liquid
            if solid > zero(solid) && rho > EPS_TINY
                thickness_local += solid / rho
            end
        end
        thickness[idx] = thickness_local
        wet_mass[idx] = wet_mass_local
        bulk_density[idx] = thickness_local > EPS_TINY ? solid_mass_local / thickness_local : zero(eltype(bulk_density))
        base_mass[idx] = mass_base_state[idx]
    end
end

"""
    summarize_domain_state!(thickness, wet_mass, bulk_density, base_mass, smb_ice, liquid_water, runoff, domain)

Fill preallocated summary arrays with one-column diagnostics from `domain`.
Outputs are column-wise totals or aggregates in SI-like model units.
"""
function summarize_domain_state!(
    thickness::AbstractVector,
    wet_mass::AbstractVector,
    bulk_density::AbstractVector,
    base_mass::AbstractVector,
    smb_ice::AbstractVector,
    liquid_water::AbstractVector,
    runoff::AbstractVector,
    domain::AbstractSnowpackDomain,
)
    kernel! = _summarize_domain_state_kernel!(_ka_backend(domain.mass))
    event = kernel!(
        thickness,
        wet_mass,
        bulk_density,
        base_mass,
        smb_ice,
        liquid_water,
        runoff,
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.mass_base,
        domain.smb_ice,
        domain.runoff;
        ndrange=column_count(domain),
    )
    _wait_kernel(event)
    return nothing
end

"""
    summarize_domain_state(domain)

Allocate and return a named tuple of per-column summary arrays for `domain`.
This is a convenience wrapper around `summarize_domain_state!`.
"""
function summarize_domain_state(domain::AbstractSnowpackDomain)
    ncol = column_count(domain)
    NF = eltype(domain)
    allocate() = similar(domain.mass, NF, ncol)
    thickness = allocate()
    wet_mass = allocate()
    bulk_density = allocate()
    base_mass = allocate()
    smb_ice = allocate()
    liquid_water = allocate()
    runoff = allocate()
    summarize_domain_state!(
        thickness,
        wet_mass,
        bulk_density,
        base_mass,
        smb_ice,
        liquid_water,
        runoff,
        domain,
    )
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
        smb_ice=smb_ice,
        liquid_water=liquid_water,
        runoff=runoff,
    )
end

"""
    summarize_cycle_state!(thickness, wet_mass, bulk_density, base_mass, domain)

Fill preallocated arrays with the smaller summary set used to compare
equilibrium cycles. Mutates the output arrays and returns `nothing`.
"""
function summarize_cycle_state!(
    thickness::AbstractVector,
    wet_mass::AbstractVector,
    bulk_density::AbstractVector,
    base_mass::AbstractVector,
    domain::AbstractSnowpackDomain,
)
    kernel! = _summarize_cycle_state_kernel!(_ka_backend(domain.mass))
    event = kernel!(
        thickness,
        wet_mass,
        bulk_density,
        base_mass,
        domain.N,
        domain.mass,
        domain.mass_w,
        domain.density,
        domain.mass_base;
        ndrange=column_count(domain),
    )
    _wait_kernel(event)
    return nothing
end

"""
    summarize_cycle_state(domain)

Allocate and return a named tuple of cycle-level summary arrays for `domain`.
This is the allocating counterpart to `summarize_cycle_state!`.
"""
function summarize_cycle_state(domain::AbstractSnowpackDomain)
    ncol = column_count(domain)
    NF = eltype(domain)
    allocate() = similar(domain.mass, NF, ncol)
    thickness = allocate()
    wet_mass = allocate()
    bulk_density = allocate()
    base_mass = allocate()
    summarize_cycle_state!(thickness, wet_mass, bulk_density, base_mass, domain)
    return (
        thickness=thickness,
        wet_mass=wet_mass,
        bulk_density=bulk_density,
        base_mass=base_mass,
    )
end
