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

"""
    variables(domain)

Return metadata describing the published prognostic and auxiliary fields of a
snowpack domain.
"""
variables(::AbstractSnowpackDomain) = (
    PrognosticVariable{:mass}("Solid snow and firn mass per layer."),
    PrognosticVariable{:mass_w}("Liquid water mass per layer."),
    PrognosticVariable{:density}("Bulk snow density per layer."),
    PrognosticVariable{:temperature}("Layer temperature."),
    AuxiliaryVariable{:snow_cover}("Diagnosed snow-cover fraction."),
    AuxiliaryVariable{:albedo_dynamic}("Surface albedo used for radiative forcing."),
)

@adapt_structure SnowpackDomain
