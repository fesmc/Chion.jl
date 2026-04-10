"""
Core snowpack constants and domain state containers.
"""

const FRESH_SNOW_DENSITY_CONSTANT = UInt8(1)
const FRESH_SNOW_DENSITY_PARAMETERIZED = UInt8(2)
const ALBEDO_CONSTANT = UInt8(1)
const ALBEDO_DYNAMIC = UInt8(2)
const LOW_DENSIFICATION_BESSI = UInt8(1)
const LOW_DENSIFICATION_HTESSEL = UInt8(2)

"""
    SnowpackPhysicalConstants{NF}

Container for physical constants, empirical coefficients, and scheme flags
used by the snowpack model.
"""
struct SnowpackPhysicalConstants{NF <: AbstractFloat}
    rho_s::NF
    rho_i::NF
    rho_w::NF
    rho_s_a::NF
    rho_s_b::NF
    rho_s_c::NF
    fresh_snow_density_scheme::UInt8
    Ki::NF
    ci::NF
    cw::NF
    Lm::NF
    D_sh::NF
    alpha_dry::NF
    alpha_wet::NF
    alpha_ice::NF
    max_lwc_albedo::NF
    albedo_scheme::UInt8
    ϵ_air::NF
    ϵ_snow::NF
    σ::NF
    R::NF
    T0::NF
    seconds_per_day::NF
    seconds_per_month::NF
    seconds_per_year::NF
    low_density_densification::UInt8
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

Base.eltype(::SnowpackPhysicalConstants{NF}) where {NF} = NF
Base.eltype(::AbstractSnowpackDomain{NF}) where {NF} = NF

"""
    number_type(c)

Return the floating-point element type used by the physical constants set `c`.
"""
@inline number_type(::SnowpackPhysicalConstants{NF}) where {NF} = NF

"""
    column_count(domain)

Return the number of snow columns stored in `domain`.
"""
@inline column_count(domain::AbstractSnowpackDomain) = domain.ncol

"""
    _normalize_low_density_densification(scheme)

Normalize a densification-scheme symbol into the internal UInt8 flag used by
`SnowpackPhysicalConstants`.
"""
@inline function _normalize_low_density_densification(scheme::Symbol)
    scheme in (:bessi, :htessel) ||
        error("Unsupported low-density densification scheme '$scheme'. Use :bessi or :htessel.")
    return scheme == :htessel ? LOW_DENSIFICATION_HTESSEL : LOW_DENSIFICATION_BESSI
end

"""
    _normalize_fresh_snow_density_scheme(scheme)

Normalize a fresh-snow density scheme symbol into the internal UInt8 flag,
including support for legacy aliases.
"""
@inline function _normalize_fresh_snow_density_scheme(scheme::Symbol)
    normalized_scheme = if scheme == :bessi
        :constant
    elseif scheme == :htessel
        :parameterized
    else
        scheme
    end
    normalized_scheme in (:constant, :parameterized) ||
        error(
            "Unsupported fresh-snow density scheme '$scheme'. " *
            "Use :constant, :parameterized, or the aliases :bessi / :htessel.",
        )
    return normalized_scheme == :constant ? FRESH_SNOW_DENSITY_CONSTANT : FRESH_SNOW_DENSITY_PARAMETERIZED
end

"""
    _normalize_albedo_scheme(scheme)

Normalize an albedo scheme symbol into the internal UInt8 flag, including
legacy aliases.
"""
@inline function _normalize_albedo_scheme(scheme::Symbol)
    normalized_scheme = if scheme in (:bessi, :legacy)
        :constant
    else
        scheme
    end
    normalized_scheme in (:constant, :dynamic) ||
        error(
            "Unsupported albedo scheme '$scheme'. " *
            "Use :constant, :dynamic, or the aliases :legacy / :bessi.",
        )
    return normalized_scheme == :constant ? ALBEDO_CONSTANT : ALBEDO_DYNAMIC
end

@inline _uses_constant_fresh_snow_density(c::SnowpackPhysicalConstants) =
    c.fresh_snow_density_scheme == FRESH_SNOW_DENSITY_CONSTANT

@inline _uses_constant_albedo(c::SnowpackPhysicalConstants) =
    c.albedo_scheme == ALBEDO_CONSTANT

@inline _uses_htessel_densification(c::SnowpackPhysicalConstants) =
    c.low_density_densification == LOW_DENSIFICATION_HTESSEL

"""
    SnowpackPhysicalConstants(::Type{NF}; kwargs...)

Construct a self-consistent set of physical constants and scheme flags using
floating-point type `NF`.
"""
function SnowpackPhysicalConstants(::Type{NF};
    rho_s::Real=315.0,
    rho_i::Real=917.0,
    rho_w::Real=1000.0,
    rho_s_a::Real=109.0,
    rho_s_b::Real=6.0,
    rho_s_c::Real=26.0,
    fresh_snow_density_scheme::Symbol=:constant,
    Ki::Real=2.1,
    ci::Real=2110.0,
    cw::Real=4181.0,
    Lm::Real=334000.0,
    D_sh::Real=10.0,
    alpha_dry::Real=0.85,
    alpha_wet::Real=0.72,
    alpha_ice::Real=0.3,
    max_lwc_albedo::Real=0.1,
    albedo_scheme::Symbol=:dynamic,
    ϵ_air::Real=0.75,
    ϵ_snow::Real=0.98,
    σ::Real=5.670373e-8,
    R::Real=8.314,
    T0::Real=273.15,
    seconds_per_day::Real=DEFAULT_SECONDS_PER_DAY,
    seconds_per_month::Real=DEFAULT_SECONDS_PER_MONTH,
    seconds_per_year::Real=DEFAULT_SECONDS_PER_YEAR,
    low_density_densification::Symbol=:bessi,
) where {NF <: AbstractFloat}
    return SnowpackPhysicalConstants(
        convert(NF, rho_s),
        convert(NF, rho_i),
        convert(NF, rho_w),
        convert(NF, rho_s_a),
        convert(NF, rho_s_b),
        convert(NF, rho_s_c),
        _normalize_fresh_snow_density_scheme(fresh_snow_density_scheme),
        convert(NF, Ki),
        convert(NF, ci),
        convert(NF, cw),
        convert(NF, Lm),
        convert(NF, D_sh),
        convert(NF, alpha_dry),
        convert(NF, alpha_wet),
        convert(NF, alpha_ice),
        convert(NF, max_lwc_albedo),
        _normalize_albedo_scheme(albedo_scheme),
        convert(NF, ϵ_air),
        convert(NF, ϵ_snow),
        convert(NF, σ),
        convert(NF, R),
        convert(NF, T0),
        convert(NF, seconds_per_day),
        convert(NF, seconds_per_month),
        convert(NF, seconds_per_year),
        _normalize_low_density_densification(low_density_densification),
    )
end

"""
    SnowpackPhysicalConstants(; kwargs...)

Convenience constructor for `SnowpackPhysicalConstants{Float64}`.
"""
SnowpackPhysicalConstants(; kwargs...) = SnowpackPhysicalConstants(Float64; kwargs...)

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

Return a copy of `domain` adapted to `CUDA.CuArray` storage. Throws if CUDA is
not functional in the current session.
"""
function gpu_domain(domain::SnowpackDomain)
    cuda_available() || error("CUDA is not functional in the current environment.")
    return adapt(CUDA.CuArray, domain)
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

@adapt_structure SnowpackPhysicalConstants
@adapt_structure SnowpackDomain
