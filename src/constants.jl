"""
Constants used across the snowpack model: numeric tolerances, default configuration
values, scheme flags, and physical constants.
"""

# ---------------------------------------------------------------------------
# Numeric tolerances
# ---------------------------------------------------------------------------
const EPS_TINY = 1.0e-12
const EPS_EMPTY_LAYER = 1.0e-10

# ---------------------------------------------------------------------------
# Time conversion defaults
# ---------------------------------------------------------------------------
const DEFAULT_SECONDS_PER_DAY = 86_400.0
const DEFAULT_SECONDS_PER_MONTH = DEFAULT_SECONDS_PER_DAY * 30.0
const DEFAULT_SECONDS_PER_YEAR = DEFAULT_SECONDS_PER_MONTH * 12.0

# ---------------------------------------------------------------------------
# Snowpack column defaults
# ---------------------------------------------------------------------------
const DEFAULT_NTOT = 2
const DEFAULT_N_ACTIVE = 0
const DEFAULT_MASS_MAX = 500.0
const DEFAULT_MASS_SPLIT = 300.0
const DEFAULT_MASS_MIN = 100.0
const DEFAULT_RHO_MAX = 900.0
const DEFAULT_DENSITY_INIT = 300.0
const DEFAULT_TEMPERATURE_INIT = 273.0

# Legacy BESSI reference uses a 15-layer column when capping total column depth.
const BESSI_REFERENCE_LAYER_COUNT = 15
const BESSI_REFERENCE_DEPTH_DENSITY = 300.0

# ---------------------------------------------------------------------------
# Scheme flag constants
# ---------------------------------------------------------------------------
const FRESH_SNOW_DENSITY_CONSTANT = UInt8(1)
const FRESH_SNOW_DENSITY_PARAMETERIZED = UInt8(2)
const ALBEDO_CONSTANT = UInt8(1)
const ALBEDO_DYNAMIC = UInt8(2)
const ALBEDO_PRESCRIBED = UInt8(3)
const LOW_DENSIFICATION_BESSI = UInt8(1)
const LOW_DENSIFICATION_HTESSEL = UInt8(2)

# ---------------------------------------------------------------------------
# Physical constants struct and constructors
# ---------------------------------------------------------------------------

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
    Lv::NF
    cp_air::NF
    latent_heat_flux_ratio::NF
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

Base.eltype(::SnowpackPhysicalConstants{NF}) where {NF} = NF

"""
    number_type(c)

Return the floating-point element type used by the physical constants set `c`.
"""
@inline number_type(::SnowpackPhysicalConstants{NF}) where {NF} = NF

@inline _uses_constant_fresh_snow_density(c::SnowpackPhysicalConstants) =
    c.fresh_snow_density_scheme == FRESH_SNOW_DENSITY_CONSTANT

@inline _uses_constant_albedo(c::SnowpackPhysicalConstants) =
    c.albedo_scheme == ALBEDO_CONSTANT

@inline _uses_prescribed_albedo(c::SnowpackPhysicalConstants) =
    c.albedo_scheme == ALBEDO_PRESCRIBED

@inline _uses_htessel_densification(c::SnowpackPhysicalConstants) =
    c.low_density_densification == LOW_DENSIFICATION_HTESSEL

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
    normalized_scheme in (:constant, :dynamic, :prescribed) ||
        error(
            "Unsupported albedo scheme '$scheme'. " *
            "Use :constant, :dynamic, :prescribed, or the aliases :legacy / :bessi.",
        )
    return normalized_scheme == :constant ? ALBEDO_CONSTANT :
        normalized_scheme == :prescribed ? ALBEDO_PRESCRIBED :
        ALBEDO_DYNAMIC
end

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
    Lv::Real=2.501e6,
    cp_air::Real=1003.0,
    latent_heat_flux_ratio::Real=1.5,
    D_sh::Real=20.0,
    alpha_dry::Real=0.81,
    alpha_wet::Real=0.70,
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
        convert(NF, Lv),
        convert(NF, cp_air),
        convert(NF, latent_heat_flux_ratio),
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

@adapt_structure SnowpackPhysicalConstants
