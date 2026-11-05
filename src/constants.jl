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

# ---------------------------------------------------------------------------
# Atmospheric defaults
# ---------------------------------------------------------------------------
const DEFAULT_SEA_LEVEL_AIR_PRESSURE = 101_325.0
const DEFAULT_GRAVITY = 9.80665
const DEFAULT_MOLAR_MASS_DRY_AIR = 0.0289644
const DEFAULT_UNIVERSAL_GAS_CONSTANT = 8.31446261815324

# ---------------------------------------------------------------------------
# Snowpack column defaults
# ---------------------------------------------------------------------------
const DEFAULT_NTOT = 15
const DEFAULT_N_ACTIVE = 0
const DEFAULT_MASS_MAX = 500.0
const DEFAULT_MASS_SPLIT = 300.0
const DEFAULT_MASS_MIN = 100.0
const DEFAULT_DENSITY_INIT = 300.0
const DEFAULT_TEMPERATURE_INIT = 273.0

# Legacy BESSI reference uses a 15-layer column when capping total column depth.
const BESSI_REFERENCE_LAYER_COUNT = 15
const BESSI_REFERENCE_DEPTH_DENSITY = 300.0
const BESSI_REFERENCE_SNOW_DEPTH_M = 22.5

# ---------------------------------------------------------------------------
# Scheme types
# ---------------------------------------------------------------------------
# Process schemes are singleton types carried as type parameters of
# `SnowpackPhysicalConstants`; process code dispatches on them, so a new scheme
# is a new type plus its methods.

abstract type AbstractFreshSnowDensityScheme end
"""Fresh snow takes the constant density `rho_s`."""
struct ConstantFreshSnowDensity <: AbstractFreshSnowDensityScheme end
"""Fresh-snow density from air temperature and wind speed."""
struct ParameterizedFreshSnowDensity <: AbstractFreshSnowDensityScheme end

abstract type AbstractAlbedoScheme end
"""Fixed dry/wet snow albedo depending on the surface temperature."""
struct ConstantAlbedo <: AbstractAlbedoScheme end
"""BESSI snow aging and wetness relaxation."""
struct DynamicAlbedo <: AbstractAlbedoScheme end
"""Albedo from the `prescribed_albedo` forcing (dynamic where it is absent)."""
struct PrescribedAlbedo <: AbstractAlbedoScheme end
"""Exponential decay with the time since the last snowfall."""
struct AgingAlbedo <: AbstractAlbedoScheme end
"""Two-band SEMIX snow albedo."""
struct SemixAlbedo <: AbstractAlbedoScheme end

abstract type AbstractDensificationScheme end
"""BESSI low-density densification."""
struct BESSIDensification <: AbstractDensificationScheme end
"""HTESSEL low-density densification with liquid-water compaction."""
struct HTESSELDensification <: AbstractDensificationScheme end

const _FRESH_SNOW_DENSITY_SCHEMES = (
    constant=ConstantFreshSnowDensity(),
    parameterized=ParameterizedFreshSnowDensity(),
)
const _ALBEDO_SCHEMES = (
    constant=ConstantAlbedo(),
    dynamic=DynamicAlbedo(),
    prescribed=PrescribedAlbedo(),
    aging=AgingAlbedo(),
    semix=SemixAlbedo(),
)
const _DENSIFICATION_SCHEMES = (bessi=BESSIDensification(), htessel=HTESSELDensification())

# ---------------------------------------------------------------------------
# Runtime scheme flags
# ---------------------------------------------------------------------------
const SEB_BESSI = UInt8(1)
const SEB_SEMIX = UInt8(2)
const LONGWAVE_GRAYBODY = UInt8(1)
const LONGWAVE_CLOUD_PROXY = UInt8(2)
const SEMIX_ALBEDO_WW = UInt8(1)
const SEMIX_ALBEDO_DANG = UInt8(2)

# ---------------------------------------------------------------------------
# Physical constants struct and constructors
# ---------------------------------------------------------------------------

"""
    SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification}

Container for physical constants, empirical coefficients, and scheme flags
used by the snowpack model. The fresh-snow density, albedo and low-density
densification schemes are type parameters (see `AbstractAlbedoScheme` etc.).
"""
struct SnowpackPhysicalConstants{
        NF <: AbstractFloat,
        FreshSnowDensity <: AbstractFreshSnowDensityScheme,
        Albedo <: AbstractAlbedoScheme,
        Densification <: AbstractDensificationScheme,
    }
    rho_s::NF
    rho_i::NF
    rho_w::NF
    rho_s_a::NF
    rho_s_b::NF
    rho_s_c::NF
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
    aging_cold_timescale_days::NF
    aging_melting_timescale_days::NF
    ϵ_air::NF
    longwave_scheme::UInt8
    lw_clear_sky_transmissivity::NF
    lw_clear_sky_transmissivity_per_km::NF
    lw_clear_sky_transmissivity_insolation::NF
    lw_cloud_transmission_coefficient::NF
    lw_cloud_transmission_exponent::NF
    lw_cloud_transmission_decay_per_km::NF
    lw_clear_sky_emissivity::NF
    lw_clear_sky_emissivity_temperature_slope::NF
    lw_overcast_emissivity::NF
    lw_overcast_emissivity_temperature_slope::NF
    lw_cloud_emissivity_exponent::NF
    lw_night_cloud_fraction::NF
    lw_emissivity_max::NF
    ϵ_snow::NF
    σ::NF
    T0::NF
    seconds_per_day::NF
    seb_scheme::UInt8
    turbulent_flux_scheme::UInt8
    eps_ice::NF
    semix_karman::NF
    semix_surface_height::NF
    semix_z0m_snow::NF
    semix_z0m_ice::NF
    semix_zm_to_zh::NF
    semix_sensible_exchange_factor::NF
    semix_stable_coefficient::NF
    semix_latent_exchange_factor::NF
    semix_snow_albedo::UInt8
    semix_frac_vu::NF
    semix_alb_snow_vis_new::NF
    semix_alb_snow_nir_new::NF
    semix_snow_grain_fresh::NF
    semix_snow_grain_old::NF
    semix_d_alb_age_vis::NF
    semix_d_alb_age_nir::NF
    semix_f_age_t::NF
    semix_dT_age::NF
    semix_snow_0::NF
    semix_snow_1::NF
    semix_w_snow_dust::NF
    semix_dust_con_scale::NF
    semix_dalb_snow_vis::NF
    semix_dalb_snow_nir::NF
    semix_k_sigma_orog::NF
    semix_sigma_orog_crit::NF
end

"""Fresh-snow density scheme of `c` (an `AbstractFreshSnowDensityScheme` instance)."""
@inline fresh_snow_density_scheme(::SnowpackPhysicalConstants{NF, F}) where {NF, F} = F()
"""Albedo scheme of `c` (an `AbstractAlbedoScheme` instance)."""
@inline albedo_scheme(::SnowpackPhysicalConstants{NF, F, A}) where {NF, F, A} = A()
"""Low-density densification scheme of `c` (an `AbstractDensificationScheme` instance)."""
@inline densification_scheme(::SnowpackPhysicalConstants{NF, F, A, D}) where {NF, F, A, D} = D()

@inline function Base.getproperty(c::SnowpackPhysicalConstants, name::Symbol)
    name === :fresh_snow_density_scheme && return fresh_snow_density_scheme(c)
    name === :albedo_scheme && return albedo_scheme(c)
    name === :low_density_densification && return densification_scheme(c)
    return getfield(c, name)
end

Base.propertynames(c::SnowpackPhysicalConstants, private::Bool=false) =
    (fieldnames(typeof(c))..., :fresh_snow_density_scheme, :albedo_scheme, :low_density_densification)

Base.eltype(::SnowpackPhysicalConstants{NF}) where {NF} = NF

"""
    number_type(c)

Return the floating-point element type used by the physical constants set `c`.
"""
@inline number_type(::SnowpackPhysicalConstants{NF}) where {NF} = NF

@inline _uses_semix_seb(c::SnowpackPhysicalConstants) = c.seb_scheme == SEB_SEMIX

"""Whether turbulent sensible and latent heat use the SEMIX bulk formulation."""
@inline _uses_semix_turbulence(c::SnowpackPhysicalConstants) =
    c.turbulent_flux_scheme == SEB_SEMIX

@inline _initial_snow_albedo(c::SnowpackPhysicalConstants) =
    c.alpha_dry

"""
    _normalize_low_density_densification(scheme)

Validate a low-density densification scheme (`:bessi`, `:htessel` or an
`AbstractDensificationScheme`) and return the scheme instance.
"""
@inline _normalize_low_density_densification(scheme::AbstractDensificationScheme) = scheme
@inline function _normalize_low_density_densification(scheme::Symbol)
    haskey(_DENSIFICATION_SCHEMES, scheme) ||
        error("Unsupported low-density densification scheme '$scheme'. Use :bessi or :htessel.")
    return _DENSIFICATION_SCHEMES[scheme]
end

"""
    _normalize_fresh_snow_density_scheme(scheme)

Validate a fresh-snow density scheme (`:constant`, `:parameterized` or an
`AbstractFreshSnowDensityScheme`) and return the scheme instance.
"""
@inline _normalize_fresh_snow_density_scheme(scheme::AbstractFreshSnowDensityScheme) = scheme
@inline function _normalize_fresh_snow_density_scheme(scheme::Symbol)
    haskey(_FRESH_SNOW_DENSITY_SCHEMES, scheme) ||
        error(
            "Unsupported fresh-snow density scheme '$scheme'. " *
            "Use :constant or :parameterized.",
        )
    return _FRESH_SNOW_DENSITY_SCHEMES[scheme]
end

"""
    _normalize_albedo_scheme(scheme)

Validate an albedo scheme (`:constant`, `:dynamic`, `:prescribed`, `:aging`,
`:semix` or an `AbstractAlbedoScheme`) and return the scheme instance.
"""
@inline _normalize_albedo_scheme(scheme::AbstractAlbedoScheme) = scheme
@inline function _normalize_albedo_scheme(scheme::Symbol)
    haskey(_ALBEDO_SCHEMES, scheme) ||
        error(
            "Unsupported albedo scheme '$scheme'. " *
            "Use :constant, :dynamic, :prescribed, :aging, or :semix.",
        )
    return _ALBEDO_SCHEMES[scheme]
end

@inline function _normalize_longwave_scheme(scheme::Symbol)
    scheme in (:graybody, :cloud_proxy) ||
        error("Unsupported longwave scheme '$scheme'. Use :graybody or :cloud_proxy.")
    return scheme == :cloud_proxy ? LONGWAVE_CLOUD_PROXY : LONGWAVE_GRAYBODY
end

@inline function _normalize_seb_scheme(scheme::Symbol)
    scheme in (:bessi, :semix) ||
        error("Unsupported surface-energy-balance scheme '$scheme'. Use :bessi or :semix.")
    return scheme == :semix ? SEB_SEMIX : SEB_BESSI
end

@inline function _normalize_semix_snow_albedo(scheme::Symbol)
    scheme in (:warren_wiscombe, :dang) ||
        error("Unsupported `semix_snow_albedo` '$scheme'. Use :warren_wiscombe or :dang.")
    return scheme == :dang ? SEMIX_ALBEDO_DANG : SEMIX_ALBEDO_WW
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
    fresh_snow_density_scheme::Union{Symbol, AbstractFreshSnowDensityScheme}=:constant,
    Ki::Real=2.1,
    ci::Real=2110.0,
    cw::Real=4181.0,
    Lm::Real=334000.0,
    Lv::Real=2.501e6,
    cp_air::Real=1003.0,
    latent_heat_flux_ratio::Real=1,
    D_sh::Real=10.0,
    alpha_dry::Real=0.81,
    alpha_wet::Real=0.70,
    alpha_ice::Real=0.40,
    max_lwc_albedo::Real=0.1,
    aging_cold_timescale_days::Real=20.0,
    aging_melting_timescale_days::Real=2.0,
    albedo_scheme::Union{Symbol, AbstractAlbedoScheme}=:dynamic,
    ϵ_air::Real=0.8,
    longwave_scheme::Symbol=:cloud_proxy,
    lw_clear_sky_transmissivity::Real=0.718,
    lw_clear_sky_transmissivity_per_km::Real=0.0298,
    lw_clear_sky_transmissivity_insolation::Real=0.0569,
    lw_cloud_transmission_coefficient::Real=0.447,
    lw_cloud_transmission_exponent::Real=2.75,
    lw_cloud_transmission_decay_per_km::Real=0.225,
    lw_clear_sky_emissivity::Real=0.716,
    lw_clear_sky_emissivity_temperature_slope::Real=0.00239,
    lw_overcast_emissivity::Real=0.952,
    lw_overcast_emissivity_temperature_slope::Real=0.00362,
    lw_cloud_emissivity_exponent::Real=3.0,
    lw_night_cloud_fraction::Real=0.85,
    lw_emissivity_max::Real=1.0,
    ϵ_snow::Real=0.98,
    σ::Real=5.670373e-8,
    T0::Real=273.15,
    seconds_per_day::Real=DEFAULT_SECONDS_PER_DAY,
    low_density_densification::Union{Symbol, AbstractDensificationScheme}=:bessi,
    seb_scheme::Symbol=:semix,
    turbulent_flux_scheme::Symbol=:semix,
    eps_ice::Real=0.98,
    semix_karman::Real=0.4,
    semix_surface_height::Real=10.0,
    semix_z0m_snow::Real=0.001,
    semix_z0m_ice::Real=0.01,
    semix_zm_to_zh::Real=10.0,
    semix_sensible_exchange_factor::Real=3.0,
    semix_stable_coefficient::Real=40.0,
    semix_latent_exchange_factor::Real=1.0,
    semix_snow_albedo::Symbol=:dang,
    semix_frac_vu::Real=0.45,
    semix_alb_snow_vis_new::Real=0.99,
    semix_alb_snow_nir_new::Real=0.65,
    semix_snow_grain_fresh::Real=50.0,
    semix_snow_grain_old::Real=1000.0,
    semix_d_alb_age_vis::Real=0.05,
    semix_d_alb_age_nir::Real=0.25,
    semix_f_age_t::Real=0.1,
    semix_dT_age::Real=0.0,
    semix_snow_0::Real=1.0,
    semix_snow_1::Real=0.5,
    semix_w_snow_dust::Real=10.0,
    semix_dust_con_scale::Real=1.0,
    semix_dalb_snow_vis::Real=0.0,
    semix_dalb_snow_nir::Real=0.0,
    semix_k_sigma_orog::Real=0.0,
    semix_sigma_orog_crit::Real=1000.0,
) where {NF <: AbstractFloat}
    resolved_albedo_scheme = _normalize_albedo_scheme(albedo_scheme)
    if resolved_albedo_scheme isa AgingAlbedo
        0 <= alpha_wet <= alpha_dry <= 1 || error(
            "Aging albedos must satisfy 0 <= alpha_wet <= alpha_dry <= 1.",
        )
    end
    aging_cold_timescale_days > 0 || error("`aging_cold_timescale_days` must be positive.")
    aging_melting_timescale_days > 0 || error("`aging_melting_timescale_days` must be positive.")
    resolved_fresh_snow_density_scheme = _normalize_fresh_snow_density_scheme(fresh_snow_density_scheme)
    resolved_densification_scheme = _normalize_low_density_densification(low_density_densification)
    seb_flag = _normalize_seb_scheme(seb_scheme)
    turbulent_flux_flag = _normalize_seb_scheme(turbulent_flux_scheme)
    semix_albedo_flag = _normalize_semix_snow_albedo(semix_snow_albedo)
    longwave_flag = _normalize_longwave_scheme(longwave_scheme)
    lw_clear_sky_transmissivity > 0 || error("`lw_clear_sky_transmissivity` must be positive.")
    lw_cloud_transmission_coefficient > 0 || error("`lw_cloud_transmission_coefficient` must be positive.")
    lw_cloud_transmission_exponent > 0 || error("`lw_cloud_transmission_exponent` must be positive.")
    lw_cloud_emissivity_exponent > 0 || error("`lw_cloud_emissivity_exponent` must be positive.")
    0 <= lw_night_cloud_fraction <= 1 || error("`lw_night_cloud_fraction` must be in [0, 1].")
    0.4 < lw_emissivity_max <= 1 || error("`lw_emissivity_max` must be in (0.4, 1].")
    semix_karman > 0 || error("`semix_karman` must be positive.")
    semix_surface_height > 0 || error("`semix_surface_height` must be positive.")
    semix_z0m_snow > 0 || error("`semix_z0m_snow` must be positive.")
    semix_z0m_ice > 0 || error("`semix_z0m_ice` must be positive.")
    semix_zm_to_zh > 0 || error("`semix_zm_to_zh` must be positive.")
    semix_sensible_exchange_factor > 0 || error("`semix_sensible_exchange_factor` must be positive.")
    semix_stable_coefficient >= 0 || error("`semix_stable_coefficient` must be nonnegative.")
    semix_latent_exchange_factor > 0 || error("`semix_latent_exchange_factor` must be positive.")
    return SnowpackPhysicalConstants{
        NF,
        typeof(resolved_fresh_snow_density_scheme),
        typeof(resolved_albedo_scheme),
        typeof(resolved_densification_scheme),
    }(
        convert(NF, rho_s),
        convert(NF, rho_i),
        convert(NF, rho_w),
        convert(NF, rho_s_a),
        convert(NF, rho_s_b),
        convert(NF, rho_s_c),
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
        convert(NF, aging_cold_timescale_days),
        convert(NF, aging_melting_timescale_days),
        convert(NF, ϵ_air),
        longwave_flag,
        convert(NF, lw_clear_sky_transmissivity),
        convert(NF, lw_clear_sky_transmissivity_per_km),
        convert(NF, lw_clear_sky_transmissivity_insolation),
        convert(NF, lw_cloud_transmission_coefficient),
        convert(NF, lw_cloud_transmission_exponent),
        convert(NF, lw_cloud_transmission_decay_per_km),
        convert(NF, lw_clear_sky_emissivity),
        convert(NF, lw_clear_sky_emissivity_temperature_slope),
        convert(NF, lw_overcast_emissivity),
        convert(NF, lw_overcast_emissivity_temperature_slope),
        convert(NF, lw_cloud_emissivity_exponent),
        convert(NF, lw_night_cloud_fraction),
        convert(NF, lw_emissivity_max),
        convert(NF, ϵ_snow),
        convert(NF, σ),
        convert(NF, T0),
        convert(NF, seconds_per_day),
        seb_flag,
        turbulent_flux_flag,
        convert(NF, eps_ice),
        convert(NF, semix_karman),
        convert(NF, semix_surface_height),
        convert(NF, semix_z0m_snow),
        convert(NF, semix_z0m_ice),
        convert(NF, semix_zm_to_zh),
        convert(NF, semix_sensible_exchange_factor),
        convert(NF, semix_stable_coefficient),
        convert(NF, semix_latent_exchange_factor),
        semix_albedo_flag,
        convert(NF, semix_frac_vu),
        convert(NF, semix_alb_snow_vis_new),
        convert(NF, semix_alb_snow_nir_new),
        convert(NF, semix_snow_grain_fresh),
        convert(NF, semix_snow_grain_old),
        convert(NF, semix_d_alb_age_vis),
        convert(NF, semix_d_alb_age_nir),
        convert(NF, semix_f_age_t),
        convert(NF, semix_dT_age),
        convert(NF, semix_snow_0),
        convert(NF, semix_snow_1),
        convert(NF, semix_w_snow_dust),
        convert(NF, semix_dust_con_scale),
        convert(NF, semix_dalb_snow_vis),
        convert(NF, semix_dalb_snow_nir),
        convert(NF, semix_k_sigma_orog),
        convert(NF, semix_sigma_orog_crit),
    )
end

"""
    SnowpackPhysicalConstants(; kwargs...)

Convenience constructor for `SnowpackPhysicalConstants{Float64}`.
"""
SnowpackPhysicalConstants(; kwargs...) = SnowpackPhysicalConstants(Float64; kwargs...)

# Only scalars and flags: the constants pass to GPU kernels unchanged. (A generic
# `@adapt_structure` would rebuild them without their scheme type parameters.)
Adapt.adapt_structure(to, c::SnowpackPhysicalConstants) = c
