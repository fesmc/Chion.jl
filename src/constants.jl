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

# ---------------------------------------------------------------------------
# Scheme flag constants
# ---------------------------------------------------------------------------
const FRESH_SNOW_DENSITY_CONSTANT = UInt8(1)
const FRESH_SNOW_DENSITY_PARAMETERIZED = UInt8(2)
const ALBEDO_CONSTANT = UInt8(1)
const ALBEDO_DYNAMIC = UInt8(2)
const ALBEDO_PRESCRIBED = UInt8(3)
const ALBEDO_AGING = UInt8(4)
const ALBEDO_SEMIX = UInt8(5)
const LOW_DENSIFICATION_BESSI = UInt8(1)
const LOW_DENSIFICATION_HTESSEL = UInt8(2)
const SEB_BESSI = UInt8(1)
const SEB_SEMIX = UInt8(2)
const SEMIX_ALBEDO_WW = UInt8(1)
const SEMIX_ALBEDO_DANG = UInt8(2)

# ---------------------------------------------------------------------------
# Physical constants struct and constructors
# ---------------------------------------------------------------------------

"""
    SnowpackPhysicalConstants{NF}

Container for physical constants, empirical coefficients, and scheme flags
used by the snowpack model.
"""
struct SnowpackPhysicalConstants{
        NF <: AbstractFloat,
        FreshSnowDensity,
        Albedo,
        Densification,
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

@inline _fresh_snow_density_flag(
    ::SnowpackPhysicalConstants{NF, Scheme, Albedo, Densification},
) where {NF, Scheme, Albedo, Densification} = Scheme === :constant ? FRESH_SNOW_DENSITY_CONSTANT : FRESH_SNOW_DENSITY_PARAMETERIZED

@inline _albedo_flag(
    ::SnowpackPhysicalConstants{NF, Fresh, Scheme, Densification},
) where {NF, Fresh, Scheme, Densification} = Scheme === :constant ? ALBEDO_CONSTANT :
    Scheme === :prescribed ? ALBEDO_PRESCRIBED :
    Scheme === :aging ? ALBEDO_AGING :
    Scheme === :semix ? ALBEDO_SEMIX : ALBEDO_DYNAMIC

@inline _densification_flag(
    ::SnowpackPhysicalConstants{NF, Fresh, Albedo, Scheme},
) where {NF, Fresh, Albedo, Scheme} = Scheme === :htessel ? LOW_DENSIFICATION_HTESSEL : LOW_DENSIFICATION_BESSI

@inline function Base.getproperty(c::SnowpackPhysicalConstants, name::Symbol)
    name === :fresh_snow_density_scheme && return _fresh_snow_density_flag(c)
    name === :albedo_scheme && return _albedo_flag(c)
    name === :low_density_densification && return _densification_flag(c)
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

@inline _uses_constant_fresh_snow_density(
    ::SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification},
) where {NF, FreshSnowDensity, Albedo, Densification} = FreshSnowDensity === :constant

@inline _uses_constant_albedo(
    ::SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification},
) where {NF, FreshSnowDensity, Albedo, Densification} = Albedo === :constant

@inline _uses_prescribed_albedo(
    ::SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification},
) where {NF, FreshSnowDensity, Albedo, Densification} = Albedo === :prescribed

@inline _uses_aging_albedo(
    ::SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification},
) where {NF, FreshSnowDensity, Albedo, Densification} = Albedo === :aging

@inline _uses_semix_albedo(
    ::SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification},
) where {NF, FreshSnowDensity, Albedo, Densification} = Albedo === :semix

@inline _uses_semix_seb(c::SnowpackPhysicalConstants) = c.seb_scheme == SEB_SEMIX

"""Whether turbulent sensible and latent heat use the SEMIX bulk formulation."""
@inline _uses_semix_turbulence(c::SnowpackPhysicalConstants) =
    c.turbulent_flux_scheme == SEB_SEMIX

@inline _initial_snow_albedo(c::SnowpackPhysicalConstants) =
    c.alpha_dry

@inline _uses_htessel_densification(
    ::SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification},
) where {NF, FreshSnowDensity, Albedo, Densification} = Densification === :htessel

@inline _fresh_snow_density_tag(::SnowpackPhysicalConstants{NF, Scheme, Albedo, Densification}) where {NF, Scheme, Albedo, Densification} = Val(Scheme)
@inline _albedo_tag(::SnowpackPhysicalConstants{NF, Fresh, Scheme, Densification}) where {NF, Fresh, Scheme, Densification} = Val(Scheme)
@inline _densification_tag(::SnowpackPhysicalConstants{NF, Fresh, Albedo, Scheme}) where {NF, Fresh, Albedo, Scheme} = Val(Scheme)

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
    normalized_scheme in (:constant, :dynamic, :prescribed, :aging, :semix) ||
        error(
            "Unsupported albedo scheme '$scheme'. " *
            "Use :constant, :dynamic, :prescribed, :aging, :semix, or the aliases :legacy / :bessi.",
        )
    return normalized_scheme == :constant ? ALBEDO_CONSTANT :
        normalized_scheme == :prescribed ? ALBEDO_PRESCRIBED :
        normalized_scheme == :aging ? ALBEDO_AGING :
        normalized_scheme == :semix ? ALBEDO_SEMIX :
        ALBEDO_DYNAMIC
end

@inline function _normalize_seb_scheme(scheme::Symbol)
    scheme in (:bessi, :semix) ||
        error("Unsupported surface-energy-balance scheme '$scheme'. Use :bessi or :semix.")
    return scheme == :semix ? SEB_SEMIX : SEB_BESSI
end

@inline function _normalize_semix_snow_albedo(scheme::Symbol)
    scheme in (:ww, :warren, :warren_wiscombe, :dang) ||
        error("Unsupported `semix_snow_albedo` '$scheme'. Use :ww or :dang.")
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
    fresh_snow_density_scheme::Symbol=:constant,
    Ki::Real=2.1,
    ci::Real=2110.0,
    cw::Real=4181.0,
    Lm::Real=334000.0,
    Lv::Real=2.501e6,
    cp_air::Real=1003.0,
    latent_heat_flux_ratio::Real=1,
    D_sh::Real=10.0,
    alpha_dry::Real=0.81,
    alpha_wet::Real=0.60,
    alpha_ice::Real=0.3,
    max_lwc_albedo::Real=0.1,
    aging_cold_timescale_days::Real=20.0,
    aging_melting_timescale_days::Real=2.0,
    albedo_scheme::Symbol=:dynamic,
    ϵ_air::Real=0.8,
    ϵ_snow::Real=0.98,
    σ::Real=5.670373e-8,
    T0::Real=273.15,
    seconds_per_day::Real=DEFAULT_SECONDS_PER_DAY,
    low_density_densification::Symbol=:bessi,
    seb_scheme::Symbol=:bessi,
    turbulent_flux_scheme::Symbol=:bessi,
    eps_ice::Real=0.98,
    semix_karman::Real=0.4,
    semix_surface_height::Real=10.0,
    semix_z0m_snow::Real=0.001,
    semix_z0m_ice::Real=0.01,
    semix_zm_to_zh::Real=10.0,
    semix_sensible_exchange_factor::Real=1.0,
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
    if resolved_albedo_scheme == ALBEDO_AGING
        0 <= alpha_wet <= alpha_dry <= 1 || error(
            "Aging albedos must satisfy 0 <= alpha_wet <= alpha_dry <= 1.",
        )
    end
    aging_cold_timescale_days > 0 || error("`aging_cold_timescale_days` must be positive.")
    aging_melting_timescale_days > 0 || error("`aging_melting_timescale_days` must be positive.")
    fresh_snow_density_flag = _normalize_fresh_snow_density_scheme(fresh_snow_density_scheme)
    albedo_flag = _normalize_albedo_scheme(albedo_scheme)
    densification_flag = _normalize_low_density_densification(low_density_densification)
    seb_flag = _normalize_seb_scheme(seb_scheme)
    turbulent_flux_flag = _normalize_seb_scheme(turbulent_flux_scheme)
    semix_albedo_flag = _normalize_semix_snow_albedo(semix_snow_albedo)
    semix_karman > 0 || error("`semix_karman` must be positive.")
    semix_surface_height > 0 || error("`semix_surface_height` must be positive.")
    semix_z0m_snow > 0 || error("`semix_z0m_snow` must be positive.")
    semix_z0m_ice > 0 || error("`semix_z0m_ice` must be positive.")
    semix_zm_to_zh > 0 || error("`semix_zm_to_zh` must be positive.")
    semix_sensible_exchange_factor > 0 || error("`semix_sensible_exchange_factor` must be positive.")
    semix_latent_exchange_factor > 0 || error("`semix_latent_exchange_factor` must be positive.")
    fresh_snow_density_tag = fresh_snow_density_flag == FRESH_SNOW_DENSITY_CONSTANT ? :constant : :parameterized
    albedo_tag = albedo_flag == ALBEDO_CONSTANT ? :constant :
                 albedo_flag == ALBEDO_PRESCRIBED ? :prescribed :
                 albedo_flag == ALBEDO_AGING ? :aging :
                 albedo_flag == ALBEDO_SEMIX ? :semix : :dynamic
    densification_tag = densification_flag == LOW_DENSIFICATION_HTESSEL ? :htessel : :bessi
    return SnowpackPhysicalConstants{NF, fresh_snow_density_tag, albedo_tag, densification_tag}(
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

@adapt_structure SnowpackPhysicalConstants
