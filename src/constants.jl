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

abstract type AbstractSurfaceEnergyBalanceScheme end
"""SEMIX graybody exchange: the surface absorbs and emits longwave radiation with the snow or
ice emissivity ε, Q_LW = ε (L↓ − σ T_s⁴)."""
struct SemixSurfaceEnergyBalance <: AbstractSurfaceEnergyBalanceScheme end
"""Original BESSI: all incident longwave is absorbed and the snow emissivity enters only the
emission, Q_LW = L↓ − ε_snow σ T_s⁴."""
struct BESSISurfaceEnergyBalance <: AbstractSurfaceEnergyBalanceScheme end

abstract type AbstractTurbulentFluxScheme end
"""SEMIX bulk sensible and latent heat with Richardson-number stability
(`_semix_turbulent_flux_linearized`)."""
struct SemixTurbulentFluxes <: AbstractTurbulentFluxScheme end
"""BESSI sensible heat `D_sh (T_a − T_s)` and vapour-gradient latent heat."""
struct BESSITurbulentFluxes <: AbstractTurbulentFluxScheme end

abstract type AbstractLongwaveScheme end
"""Longwave down from a cloud fraction diagnosed from the shortwave transmissivity
(Konzelmann et al., 1994), resolved once per daily forcing step."""
struct CloudProxyLongwave <: AbstractLongwaveScheme end
"""Longwave down of a graybody sky, `ϵ_air σ T_a⁴`."""
struct GraybodyLongwave <: AbstractLongwaveScheme end

abstract type AbstractSemixSnowAlbedo end
"""Dang et al. (2015) spectral snow albedo bands, including surface roughness."""
struct DangSnowAlbedo <: AbstractSemixSnowAlbedo end
"""Warren and Wiscombe (1980) spectral snow albedo bands."""
struct WarrenWiscombeSnowAlbedo <: AbstractSemixSnowAlbedo end

# ---------------------------------------------------------------------------
# Scheme selection: names (keywords, TOML) are converted once, by dispatch
# ---------------------------------------------------------------------------

"""
    _scheme(kind, scheme)

The scheme instance of abstract type `kind` that `scheme` names: an instance of `kind` is returned
as it is, a `Symbol` or string is converted by the `_scheme(kind, ::Val{name})` methods. An
unknown name lists the available ones.
"""
@inline _scheme(::Type{K}, scheme::K) where {K} = scheme
_scheme(::Type{K}, name::Union{Symbol, AbstractString}) where {K} = _scheme(K, Val(Symbol(name)))
function _scheme(::Type{K}, ::Val{name}) where {K, name}
    names = join((":" * String(n) for n in _scheme_names(K)), ", ")
    throw(ArgumentError("Unsupported $(_scheme_description(K)) '$name'. Use $names."))
end
function _scheme(::Type{K}, scheme) where {K}
    throw(ArgumentError("Expected a $(_scheme_description(K)) (a name or an `$(nameof(K))`), got $(repr(scheme))."))
end

# One table of the scheme names: it defines the `_scheme(kind, ::Val{name})` methods and the
# names listed in errors.
const _SCHEME_TABLE = (
    (AbstractFreshSnowDensityScheme, "fresh-snow density scheme",
        (constant = ConstantFreshSnowDensity(), parameterized = ParameterizedFreshSnowDensity())),
    (AbstractAlbedoScheme, "albedo scheme",
        (constant = ConstantAlbedo(), dynamic = DynamicAlbedo(), prescribed = PrescribedAlbedo(),
         aging = AgingAlbedo(), semix = SemixAlbedo())),
    (AbstractDensificationScheme, "low-density densification scheme",
        (bessi = BESSIDensification(), htessel = HTESSELDensification())),
    (AbstractSurfaceEnergyBalanceScheme, "surface-energy-balance scheme",
        (semix = SemixSurfaceEnergyBalance(), bessi = BESSISurfaceEnergyBalance())),
    (AbstractTurbulentFluxScheme, "turbulent flux scheme",
        (semix = SemixTurbulentFluxes(), bessi = BESSITurbulentFluxes())),
    (AbstractLongwaveScheme, "longwave scheme",
        (cloud_proxy = CloudProxyLongwave(), graybody = GraybodyLongwave())),
    (AbstractSemixSnowAlbedo, "`semix_snow_albedo`",
        (dang = DangSnowAlbedo(), warren_wiscombe = WarrenWiscombeSnowAlbedo())),
)

for (kind, description, schemes) in _SCHEME_TABLE
    @eval _scheme_names(::Type{$kind}) = $(keys(schemes))
    @eval _scheme_description(::Type{$kind}) = $description
    for (name, scheme) in pairs(schemes)
        @eval _scheme(::Type{$kind}, ::Val{$(QuoteNode(name))}) = $scheme
    end
end

"""
    _scheme_name(scheme)

The keyword/TOML name of the scheme instance `scheme` (e.g. `:semix` for
`SemixSurfaceEnergyBalance()`), or its type name for a scheme without one. Host side only.
"""
function _scheme_name(scheme)
    for (kind, _, schemes) in _SCHEME_TABLE
        scheme isa kind || continue
        for (name, candidate) in pairs(schemes)
            candidate === scheme && return name
        end
    end
    return nameof(typeof(scheme))
end

"""Check scheme-specific constraints on the constants; the default accepts everything."""
_validate_albedo_scheme(::AbstractAlbedoScheme, alpha_dry, alpha_wet) = nothing
_validate_albedo_scheme(::AgingAlbedo, alpha_dry, alpha_wet) =
    0 <= alpha_wet <= alpha_dry <= 1 ||
    throw(ArgumentError("Aging albedos must satisfy 0 <= alpha_wet <= alpha_dry <= 1."))

# ---------------------------------------------------------------------------
# Physical constants struct and constructors
# ---------------------------------------------------------------------------

"""
    SnowpackPhysicalConstants{NF, FreshSnowDensity, Albedo, Densification, SurfaceEnergyBalance,
                              TurbulentFluxes, Longwave, SemixSnowAlbedo}

Container for physical constants and empirical coefficients used by the snowpack model. Every
scheme choice is a type parameter (see `AbstractAlbedoScheme` etc.), so process code selects
the scheme by dispatch, at compile time.
"""
struct SnowpackPhysicalConstants{
        NF <: AbstractFloat,
        FreshSnowDensity <: AbstractFreshSnowDensityScheme,
        Albedo <: AbstractAlbedoScheme,
        Densification <: AbstractDensificationScheme,
        SurfaceEnergyBalance <: AbstractSurfaceEnergyBalanceScheme,
        TurbulentFluxes <: AbstractTurbulentFluxScheme,
        Longwave <: AbstractLongwaveScheme,
        SemixSnowAlbedo <: AbstractSemixSnowAlbedo,
    }
    rho_s::NF
    rho_i::NF
    rho_w::NF
    rho_s_a::NF
    rho_s_b::NF
    rho_s_c::NF
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
    eps_ice::NF
    semix_karman::NF
    semix_surface_height::NF
    semix_z0m_snow::NF
    semix_z0m_ice::NF
    semix_zm_to_zh::NF
    semix_sensible_exchange_factor::NF
    semix_stable_coefficient::NF
    semix_latent_exchange_factor::NF
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
"""Surface-energy-balance scheme of `c` (an `AbstractSurfaceEnergyBalanceScheme` instance)."""
@inline surface_energy_balance_scheme(::SnowpackPhysicalConstants{NF, F, A, D, S}) where {NF, F, A, D, S} = S()
"""Turbulent flux scheme of `c` (an `AbstractTurbulentFluxScheme` instance)."""
@inline turbulent_flux_scheme(::SnowpackPhysicalConstants{NF, F, A, D, S, T}) where {NF, F, A, D, S, T} = T()
"""Longwave scheme of `c` (an `AbstractLongwaveScheme` instance)."""
@inline longwave_scheme(::SnowpackPhysicalConstants{NF, F, A, D, S, T, L}) where {NF, F, A, D, S, T, L} = L()
"""SEMIX snow albedo bands of `c` (an `AbstractSemixSnowAlbedo` instance)."""
@inline semix_snow_albedo(::SnowpackPhysicalConstants{NF, F, A, D, S, T, L, W}) where {NF, F, A, D, S, T, L, W} = W()

const _SCHEME_PROPERTIES = (:fresh_snow_density_scheme, :albedo_scheme, :low_density_densification,
    :seb_scheme, :turbulent_flux_scheme, :longwave_scheme, :semix_snow_albedo)

@inline function Base.getproperty(c::SnowpackPhysicalConstants, name::Symbol)
    name === :fresh_snow_density_scheme && return fresh_snow_density_scheme(c)
    name === :albedo_scheme && return albedo_scheme(c)
    name === :low_density_densification && return densification_scheme(c)
    name === :seb_scheme && return surface_energy_balance_scheme(c)
    name === :turbulent_flux_scheme && return turbulent_flux_scheme(c)
    name === :longwave_scheme && return longwave_scheme(c)
    name === :semix_snow_albedo && return semix_snow_albedo(c)
    return getfield(c, name)
end

Base.propertynames(c::SnowpackPhysicalConstants, private::Bool=false) =
    (fieldnames(typeof(c))..., _SCHEME_PROPERTIES...)

Base.eltype(::SnowpackPhysicalConstants{NF}) where {NF} = NF

"""
    number_type(c)

Return the floating-point element type used by the physical constants set `c`.
"""
@inline number_type(::SnowpackPhysicalConstants{NF}) where {NF} = NF

@inline _initial_snow_albedo(c::SnowpackPhysicalConstants) =
    c.alpha_dry

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
    longwave_scheme::Union{Symbol, AbstractLongwaveScheme}=:cloud_proxy,
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
    seb_scheme::Union{Symbol, AbstractSurfaceEnergyBalanceScheme}=:semix,
    turbulent_flux_scheme::Union{Symbol, AbstractTurbulentFluxScheme}=:semix,
    eps_ice::Real=0.98,
    semix_karman::Real=0.4,
    semix_surface_height::Real=10.0,
    semix_z0m_snow::Real=0.001,
    semix_z0m_ice::Real=0.01,
    semix_zm_to_zh::Real=10.0,
    semix_sensible_exchange_factor::Real=3.0,
    semix_stable_coefficient::Real=40.0,
    semix_latent_exchange_factor::Real=1.0,
    semix_snow_albedo::Union{Symbol, AbstractSemixSnowAlbedo}=:dang,
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
    schemes = (
        _scheme(AbstractFreshSnowDensityScheme, fresh_snow_density_scheme),
        _scheme(AbstractAlbedoScheme, albedo_scheme),
        _scheme(AbstractDensificationScheme, low_density_densification),
        _scheme(AbstractSurfaceEnergyBalanceScheme, seb_scheme),
        _scheme(AbstractTurbulentFluxScheme, turbulent_flux_scheme),
        _scheme(AbstractLongwaveScheme, longwave_scheme),
        _scheme(AbstractSemixSnowAlbedo, semix_snow_albedo),
    )
    _validate_albedo_scheme(schemes[2], alpha_dry, alpha_wet)
    aging_cold_timescale_days > 0 || error("`aging_cold_timescale_days` must be positive.")
    aging_melting_timescale_days > 0 || error("`aging_melting_timescale_days` must be positive.")
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
    return SnowpackPhysicalConstants{NF, map(typeof, schemes)...}(
        convert(NF, rho_s),
        convert(NF, rho_i),
        convert(NF, rho_w),
        convert(NF, rho_s_a),
        convert(NF, rho_s_b),
        convert(NF, rho_s_c),
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
        convert(NF, eps_ice),
        convert(NF, semix_karman),
        convert(NF, semix_surface_height),
        convert(NF, semix_z0m_snow),
        convert(NF, semix_z0m_ice),
        convert(NF, semix_zm_to_zh),
        convert(NF, semix_sensible_exchange_factor),
        convert(NF, semix_stable_coefficient),
        convert(NF, semix_latent_exchange_factor),
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
