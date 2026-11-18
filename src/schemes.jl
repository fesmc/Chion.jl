#=
Exchangeable schemes. Every choice of a `BESSIModel` is a concrete subtype of an abstract scheme
type and is passed as an instance to `Params`, e.g. `Params(albedo_scheme = SemixAlbedo())`. A
scheme owns its coefficients (`AgingAlbedo(cold_timescale_days = 15)`); constants shared by several
schemes are fields of `Params`. The column step dispatches on the scheme instances, so a new scheme
is a new type plus the methods of its process (`src/processes/`).
=#

# ---------------------------------------------------------------------------
# Fresh-snow density
# ---------------------------------------------------------------------------

"""
Abstract supertype of the density of fresh snow. Pass an instance as
`Params(fresh_snow_density_scheme = …)`.

Available subtypes: [`ConstantFreshSnowDensity`](@ref), [`ParameterizedFreshSnowDensity`](@ref).
"""
abstract type AbstractFreshSnowDensityScheme end

"""
$(TYPEDEF)

Fresh snow takes a constant density, clamped to [50 kg m⁻³, `rho_i`].

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct ConstantFreshSnowDensity{NF <: AbstractFloat} <: AbstractFreshSnowDensityScheme
    "density of fresh snow (kg m⁻³, default `315`)"
    density::NF = 315.0
end

"""
$(TYPEDEF)

Fresh-snow density from the air temperature T_a and the wind speed U,
ρ_s = a + b (T_a − T₀) + c √U, clamped to [50 kg m⁻³, `rho_i`] (as in Crocus, Vionnet et al.,
2012).

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct ParameterizedFreshSnowDensity{NF <: AbstractFloat} <: AbstractFreshSnowDensityScheme
    "density at the melting point without wind, a (kg m⁻³, default `109`)"
    intercept::NF = 109.0
    "temperature coefficient b (kg m⁻³ K⁻¹, default `6`)"
    temperature_coefficient::NF = 6.0
    "wind coefficient c (kg m⁻³ (m s⁻¹)^-½, default `26`)"
    wind_coefficient::NF = 26.0
end

# ---------------------------------------------------------------------------
# Albedo
# ---------------------------------------------------------------------------

"""
Abstract supertype of the snow albedo. Pass an instance as `Params(albedo_scheme = …)`. The dry,
wet and ice albedos `alpha_dry`, `alpha_wet`, `alpha_ice` are fields of [`Params`](@ref), since
several schemes share them.

Available subtypes: [`DynamicAlbedo`](@ref), [`ConstantAlbedo`](@ref), [`AgingAlbedo`](@ref),
[`SemixAlbedo`](@ref), [`PrescribedAlbedo`](@ref).
"""
abstract type AbstractAlbedoScheme end

"""
$(TYPEDEF)

Fixed albedos: `alpha_dry` for a subfreezing snow surface, `alpha_wet` for a melting one and
`alpha_ice` for bare ice.
"""
struct ConstantAlbedo <: AbstractAlbedoScheme end

"""
$(TYPEDEF)

BESSI snow aging: the albedo decreases by (1.35·10⁻³ (T − T₀) + 0.0278) per day, relaxes toward
`alpha_wet` with the liquid-water content of the surface layer (saturated at `max_lwc_albedo`),
and brightens with snowfall, by (`alpha_dry` − `alpha_wet`)(1 − exp(−m_s / 3 kg m⁻²)).
"""
struct DynamicAlbedo <: AbstractAlbedoScheme end

"""
$(TYPEDEF)

The surface albedo of the `prescribed_albedo` forcing where it is given, and [`DynamicAlbedo`](@ref)
elsewhere.
"""
struct PrescribedAlbedo <: AbstractAlbedoScheme end

"""
$(TYPEDEF)

Exponential decay from `alpha_dry` toward `alpha_wet` with the time since the latest snowfall,
α = α_wet + (α − α_wet) exp(−Δt / τ), with τ the cold or the melting time scale (Oerlemans and
Knap, 1998). Snowfall resets the albedo to `alpha_dry`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct AgingAlbedo{NF <: AbstractFloat} <: AbstractAlbedoScheme
    "e-folding time of the albedo decay of a subfreezing surface (day, default `20`)"
    cold_timescale_days::NF = 20.0
    "e-folding time of the albedo decay of a melting surface (day, default `2`)"
    melting_timescale_days::NF = 2.0
end

"""
Abstract supertype of the spectral snow albedo of [`SemixAlbedo`](@ref).

Available subtypes: [`DangSnowAlbedo`](@ref), [`WarrenWiscombeSnowAlbedo`](@ref).
"""
abstract type AbstractSemixSnowAlbedo end

"""Dang et al. (2015) spectral snow albedo bands, including surface roughness."""
struct DangSnowAlbedo <: AbstractSemixSnowAlbedo end

"""Warren and Wiscombe (1980) spectral snow albedo bands."""
struct WarrenWiscombeSnowAlbedo <: AbstractSemixSnowAlbedo end

"""
$(TYPEDEF)

Two-band snow albedo of the SEMIX surface scheme: visible and near-infrared
albedos from the snow grain size, dust content and solar zenith angle, combined with the visible
fraction and the cloud cover. Bare ice takes `prescribed_ice_albedo` where the forcing gives it,
else `alpha_ice`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct SemixAlbedo{NF <: AbstractFloat, B <: AbstractSemixSnowAlbedo} <: AbstractAlbedoScheme
    "spectral snow albedo bands (default `DangSnowAlbedo()`)"
    snow_bands::B = DangSnowAlbedo()
    "visible fraction of the shortwave down under a clear sky (–, default `0.45`)"
    frac_vu::NF = 0.45
    "visible albedo of fresh snow (–, default `0.99`)"
    alb_snow_vis_new::NF = 0.99
    "near-infrared albedo of fresh snow (–, default `0.65`)"
    alb_snow_nir_new::NF = 0.65
    "grain size of fresh snow (µm, default `50`)"
    snow_grain_fresh::NF = 50.0
    "grain size of old snow (µm, default `1000`)"
    snow_grain_old::NF = 1000.0
    "visible albedo decrease of old snow (–, default `0.05`)"
    d_alb_age_vis::NF = 0.05
    "near-infrared albedo decrease of old snow (–, default `0.25`)"
    d_alb_age_nir::NF = 0.25
    "temperature factor of the aging rate (K⁻¹, default `0.1`)"
    f_age_t::NF = 0.1
    "temperature offset of the aging rate (K, default `0`)"
    dT_age::NF = 0.0
    "snowfall that refreshes the surface (kg m⁻² day⁻¹, default `1`)"
    snow_0::NF = 1.0
    "exponent of the snowfall refresh (–, default `0.5`)"
    snow_1::NF = 0.5
    "water equivalent over which melt concentrates dust (kg m⁻², default `10`)"
    w_snow_dust::NF = 10.0
    "scale of the dust concentration (–, default `1`)"
    dust_con_scale::NF = 1.0
    "visible albedo offset (–, default `0`)"
    dalb_snow_vis::NF = 0.0
    "near-infrared albedo offset (–, default `0`)"
    dalb_snow_nir::NF = 0.0
    "albedo reduction by sub-grid orography (–, default `0`)"
    k_sigma_orog::NF = 0.0
    "orography standard deviation of the half reduction (m, default `1000`)"
    sigma_orog_crit::NF = 1000.0
end

# ---------------------------------------------------------------------------
# Densification
# ---------------------------------------------------------------------------

"""
Abstract supertype of the densification of low-density snow. Pass an instance as
`Params(densification_scheme = …)`.

Available subtypes: [`BESSIDensification`](@ref), [`HTESSELDensification`](@ref).
"""
abstract type AbstractDensificationScheme end

"""$(TYPEDEF)

BESSI densification of snow and firn (Born et al., 2019)."""
struct BESSIDensification <: AbstractDensificationScheme end

"""$(TYPEDEF)

HTESSEL low-density densification (Dutra et al., 2010) with compaction by liquid water."""
struct HTESSELDensification <: AbstractDensificationScheme end

# ---------------------------------------------------------------------------
# Surface energy balance, turbulent fluxes, longwave
# ---------------------------------------------------------------------------

"""
Abstract supertype of the longwave exchange of the surface. Pass an instance as
`Params(surface_energy_balance_scheme = …)`.

Available subtypes: [`SemixSurfaceEnergyBalance`](@ref), [`BESSISurfaceEnergyBalance`](@ref).
"""
abstract type AbstractSurfaceEnergyBalanceScheme end

"""$(TYPEDEF)

SEMIX graybody exchange: the surface absorbs and emits longwave radiation with the snow or ice
emissivity ε, Q_LW = ε (L↓ − σ T_s⁴)."""
struct SemixSurfaceEnergyBalance <: AbstractSurfaceEnergyBalanceScheme end

"""$(TYPEDEF)

Original BESSI: all incident longwave is absorbed and the snow emissivity enters only the
emission, Q_LW = L↓ − ε_snow σ T_s⁴."""
struct BESSISurfaceEnergyBalance <: AbstractSurfaceEnergyBalanceScheme end

"""
Abstract supertype of the turbulent sensible and latent heat fluxes. Pass an instance as
`Params(turbulent_flux_scheme = …)`.

Available subtypes: [`SemixTurbulentFluxes`](@ref), [`BESSITurbulentFluxes`](@ref).
"""
abstract type AbstractTurbulentFluxScheme end

"""
$(TYPEDEF)

SEMIX bulk sensible and latent heat fluxes: a neutral transfer coefficient
C_H = κ² / (ln(z/z₀ₘ) ln(z/z₀ₕ)) corrected with the bulk Richardson number Ri
(1/(1 + b Ri) for stable, √(1 − 16 Ri) for unstable stratification), Q_H = f_H ρ c_p C_H U (T_a − T_s)
and Q_E = f_E L ρ C_H U (q_a − q_s).

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct SemixTurbulentFluxes{NF <: AbstractFloat} <: AbstractTurbulentFluxScheme
    "von Kármán constant κ (–, default `0.4`)"
    karman::NF = 0.4
    "height z of the air temperature, humidity and wind (m, default `10`)"
    reference_height::NF = 10.0
    "momentum roughness length of snow (m, default `0.001`)"
    z0m_snow::NF = 0.001
    "momentum roughness length of bare ice (m, default `0.01`)"
    z0m_ice::NF = 0.01
    "ratio of the momentum to the heat roughness length z₀ₘ/z₀ₕ (–, default `10`)"
    zm_to_zh::NF = 10.0
    "factor f_H of the sensible heat flux (–, default `3`, the GrIS calibration)"
    sensible_exchange_factor::NF = 3.0
    "stability coefficient b of stable stratification (–, default `40`)"
    stable_coefficient::NF = 40.0
    "factor f_E of the latent heat flux (–, default `1`)"
    latent_exchange_factor::NF = 1.0
end

"""
$(TYPEDEF)

BESSI turbulent fluxes: sensible heat D_sh (T_a − T_s) and latent heat from the vapour-pressure
gradient, L (r D_sh / c_p) (0.622/p) (e_a − e_s(T_s)).

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct BESSITurbulentFluxes{NF <: AbstractFloat} <: AbstractTurbulentFluxScheme
    "sensible heat transfer coefficient D_sh (W m⁻² K⁻¹, default `10`)"
    D_sh::NF = 10.0
    "ratio r of the latent to the sensible transfer coefficient (–, default `1`)"
    latent_heat_flux_ratio::NF = 1.0
end

"""
Abstract supertype of the longwave down where the forcing does not give it. Pass an instance as
`Params(longwave_scheme = …)`.

Available subtypes: [`CloudProxyLongwave`](@ref), [`GraybodyLongwave`](@ref).
"""
abstract type AbstractLongwaveScheme end

"""
$(TYPEDEF)

Longwave down ε σ T_a⁴ with an effective emissivity weighted between clear sky and overcast,
ε = ε_cs (1 − nᵖ) + ε_oc nᵖ, both linear in T_a − T₀, and the cloud fraction n diagnosed from the
shortwave transmissivity: the clear-sky transmissivity
τ_cs = min(1, τ₀ + τ_z z + τ_S ln(SW_TOA / 300 W m⁻²)) is reduced by clouds by the factor
1 − a nᵇ exp(−d z) (Konzelmann et al., 1994). Resolved once per daily forcing step; in the polar
night and for sub-daily forcing n is `night_cloud_fraction`.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct CloudProxyLongwave{NF <: AbstractFloat} <: AbstractLongwaveScheme
    "clear-sky transmissivity τ₀ at sea level and SW_TOA = 300 W m⁻² (–, default `0.718`)"
    clear_sky_transmissivity_sea_level::NF = 0.718
    "increase τ_z of the clear-sky transmissivity with height (km⁻¹, default `0.0298`)"
    clear_sky_transmissivity_per_km::NF = 0.0298
    "insolation coefficient τ_S of the clear-sky transmissivity (–, default `0.0569`)"
    clear_sky_transmissivity_insolation::NF = 0.0569
    "cloud transmission coefficient a (–, default `0.447`)"
    cloud_transmission_coefficient::NF = 0.447
    "cloud transmission exponent b (–, default `2.75`)"
    cloud_transmission_exponent::NF = 2.75
    "decay d of the cloud transmission with height (km⁻¹, default `0.225`)"
    cloud_transmission_decay_per_km::NF = 0.225
    "clear-sky emissivity at T₀ (–, default `0.716`)"
    clear_sky_emissivity::NF = 0.716
    "temperature slope of the clear-sky emissivity (K⁻¹, default `0.00239`)"
    clear_sky_emissivity_temperature_slope::NF = 0.00239
    "overcast emissivity at T₀ (–, default `0.952`)"
    overcast_emissivity::NF = 0.952
    "temperature slope of the overcast emissivity (K⁻¹, default `0.00362`)"
    overcast_emissivity_temperature_slope::NF = 0.00362
    "exponent p of the cloud fraction in the emissivity weight (–, default `3`)"
    cloud_emissivity_exponent::NF = 3.0
    "cloud fraction in the polar night and for sub-daily forcing (–, default `0.85`)"
    night_cloud_fraction::NF = 0.85
    "upper bound of the effective emissivity (–, default `1`)"
    emissivity_max::NF = 1.0
end

"""$(TYPEDEF)

Longwave down of a graybody sky, ϵ_air σ T_a⁴ (`ϵ_air` is a field of [`Params`](@ref))."""
struct GraybodyLongwave <: AbstractLongwaveScheme end

# ---------------------------------------------------------------------------
# Diurnal cycle
# ---------------------------------------------------------------------------

"""
Abstract supertype of the diurnal cycle within a daily forcing step. Pass an instance as
`Params(diurnal_cycle = …)`.

Available subtypes: [`DiurnalSubsteps`](@ref), [`NoDiurnalCycle`](@ref).
"""
abstract type AbstractDiurnalCycle end

"""$(TYPEDEF)

One column step per forcing step, with the daily means of the forcing."""
struct NoDiurnalCycle <: AbstractDiurnalCycle end

"""
Abstract supertype of the air-temperature cycle within the substeps of [`DiurnalSubsteps`](@ref).

Available subtypes: [`DiurnalTemperatureCycle`](@ref), [`NoDiurnalTemperatureCycle`](@ref).
"""
abstract type AbstractDiurnalTemperatureCycle end

"""$(TYPEDEF)

The substeps take the daily-mean air temperature."""
struct NoDiurnalTemperatureCycle <: AbstractDiurnalTemperatureCycle end

"""
$(TYPEDEF)

A sinusoidal air-temperature cycle of amplitude A = min(A₀ + γ max(z − z_ref, 0), A_max) about the
daily mean, peaking in the afternoon; each substep takes its interval average.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct DiurnalTemperatureCycle{NF <: AbstractFloat} <: AbstractDiurnalTemperatureCycle
    "amplitude A₀ (K, default `1`)"
    amplitude::NF = 1.0
    "increase γ of the amplitude with the surface height above `reference_height` (K m⁻¹, default `0`)"
    amplitude_gradient::NF = 0.0
    "height z_ref above which the amplitude increases (m, default `0`)"
    reference_height::NF = 0.0
    "upper bound A_max of the amplitude (K, default `1`)"
    amplitude_max::NF = 1.0
end

"""
$(TYPEDEF)

Up to `max_substeps` column steps per daily forcing step, each with the interval average of the
shortwave down reconstructed from the daily mean and the solar geometry (and of the air
temperature, with a [`DiurnalTemperatureCycle`](@ref)). Steps whose shortwave is below `threshold`
or whose air temperature is below `min_air_temperature` take one step.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct DiurnalSubsteps{NF <: AbstractFloat, T <: AbstractDiurnalTemperatureCycle} <: AbstractDiurnalCycle
    "maximum number of substeps per forcing step (1 to 24, default `8`)"
    max_substeps::Int = 8
    "daily shortwave below which a step is not split (W m⁻², default `0`)"
    threshold::NF = 0.0
    "air temperature below which a step is not split (K, default `273.15 − 8`)"
    min_air_temperature::NF = 273.15 - 8.0
    "air-temperature cycle of the substeps (default `DiurnalTemperatureCycle()`)"
    temperature_cycle::T = DiurnalTemperatureCycle()
end

# ---------------------------------------------------------------------------
# Bulk models
# ---------------------------------------------------------------------------

"""
Abstract supertype of the positive degree days of a [`PDDModel`](@ref). Pass an instance as
`PDDParams(method = …)`.

Available subtypes: [`SimplePDD`](@ref), [`PISMPDD`](@ref).
"""
abstract type AbstractPDDMethod end

"""$(TYPEDEF)

Positive degree days from the positive part of the mean temperature."""
struct SimplePDD <: AbstractPDDMethod end

"""$(TYPEDEF)

Calov–Greve expectation integral of a normally distributed temperature (as in PISM)."""
struct PISMPDD <: AbstractPDDMethod end

# ---------------------------------------------------------------------------
# Coefficients given as any real numbers (e.g. integer literals) become floats
# ---------------------------------------------------------------------------

for T in (:ConstantFreshSnowDensity, :ParameterizedFreshSnowDensity, :AgingAlbedo, :SemixTurbulentFluxes,
        :BESSITurbulentFluxes, :CloudProxyLongwave, :DiurnalTemperatureCycle)
    @eval $T(values::Real...) = $T(promote(map(float, values)...)...)
end
SemixAlbedo(snow_bands::AbstractSemixSnowAlbedo, values::Real...) =
    SemixAlbedo(snow_bands, promote(map(float, values)...)...)
DiurnalSubsteps(max_substeps::Integer, threshold::Real, min_air_temperature::Real,
        temperature_cycle::AbstractDiurnalTemperatureCycle) =
    DiurnalSubsteps(Int(max_substeps), promote(float(threshold), float(min_air_temperature))..., temperature_cycle)

# ---------------------------------------------------------------------------
# Validation and display
# ---------------------------------------------------------------------------

"""
    validate(scheme)

Throw an `ArgumentError` if the coefficients of `scheme` are out of range. Called by the `Params`
constructor; the default accepts everything.
"""
validate(scheme) = nothing

_require(condition::Bool, message) = condition || throw(ArgumentError(message))

function validate(s::AgingAlbedo)
    _require(s.cold_timescale_days > 0, "`cold_timescale_days` must be positive.")
    _require(s.melting_timescale_days > 0, "`melting_timescale_days` must be positive.")
    return nothing
end

function validate(s::SemixTurbulentFluxes)
    _require(s.karman > 0, "`karman` must be positive.")
    _require(s.reference_height > 0, "`reference_height` must be positive.")
    _require(s.z0m_snow > 0, "`z0m_snow` must be positive.")
    _require(s.z0m_ice > 0, "`z0m_ice` must be positive.")
    _require(s.zm_to_zh > 0, "`zm_to_zh` must be positive.")
    _require(s.sensible_exchange_factor > 0, "`sensible_exchange_factor` must be positive.")
    _require(s.stable_coefficient >= 0, "`stable_coefficient` must be nonnegative.")
    _require(s.latent_exchange_factor > 0, "`latent_exchange_factor` must be positive.")
    return nothing
end

function validate(s::CloudProxyLongwave)
    _require(s.clear_sky_transmissivity_sea_level > 0, "`clear_sky_transmissivity_sea_level` must be positive.")
    _require(s.cloud_transmission_coefficient > 0, "`cloud_transmission_coefficient` must be positive.")
    _require(s.cloud_transmission_exponent > 0, "`cloud_transmission_exponent` must be positive.")
    _require(s.cloud_emissivity_exponent > 0, "`cloud_emissivity_exponent` must be positive.")
    _require(0 <= s.night_cloud_fraction <= 1, "`night_cloud_fraction` must be in [0, 1].")
    _require(0.4 < s.emissivity_max <= 1, "`emissivity_max` must be in (0.4, 1].")
    return nothing
end

function validate(s::DiurnalSubsteps)
    _require(1 <= s.max_substeps <= 24, "`max_substeps` must be between 1 and 24.")
    _require(s.threshold >= 0, "`threshold` must be non-negative.")
    _require(isfinite(s.min_air_temperature), "`min_air_temperature` must be finite.")
    validate(s.temperature_cycle)
    return nothing
end

function validate(s::DiurnalTemperatureCycle)
    _require(s.amplitude >= 0, "`amplitude` must be non-negative.")
    _require(isfinite(s.amplitude_gradient), "`amplitude_gradient` must be finite.")
    _require(s.amplitude_max >= s.amplitude, "`amplitude_max` must be at least `amplitude`.")
    return nothing
end

# Compact display: the type name and the coefficients that differ from the defaults.
function Base.show(io::IO, scheme::Union{ConstantFreshSnowDensity, ParameterizedFreshSnowDensity, AgingAlbedo,
        SemixAlbedo, SemixTurbulentFluxes, BESSITurbulentFluxes, CloudProxyLongwave, DiurnalTemperatureCycle,
        DiurnalSubsteps})
    T = typeof(scheme)
    default = _promote_param(Base.typename(T).wrapper(), _scheme_float_type(scheme))
    changed = [name for name in fieldnames(T) if !isequal(getfield(scheme, name), getfield(default, name))]
    print(io, nameof(T), "(")
    join(io, ("$name = $(getfield(scheme, name))" for name in changed), ", ")
    print(io, ")")
    return nothing
end

# Floating-point type of the first float field of a scheme (`Float64` if it has none).
function _scheme_float_type(scheme)
    for name in fieldnames(typeof(scheme))
        value = getfield(scheme, name)
        value isa AbstractFloat && return typeof(value)
    end
    return Float64
end
