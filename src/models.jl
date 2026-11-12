"""
Model definitions exposed by Chion's simulation-first API.
"""

abstract type AbstractSnowModel end

"""Immutable BESSI configuration shared by the model and its state containers.
Albedo schemes include `:constant`, `:dynamic`, `:prescribed`, `:aging`, and `:semix`.
"""
struct BESSIParameters{C <: SnowpackPhysicalConstants, NF <: AbstractFloat}
    c::C
    Ntot::Int
    mass_max::NF
    mass_split::NF
    mass_min::NF
    density_init::NF
    temperature_init::NF
    refreezing_correction::NF
    diurnal_shortwave_threshold::NF
    diurnal_shortwave_max_substeps::Int
    diurnal_shortwave_min_air_temperature::NF
    diurnal_temperature_amplitude::NF
    diurnal_temperature_amplitude_gradient::NF
    diurnal_temperature_amplitude_reference_height::NF
    diurnal_temperature_amplitude_max::NF
    near_surface_layer_max_thicknesses_m::NTuple{4, NF}
    ice_substrate_layers::Int
    ice_substrate_top_thickness_m::NF
end

const _BESSI_MODEL_TAG_PROPERTIES = (:diurnal_shortwave_substeps, :diurnal_temperature_cycle)

"""
    BESSIModel(grid; albedo=:dynamic, seb_scheme=:semix,
               turbulent_flux_scheme=:semix, densification=:bessi, ...)

Configuration for the layered BESSI snowpack model. Evolving state is stored in
`BESSIState` and owned by `Simulation.now`. Turbulent sensible and latent heat
use the SEMIX formulation by default; set `turbulent_flux_scheme=:bessi` to use
the BESSI formulation. `seb_scheme` independently controls the remaining
surface-energy-balance details and defaults to `:semix` for consistent graybody
longwave exchange. Set `seb_scheme=:bessi` for the original BESSI formulation.

The defaults are the GrIS-calibrated setup for daily forcing: dynamic albedo,
cloud-proxy longwave, a 5-layer thermal ice substrate, fine near-surface layers
and up to 8 diurnal substeps with a 1 K air-temperature cycle.
"""
struct BESSIModel{
        DiurnalShortwave,
        DiurnalTemperature,
        G <: SnowpackGrid,
        P <: BESSIParameters,
    } <: AbstractSnowModel
    grid::G
    parameters::P
end

@inline _diurnal_shortwave_enabled(::BESSIModel{Enabled}) where {Enabled} = Enabled
@inline _diurnal_temperature_enabled(::BESSIModel{Shortwave, Enabled}) where {Shortwave, Enabled} = Enabled

@inline function Base.getproperty(model::BESSIModel, name::Symbol)
    name === :grid && return getfield(model, :grid)
    name === :parameters && return getfield(model, :parameters)
    name === :diurnal_shortwave_substeps && return _diurnal_shortwave_enabled(model)
    name === :diurnal_temperature_cycle && return _diurnal_temperature_enabled(model)
    return getproperty(getfield(model, :parameters), name)
end

function Base.propertynames(model::BESSIModel, private::Bool=false)
    public = (:grid, fieldnames(typeof(getfield(model, :parameters)))..., _BESSI_MODEL_TAG_PROPERTIES...)
    return private ? (public..., :parameters) : public
end

function BESSIModel(
    grid::SnowpackGrid;
    albedo::Union{Symbol, AbstractAlbedoScheme}=:dynamic,
    seb_scheme::Symbol=:semix,
    turbulent_flux_scheme::Symbol=:semix,
    densification::Union{Symbol, AbstractDensificationScheme}=:bessi,
    fresh_snow_density::Union{Symbol, AbstractFreshSnowDensityScheme}=:constant,
    Ntot::Int=DEFAULT_NTOT,
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    density_init::Real=DEFAULT_DENSITY_INIT,
    temperature_init::Real=DEFAULT_TEMPERATURE_INIT,
    refreezing_correction::Real=1.0,
    diurnal_shortwave_substeps::Bool=true,
    diurnal_shortwave_threshold::Real=0.0,
    diurnal_shortwave_max_substeps::Integer=8,
    diurnal_shortwave_min_air_temperature_c::Real=-8.0,
    diurnal_temperature_cycle::Bool=true,
    diurnal_temperature_amplitude_c::Real=1.0,
    diurnal_temperature_amplitude_gradient_c_per_km::Real=0.0,
    diurnal_temperature_amplitude_reference_height_m::Real=0.0,
    diurnal_temperature_amplitude_max_c::Real=1.0,
    near_surface_layer_max_thicknesses_m::NTuple{4, <:Real}=(0.02, 0.05, 0.10, 0.30),
    ice_substrate_layers::Integer=5,
    ice_substrate_top_thickness_m::Real=0.05,
    number_type::Type{<:AbstractFloat}=Float64,
    kwargs...,
)
    isfinite(refreezing_correction) && refreezing_correction > 0 ||
        error("`refreezing_correction` must be finite and positive.")
    diurnal_shortwave_threshold >= 0 || error("`diurnal_shortwave_threshold` must be non-negative.")
    1 <= diurnal_shortwave_max_substeps <= 24 || error("`diurnal_shortwave_max_substeps` must be between 1 and 24.")
    isfinite(diurnal_shortwave_min_air_temperature_c) || error("`diurnal_shortwave_min_air_temperature_c` must be finite.")
    diurnal_temperature_amplitude_c >= 0 || error("`diurnal_temperature_amplitude_c` must be non-negative.")
    isfinite(diurnal_temperature_amplitude_gradient_c_per_km) ||
        error("`diurnal_temperature_amplitude_gradient_c_per_km` must be finite.")
    diurnal_temperature_amplitude_max_c >= diurnal_temperature_amplitude_c ||
        error("`diurnal_temperature_amplitude_max_c` must be at least the base amplitude.")
    all(thickness -> thickness > 0, near_surface_layer_max_thicknesses_m) ||
        error("`near_surface_layer_max_thicknesses_m` must contain positive values.")
    ice_substrate_layers >= 0 || error("`ice_substrate_layers` must be non-negative.")
    ice_substrate_top_thickness_m > 0 || error("`ice_substrate_top_thickness_m` must be positive.")
    c = SnowpackPhysicalConstants(
        number_type;
        albedo_scheme=albedo,
        seb_scheme=seb_scheme,
        turbulent_flux_scheme=turbulent_flux_scheme,
        low_density_densification=densification,
        fresh_snow_density_scheme=fresh_snow_density,
        kwargs...,
    )
    mass_max, mass_split, mass_min =
        _domain_thresholds(number_type, mass_max, mass_split, mass_min)
    parameters = BESSIParameters(
        c,
        Ntot,
        mass_max,
        mass_split,
        mass_min,
        number_type(density_init),
        number_type(temperature_init),
        number_type(refreezing_correction),
        number_type(diurnal_shortwave_threshold),
        Int(diurnal_shortwave_max_substeps),
        number_type(diurnal_shortwave_min_air_temperature_c + 273.15),
        number_type(diurnal_temperature_amplitude_c),
        number_type(diurnal_temperature_amplitude_gradient_c_per_km / 1000),
        number_type(diurnal_temperature_amplitude_reference_height_m),
        number_type(diurnal_temperature_amplitude_max_c),
        ntuple(i -> number_type(near_surface_layer_max_thicknesses_m[i]), 4),
        Int(ice_substrate_layers),
        number_type(ice_substrate_top_thickness_m),
    )
    return BESSIModel{
        diurnal_shortwave_substeps,
        diurnal_temperature_cycle,
        typeof(grid),
        typeof(parameters),
    }(
        grid,
        parameters,
    )
end

abstract type AbstractPDDMethod end
"""Positive degree days from the positive part of the mean temperature."""
struct SimplePDD <: AbstractPDDMethod end
"""Calov-Greve expectation integral of a normally distributed temperature (as in PISM)."""
struct PISMPDD <: AbstractPDDMethod end

const _PDD_METHODS = (simple=SimplePDD(), pism=PISMPDD())

@inline _normalize_pdd_method(method::AbstractPDDMethod) = method
@inline function _normalize_pdd_method(method::Symbol)
    haskey(_PDD_METHODS, method) || error("Unsupported PDD method '$method'. Use :simple or :pism.")
    return _PDD_METHODS[method]
end

"""
    PDDModel(grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.6,
             temperature_sigma=5.0, H_snow_max=5000.0, pdd_method=:simple,
             number_type=Float64)

Bulk positive-degree-day surface mass-balance model configuration. The model
uses a capped one-layer snow reservoir: refrozen water and snow above
`H_snow_max` become ice, while `smb_ice` records only ice-facing mass changes.
Set `pdd_method=:pism` to use the Calov-Greve expectation integral for every
timestep; `:simple` uses positive mean temperature directly.
`number_type` sets the floating-point type of the parameters and state.
"""
struct PDDModel{
        Method <: AbstractPDDMethod,
        NF <: AbstractFloat,
        G <: SnowpackGrid,
        C <: SnowpackPhysicalConstants{NF},
    } <: AbstractSnowModel
    grid::G
    c::C
    ddf_snow::NF
    ddf_ice::NF
    refreezing_fraction::NF
    temperature_sigma::NF
    H_snow_max::NF
end
"""PDD method of `model` (an `AbstractPDDMethod` instance)."""
@inline pdd_method(::PDDModel{Method}) where {Method} = Method()

@inline function Base.getproperty(model::PDDModel, name::Symbol)
    name === :pdd_method && return pdd_method(model)
    return getfield(model, name)
end

Base.propertynames(model::PDDModel, private::Bool=false) =
    (fieldnames(typeof(model))..., :pdd_method)

"""
    ITMModel(grid; kwargs...)

Insolation-temperature-melt bulk surface-mass-balance model, ported from the
Fortran `fesmc/chion`/smbpal implementation. ITM requires explicit
`surface_height`, `ice_thickness`, `annual_pdd`, and `latitude_deg` fields in
`SnowpackForcing`; `q_sw_net` is used as its insolation driver when supplied,
otherwise `shortwave_down` is used.

The ITM offset depends on latitude as `itm_c + itm_b*(|lat| - |itm_lat0|)`, the
same in both hemispheres: `itm_b` [W m^-2 deg^-1] is the slope in |latitude| and
`itm_lat0` [deg] the reference latitude (its absolute value; either hemisphere).
`|itm_lat0| >= 90` disables the latitude dependence. The defaults
(`itm_c = -45`, `itm_b = -2`, `itm_lat0 = 65`) are Greenland-tuned.
"""
struct ITMModel{
        NF <: AbstractFloat,
        G <: SnowpackGrid,
        C <: SnowpackPhysicalConstants{NF},
    } <: AbstractSnowModel
    grid::G
    c::C
    trans_a::NF
    trans_b::NF
    itm_c::NF
    itm_t::NF
    itm_b::NF
    itm_lat0::NF
    H_snow_max::NF
    Pmaxfrac::NF
    H_snow_crit_desert::NF
    H_snow_crit_forest::NF
    melt_crit::NF
    alb_ocean::NF
    alb_land::NF
    alb_forest::NF
    alb_ice::NF
    alb_snow_dry::NF
    alb_snow_wet::NF
    firn_fac::NF
end

function ITMModel(
    grid::SnowpackGrid;
    trans_a::Real=0.46,
    trans_b::Real=6e-5,
    itm_c::Real=-45.0,
    itm_t::Real=10.0,
    itm_b::Real=-2.0,
    itm_lat0::Real=65.0,
    H_snow_max::Real=5000.0,
    Pmaxfrac::Real=0.6,
    H_snow_crit_desert::Real=10.0,
    H_snow_crit_forest::Real=100.0,
    melt_crit::Real=0.5,
    alb_ocean::Real=0.1,
    alb_land::Real=0.2,
    alb_forest::Real=0.1,
    alb_ice::Real=0.4,
    alb_snow_dry::Real=0.8,
    alb_snow_wet::Real=0.65,
    firn_fac::Real=0.0266,
    number_type::Type{<:AbstractFloat}=Float64,
    kwargs...,
)
    H_snow_max > 0 || error("`H_snow_max` must be positive.")
    H_snow_crit_desert > 0 || error("`H_snow_crit_desert` must be positive.")
    H_snow_crit_forest > 0 || error("`H_snow_crit_forest` must be positive.")
    0 <= Pmaxfrac <= 1 || error("`Pmaxfrac` must be between 0 and 1.")
    all(value -> 0 <= value <= 1, (alb_ocean, alb_land, alb_forest, alb_ice, alb_snow_dry, alb_snow_wet)) ||
        error("ITM albedos must be between 0 and 1.")
    NF = number_type
    return ITMModel(grid, SnowpackPhysicalConstants(NF; kwargs...), NF(trans_a),
        NF(trans_b), NF(itm_c), NF(itm_t), NF(itm_b),
        NF(itm_lat0), NF(H_snow_max), NF(Pmaxfrac),
        NF(H_snow_crit_desert), NF(H_snow_crit_forest), NF(melt_crit),
        NF(alb_ocean), NF(alb_land), NF(alb_forest), NF(alb_ice),
        NF(alb_snow_dry), NF(alb_snow_wet), NF(firn_fac))
end

function PDDModel(
    grid::SnowpackGrid;
    ddf_snow::Real=3.0,
    ddf_ice::Real=8.0,
    refreezing_fraction::Real=0.6,
    temperature_sigma::Real=5.0,
    H_snow_max::Real=5000.0,
    pdd_method::Union{Symbol, AbstractPDDMethod}=:simple,
    number_type::Type{<:AbstractFloat}=Float64,
    kwargs...,
)
    ddf_snow > 0 || error("`ddf_snow` must be positive.")
    ddf_ice >= 0 || error("`ddf_ice` must be non-negative.")
    0.0 <= refreezing_fraction <= 1.0 || error("`refreezing_fraction` must be between 0 and 1.")
    temperature_sigma > 0 || error("`temperature_sigma` must be positive.")
    H_snow_max > 0 || error("`H_snow_max` must be positive.")

    method = _normalize_pdd_method(pdd_method)
    NF = number_type
    c = SnowpackPhysicalConstants(NF; kwargs...)
    return PDDModel{typeof(method), NF, typeof(grid), typeof(c)}(
        grid,
        c,
        NF(ddf_snow),
        NF(ddf_ice),
        NF(refreezing_fraction),
        NF(temperature_sigma),
        NF(H_snow_max),
    )
end
