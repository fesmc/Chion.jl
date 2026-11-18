#=
Physical constants and scheme choices of the models, one immutable, concretely typed struct per
model (as Laddie's `Params`): `Params` for `BESSIModel`, `PDDParams` and `ITMParams` for the bulk
models. Every keyword of the constructors is a field; unknown keywords throw.
=#

# ---------------------------------------------------------------------------
# Number-type promotion
# ---------------------------------------------------------------------------

# A struct of type `T` rebuilt from new field values: the type's name, so that the type parameters
# follow the new values (e.g. `SemixAlbedo{Float64}` → `SemixAlbedo{Float32}`).
_rebuild(::Type{T}, values...) where {T} = Base.typename(T).wrapper(values...)

"""`x` with `f` applied to every field (singletons are returned as they are)."""
function _mapfields(f, x)
    T = typeof(x)
    fieldcount(T) == 0 && return x
    return _rebuild(T, ntuple(i -> f(getfield(x, i)), fieldcount(T))...)
end

_to_nf(value::Real, ::Type{NF}) where {NF} = NF(value)
_to_nf(value::Tuple, ::Type{NF}) where {NF} = map(v -> _to_nf(v, NF), value)

# Fields declared with the number-type parameter `NF` (or as a tuple of it), read from the struct
# definition; fields declared otherwise (e.g. `Ntot::Int`) keep their values.
_is_number_type_parameter(declared) = declared isa TypeVar && declared.name === :NF
_is_number_type_field(declared) = _is_number_type_parameter(declared) ||
    (declared isa DataType && declared <: Tuple && !isempty(declared.parameters) &&
     all(_is_number_type_parameter, declared.parameters))

"""
    _promote_param(x, NF)

A scheme, boundary condition or `Layering` with every field declared with its number-type
parameter `NF` (integers included) converted to `NF`, nested schemes recursively.
"""
function _promote_param(x, ::Type{NF}) where {NF}
    T = typeof(x)
    fieldcount(T) == 0 && return x
    declared = Base.unwrap_unionall(Base.typename(T).wrapper).types
    values = ntuple(fieldcount(T)) do i
        value = getfield(x, i)
        _is_number_type_field(declared[i]) ? _to_nf(value, NF) :
            !(value isa Number) && isstructtype(typeof(value)) ? _promote_param(value, NF) : value
    end
    return _rebuild(T, values...)
end

# A field value of a parameter struct at type NF: a number converted, a scheme promoted.
_param_value(v::Integer, ::Type{NF}, ::Type{Int}) where {NF} = Int(v)
_param_value(v::Number, ::Type{NF}, ::Type) where {NF} = NF(v)
_param_value(v, ::Type{NF}, ::Type) where {NF} = _promote_param(v, NF)

"""Build a parameter struct `P` from its defaults and the keywords; unknown keywords throw."""
function _build_params(::Type{P}, defaults::NamedTuple, ::Type{NF}, kwargs) where {P, NF}
    for name in keys(kwargs)
        name in keys(defaults) || throw(ArgumentError("$(nameof(P)) has no field `$name`. " *
            "The fields are $(join(("`$n`" for n in keys(defaults)), ", "))."))
    end
    values = merge(defaults, kwargs)
    # Integer fields stay `Int` (e.g. a count); every other number becomes NF.
    integer_fields = _integer_fields(P)
    return _rebuild(P, map(name -> _param_value(values[name], NF, name in integer_fields ? Int : NF),
        keys(defaults))...)
end
_integer_fields(::Type) = ()

# ---------------------------------------------------------------------------
# BESSI
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

All physical constants and scheme choices of a [`BESSIModel`](@ref), in one immutable, concretely
typed struct. Build it with the keyword constructor, which fills in the defaults, the Greenland
calibration against MAR v3.14.3 (daily forcing):

```jldoctest
julia> params = Params(albedo_scheme = AgingAlbedo(cold_timescale_days = 15), NF = Float32);

julia> params.albedo_scheme
AgingAlbedo(cold_timescale_days = 15.0)

julia> params.alpha_dry
0.81f0
```

Every scheme slot takes an instance of its abstract type; the instance carries the coefficients of
that scheme. [`ParamsAISv1`](@ref) gives the Antarctic calibration.

# Fields
$(TYPEDFIELDS)
"""
struct Params{
        NF <: Real,
        FS <: AbstractFreshSnowDensityScheme,
        A <: AbstractAlbedoScheme,
        D <: AbstractDensificationScheme,
        S <: AbstractSurfaceEnergyBalanceScheme,
        T <: AbstractTurbulentFluxScheme,
        L <: AbstractLongwaveScheme,
        DC <: AbstractDiurnalCycle,
    }
    # Material constants
    "ice density (kg m⁻³, default `917`)"
    rho_i::NF
    "liquid water density (kg m⁻³, default `1000`)"
    rho_w::NF
    "specific heat capacity of ice (J kg⁻¹ K⁻¹, default `2110`)"
    ci::NF
    "specific heat capacity of liquid water (J kg⁻¹ K⁻¹, default `4181`)"
    cw::NF
    "latent heat of fusion (J kg⁻¹, default `3.34e5`)"
    Lm::NF
    "latent heat of vaporisation (J kg⁻¹, default `2.501e6`)"
    Lv::NF
    "specific heat capacity of air (J kg⁻¹ K⁻¹, default `1003`)"
    cp_air::NF
    "Stefan–Boltzmann constant (W m⁻² K⁻⁴, default `5.670373e-8`)"
    σ::NF
    "melting point (K, default `273.15`)"
    T0::NF
    "seconds per day (default `86400`)"
    seconds_per_day::NF
    # Surface radiative properties
    "longwave emissivity of snow (–, default `0.98`)"
    ϵ_snow::NF
    "longwave emissivity of bare ice (–, default `0.98`)"
    eps_ice::NF
    "emissivity of a graybody sky, used by `GraybodyLongwave` (–, default `0.8`)"
    ϵ_air::NF
    "albedo of dry fresh snow (–, default `0.81`)"
    alpha_dry::NF
    "albedo of wet snow (–, default `0.70`)"
    alpha_wet::NF
    "albedo of bare ice (–, default `0.40`)"
    alpha_ice::NF
    "volumetric liquid-water content at which `DynamicAlbedo` reaches `alpha_wet` (–, default `0.1`)"
    max_lwc_albedo::NF
    # Initial state and refreezing
    "density of the initial layers (kg m⁻³, default `300`)"
    density_init::NF
    "temperature of the initial layers and substrate (K, default `273`)"
    temperature_init::NF
    "factor on the refreezing capacity (–, default `1`)"
    refreezing_correction::NF
    # Schemes
    "density of fresh snow, an [`AbstractFreshSnowDensityScheme`](@ref) (default `ConstantFreshSnowDensity()`)"
    fresh_snow_density_scheme::FS
    "snow albedo, an [`AbstractAlbedoScheme`](@ref) (default `DynamicAlbedo()`)"
    albedo_scheme::A
    "densification of low-density snow, an [`AbstractDensificationScheme`](@ref) (default `BESSIDensification()`)"
    densification_scheme::D
    "longwave exchange of the surface, an [`AbstractSurfaceEnergyBalanceScheme`](@ref) (default `SemixSurfaceEnergyBalance()`)"
    surface_energy_balance_scheme::S
    "turbulent heat fluxes, an [`AbstractTurbulentFluxScheme`](@ref) (default `SemixTurbulentFluxes()`)"
    turbulent_flux_scheme::T
    "longwave down where the forcing lacks it, an [`AbstractLongwaveScheme`](@ref) (default `CloudProxyLongwave()`)"
    longwave_scheme::L
    "diurnal cycle within a forcing step, an [`AbstractDiurnalCycle`](@ref) (default `DiurnalSubsteps()`)"
    diurnal_cycle::DC
end

# The default of every field of `Params`, in field order. A function, since the schemes are
# defined after this file.
_param_defaults(::Type{<:Params}) = (
    rho_i = 917.0,
    rho_w = 1000.0,
    ci = 2110.0,
    cw = 4181.0,
    Lm = 334000.0,
    Lv = 2.501e6,
    cp_air = 1003.0,
    σ = 5.670373e-8,
    T0 = 273.15,
    seconds_per_day = DEFAULT_SECONDS_PER_DAY,
    ϵ_snow = 0.98,
    eps_ice = 0.98,
    ϵ_air = 0.8,
    alpha_dry = 0.81,
    alpha_wet = 0.70,
    alpha_ice = 0.40,
    max_lwc_albedo = 0.1,
    density_init = DEFAULT_DENSITY_INIT,
    temperature_init = DEFAULT_TEMPERATURE_INIT,
    refreezing_correction = 1.0,
    fresh_snow_density_scheme = ConstantFreshSnowDensity(),
    albedo_scheme = DynamicAlbedo(),
    densification_scheme = BESSIDensification(),
    surface_energy_balance_scheme = SemixSurfaceEnergyBalance(),
    turbulent_flux_scheme = SemixTurbulentFluxes(),
    longwave_scheme = CloudProxyLongwave(),
    diurnal_cycle = DiurnalSubsteps(),
)

"""
    Params(; NF = Float64, kwargs...)

Keyword constructor of [`Params`](@ref): every keyword is a field of the same name (see there for
units and defaults), and every number and the floating-point fields of the schemes are converted to
`NF`. `Params()` is the Greenland calibration.
"""
function Params(; NF::Type{<:Real} = Float64, kwargs...)
    params = _build_params(Params, _param_defaults(Params), NF, values(kwargs))
    validate(params)
    return params
end

"""`params` with the fields in `kwargs` replaced (schemes and numbers, converted to its `NF`)."""
Params(params::Params; kwargs...) =
    Params(; NF = number_type(params), _fields(params)..., kwargs...)

# The fields of `x` as a NamedTuple; the length is a compile-time `Val`, so the result is
# type-stable for any number of fields (and differentiable by Enzyme).
@inline _fields(x::T) where {T} = NamedTuple{fieldnames(T)}(ntuple(i -> getfield(x, i), Val(fieldcount(T))))

function validate(params::Params)
    _require(isfinite(params.refreezing_correction) && params.refreezing_correction > 0,
        "`refreezing_correction` must be finite and positive.")
    for name in fieldnames(Params)
        value = getfield(params, name)
        value isa Number || validate(value)
    end
    _validate_albedo(params.albedo_scheme, params)
    return nothing
end

_validate_albedo(::AbstractAlbedoScheme, params) = nothing
_validate_albedo(::AgingAlbedo, params) = _require(0 <= params.alpha_wet <= params.alpha_dry <= 1,
    "AgingAlbedo needs 0 ≤ alpha_wet ≤ alpha_dry ≤ 1.")

number_type(::Params{NF}) where {NF} = NF
Base.eltype(::Params{NF}) where {NF} = NF

# Only scalars and schemes: the parameters pass to GPU kernels unchanged.
Adapt.adapt_structure(to, params::Params) = params

"""
    ParamsGrISv1(; NF = Float64, kwargs...)

The Greenland calibration against MAR v3.14.3 (daily ERA5 forcing, 1980–2019), the defaults of
[`Params`](@ref). `kwargs` override single fields.
"""
ParamsGrISv1(; kwargs...) = Params(; kwargs...)

"""
    ParamsAISv1(; NF = Float64, kwargs...)

The Antarctic calibration against MAR (2010–2019): the Greenland calibration with a refitted
cloud-proxy longwave (stage 3 of `fit_cloud_proxy_longwave_antarctica.jl`), weaker stable damping
of the SEMIX turbulent exchange with a larger latent factor, and brighter snow (`alpha_dry`,
`alpha_wet`). `kwargs` override single fields.
"""
function ParamsAISv1(; kwargs...)
    return Params(;
        longwave_scheme = CloudProxyLongwave(
            clear_sky_emissivity = 0.7314, clear_sky_emissivity_temperature_slope = 0.003378,
            overcast_emissivity = 0.9877, overcast_emissivity_temperature_slope = 0.004171,
            cloud_emissivity_exponent = 2.0, night_cloud_fraction = 0.9),
        turbulent_flux_scheme = SemixTurbulentFluxes(
            sensible_exchange_factor = 2.5, stable_coefficient = 20.0, latent_exchange_factor = 3.0),
        alpha_dry = 0.87, alpha_wet = 0.82,
        kwargs...)
end

# ---------------------------------------------------------------------------
# PDD
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

Parameters of a [`PDDModel`](@ref): degree-day factors of snow and ice, the refrozen fraction of
meltwater, the temperature spread of the positive degree days, and the snow cap above which snow
turns into ice.

# Fields
$(TYPEDFIELDS)
"""
struct PDDParams{NF <: Real, M <: AbstractPDDMethod}
    "degree-day factor of snow (mm w.e. K⁻¹ day⁻¹, default `3`)"
    ddf_snow::NF
    "degree-day factor of ice (mm w.e. K⁻¹ day⁻¹, default `8`)"
    ddf_ice::NF
    "fraction of melt and rain that refreezes (–, default `0.6`)"
    refreezing_fraction::NF
    "standard deviation of the air temperature (K, default `5`)"
    temperature_sigma::NF
    "snow water equivalent above which snow becomes ice (mm w.e., default `5000`)"
    H_snow_max::NF
    "melting point (K, default `273.15`)"
    T0::NF
    "seconds per day (default `86400`)"
    seconds_per_day::NF
    "positive degree days, an [`AbstractPDDMethod`](@ref) (default `SimplePDD()`)"
    method::M
end

_param_defaults(::Type{<:PDDParams}) = (ddf_snow = 3.0, ddf_ice = 8.0, refreezing_fraction = 0.6,
    temperature_sigma = 5.0, H_snow_max = 5000.0, T0 = 273.15, seconds_per_day = DEFAULT_SECONDS_PER_DAY,
    method = SimplePDD())

"""
    PDDParams(; NF = Float64, kwargs...)

Keyword constructor of [`PDDParams`](@ref); `method = PISMPDD()` uses the Calov–Greve expectation
integral at every step.
"""
function PDDParams(; NF::Type{<:Real} = Float64, kwargs...)
    params = _build_params(PDDParams, _param_defaults(PDDParams), NF, values(kwargs))
    _require(params.ddf_snow > 0, "`ddf_snow` must be positive.")
    _require(params.ddf_ice >= 0, "`ddf_ice` must be non-negative.")
    _require(0 <= params.refreezing_fraction <= 1, "`refreezing_fraction` must be between 0 and 1.")
    _require(params.temperature_sigma > 0, "`temperature_sigma` must be positive.")
    _require(params.H_snow_max > 0, "`H_snow_max` must be positive.")
    return params
end

number_type(::PDDParams{NF}) where {NF} = NF
Adapt.adapt_structure(to, params::PDDParams) = params

# ---------------------------------------------------------------------------
# ITM
# ---------------------------------------------------------------------------

"""
$(TYPEDEF)

Parameters of an [`ITMModel`](@ref), the insolation–temperature–melt scheme of the Fortran
`fesmc/chion`/smbpal code: melt M = max(τ(z)(1 − α) S + c + t (T − T₀), 0) / (ρ_w L_m) with the
transmissivity τ = a + b √z and a latitude-dependent offset c.

# Fields
$(TYPEDFIELDS)
"""
struct ITMParams{NF <: Real}
    "transmissivity at sea level, a (–, default `0.46`)"
    trans_a::NF
    "transmissivity increase with √height, b (m^-½, default `6e-5`)"
    trans_b::NF
    "melt offset c (W m⁻², default `-45`)"
    itm_c::NF
    "temperature coefficient t (W m⁻² K⁻¹, default `10`)"
    itm_t::NF
    "latitude slope of the offset (W m⁻² degree⁻¹, default `-2`)"
    itm_b::NF
    "reference latitude of the offset (degree, default `65`; 90 or more disables the slope)"
    itm_lat0::NF
    "snow water equivalent above which snow becomes ice (mm w.e., default `5000`)"
    H_snow_max::NF
    "maximum refrozen fraction of the precipitation (–, default `0.6`)"
    Pmaxfrac::NF
    "snow depth of full snow cover on desert (mm w.e., default `10`)"
    H_snow_crit_desert::NF
    "snow depth of full snow cover on forest (mm w.e., default `100`)"
    H_snow_crit_forest::NF
    "melt rate above which snow counts as wet (mm w.e. day⁻¹, default `0.5`)"
    melt_crit::NF
    "albedo of the ocean (–, default `0.1`)"
    alb_ocean::NF
    "albedo of bare land (–, default `0.2`)"
    alb_land::NF
    "albedo of forest (–, default `0.1`)"
    alb_forest::NF
    "albedo of bare ice (–, default `0.4`)"
    alb_ice::NF
    "albedo of dry snow (–, default `0.8`)"
    alb_snow_dry::NF
    "albedo of wet snow (–, default `0.65`)"
    alb_snow_wet::NF
    "warming of the firn by refreezing (K per mm w.e. day⁻¹, default `0.0266`)"
    firn_fac::NF
    "liquid water density (kg m⁻³, default `1000`)"
    rho_w::NF
    "latent heat of fusion (J kg⁻¹, default `3.34e5`)"
    Lm::NF
    "melting point (K, default `273.15`)"
    T0::NF
    "seconds per day (default `86400`)"
    seconds_per_day::NF
end

_param_defaults(::Type{<:ITMParams}) = (trans_a = 0.46, trans_b = 6e-5, itm_c = -45.0, itm_t = 10.0,
    itm_b = -2.0, itm_lat0 = 65.0, H_snow_max = 5000.0, Pmaxfrac = 0.6, H_snow_crit_desert = 10.0,
    H_snow_crit_forest = 100.0, melt_crit = 0.5, alb_ocean = 0.1, alb_land = 0.2, alb_forest = 0.1,
    alb_ice = 0.4, alb_snow_dry = 0.8, alb_snow_wet = 0.65, firn_fac = 0.0266, rho_w = 1000.0,
    Lm = 334000.0, T0 = 273.15, seconds_per_day = DEFAULT_SECONDS_PER_DAY)

"""
    ITMParams(; NF = Float64, kwargs...)

Keyword constructor of [`ITMParams`](@ref).
"""
function ITMParams(; NF::Type{<:Real} = Float64, kwargs...)
    params = _build_params(ITMParams, _param_defaults(ITMParams), NF, values(kwargs))
    _require(params.H_snow_max > 0, "`H_snow_max` must be positive.")
    _require(params.H_snow_crit_desert > 0, "`H_snow_crit_desert` must be positive.")
    _require(params.H_snow_crit_forest > 0, "`H_snow_crit_forest` must be positive.")
    _require(0 <= params.Pmaxfrac <= 1, "`Pmaxfrac` must be between 0 and 1.")
    _require(all(name -> 0 <= getfield(params, name) <= 1,
        (:alb_ocean, :alb_land, :alb_forest, :alb_ice, :alb_snow_dry, :alb_snow_wet)), "ITM albedos must be between 0 and 1.")
    return params
end

number_type(::ITMParams{NF}) where {NF} = NF
Adapt.adapt_structure(to, params::ITMParams) = params

# ---------------------------------------------------------------------------
# Display
# ---------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", params::Union{Params, PDDParams, ITMParams})
    P = typeof(params)
    defaults = _build_params(P, _param_defaults(P), number_type(params), (;))
    println(io, nameof(P), "{", number_type(params), "}")
    for name in fieldnames(P)
        value = getfield(params, name)
        value isa Number && isequal(value, getfield(defaults, name)) && continue
        println(io, "  ", rpad(String(name), 30), value)
    end
    print(io, "  (numbers not listed are at their defaults)")
    return nothing
end
