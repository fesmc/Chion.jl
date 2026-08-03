"""
Model definitions exposed by Chion's simulation-first API.
"""

"""
    BESSIModel(grid; albedo=:dynamic, densification=:bessi, ...)

Configuration for the layered BESSI snowpack model. Evolving state is stored in
`BESSIState` and owned by `Simulation.now`.
"""
struct BESSIModel{G <: SnowpackGrid, C <: SnowpackPhysicalConstants}
    grid::G
    c::C
    Ntot::Int
    mass_max::Float64
    mass_split::Float64
    mass_min::Float64
    density_init::Float64
    temperature_init::Float64
    diurnal_shortwave_substeps::Bool
    diurnal_shortwave_threshold::Float64
    diurnal_shortwave_max_substeps::Int
    diurnal_shortwave_min_air_temperature::Float64
    diurnal_temperature_cycle::Bool
    diurnal_temperature_amplitude::Float64
end

function BESSIModel(
    grid::SnowpackGrid;
    albedo::Symbol=:dynamic,
    densification::Symbol=:bessi,
    fresh_snow_density::Symbol=:constant,
    Ntot::Int=DEFAULT_NTOT,
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    density_init::Real=DEFAULT_DENSITY_INIT,
    temperature_init::Real=DEFAULT_TEMPERATURE_INIT,
    diurnal_shortwave::Bool=false,
    diurnal_shortwave_substeps::Bool=false,
    diurnal_shortwave_threshold::Real=0.0,
    diurnal_shortwave_max_substeps::Integer=3,
    diurnal_shortwave_min_air_temperature_c::Real=-8.0,
    diurnal_temperature_cycle::Bool=false,
    diurnal_temperature_amplitude_c::Real=5.0,
    kwargs...,
)
    diurnal_shortwave_threshold >= 0 || error("`diurnal_shortwave_threshold` must be non-negative.")
    1 <= diurnal_shortwave_max_substeps <= 24 || error("`diurnal_shortwave_max_substeps` must be between 1 and 24.")
    isfinite(diurnal_shortwave_min_air_temperature_c) || error("`diurnal_shortwave_min_air_temperature_c` must be finite.")
    diurnal_temperature_amplitude_c >= 0 || error("`diurnal_temperature_amplitude_c` must be non-negative.")
    resolved_diurnal_shortwave_substeps = diurnal_shortwave_substeps || diurnal_shortwave
    c = SnowpackPhysicalConstants(
        Float64;
        albedo_scheme=albedo,
        low_density_densification=densification,
        fresh_snow_density_scheme=fresh_snow_density,
        kwargs...,
    )
    mass_max, mass_split, mass_min =
        _domain_thresholds(Float64, mass_max, mass_split, mass_min)
    return BESSIModel(
        grid,
        c,
        Ntot,
        mass_max,
        mass_split,
        mass_min,
        Float64(density_init),
        Float64(temperature_init),
        resolved_diurnal_shortwave_substeps,
        Float64(diurnal_shortwave_threshold),
        Int(diurnal_shortwave_max_substeps),
        Float64(diurnal_shortwave_min_air_temperature_c) + 273.15,
        diurnal_temperature_cycle,
        Float64(diurnal_temperature_amplitude_c),
    )
end

const PDD_METHOD_SIMPLE = UInt8(1)
const PDD_METHOD_PISM = UInt8(2)

@inline function _normalize_pdd_method(method::Symbol)
    method == :simple && return PDD_METHOD_SIMPLE
    method in (:pism, :calov_greve) && return PDD_METHOD_PISM
    error("Unsupported PDD method '$method'. Use :simple or :pism.")
end

"""
    PDDModel(grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.6,
             temperature_sigma=5.0, H_snow_max=5000.0, pdd_method=:simple)

Bulk positive-degree-day surface mass-balance model configuration. The model
uses a capped one-layer snow reservoir: refrozen water and snow above
`H_snow_max` become ice, while `smb_ice` records only ice-facing mass changes.
Set `pdd_method=:pism` to use the Calov-Greve expectation integral for every
timestep; `:simple` uses positive mean temperature directly.
"""
struct PDDModel{G <: SnowpackGrid, C <: SnowpackPhysicalConstants}
    grid::G
    c::C
    ddf_snow::Float64
    ddf_ice::Float64
    refreezing_fraction::Float64
    temperature_sigma::Float64
    H_snow_max::Float64
    pdd_method::UInt8
end

"""
    ITMModel(grid; kwargs...)

Insolation-temperature-melt bulk surface-mass-balance model, ported from the
Fortran `fesmc/chion`/smbpal implementation. ITM requires explicit
`surface_height`, `ice_thickness`, `annual_pdd`, and `latitude_deg` fields in
`SnowpackForcing`; `q_sw_net` is used as its insolation driver when supplied,
otherwise `shortwave_down` is used.
"""
struct ITMModel{G <: SnowpackGrid, C <: SnowpackPhysicalConstants}
    grid::G
    c::C
    trans_a::Float64
    trans_b::Float64
    trans_c::Float64
    itm_c::Float64
    itm_t::Float64
    itm_b::Float64
    itm_lat0::Float64
    H_snow_max::Float64
    Pmaxfrac::Float64
    H_snow_crit_desert::Float64
    H_snow_crit_forest::Float64
    melt_crit::Float64
    alb_ocean::Float64
    alb_land::Float64
    alb_forest::Float64
    alb_ice::Float64
    alb_snow_dry::Float64
    alb_snow_wet::Float64
    firn_fac::Float64
end

function ITMModel(
    grid::SnowpackGrid;
    trans_a::Real=0.46,
    trans_b::Real=6e-5,
    trans_c::Real=0.01,
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
    kwargs...,
)
    H_snow_max > 0 || error("`H_snow_max` must be positive.")
    H_snow_crit_desert > 0 || error("`H_snow_crit_desert` must be positive.")
    H_snow_crit_forest > 0 || error("`H_snow_crit_forest` must be positive.")
    0 <= Pmaxfrac <= 1 || error("`Pmaxfrac` must be between 0 and 1.")
    all(value -> 0 <= value <= 1, (alb_ocean, alb_land, alb_forest, alb_ice, alb_snow_dry, alb_snow_wet)) ||
        error("ITM albedos must be between 0 and 1.")
    return ITMModel(grid, SnowpackPhysicalConstants(Float64; kwargs...), Float64(trans_a),
        Float64(trans_b), Float64(trans_c), Float64(itm_c), Float64(itm_t), Float64(itm_b),
        Float64(itm_lat0), Float64(H_snow_max), Float64(Pmaxfrac),
        Float64(H_snow_crit_desert), Float64(H_snow_crit_forest), Float64(melt_crit),
        Float64(alb_ocean), Float64(alb_land), Float64(alb_forest), Float64(alb_ice),
        Float64(alb_snow_dry), Float64(alb_snow_wet), Float64(firn_fac))
end

function PDDModel(
    grid::SnowpackGrid;
    ddf_snow::Real=3.0,
    ddf_ice::Real=8.0,
    refreezing_fraction::Real=0.6,
    temperature_sigma::Real=5.0,
    H_snow_max::Real=5000.0,
    pdd_method::Symbol=:simple,
    kwargs...,
)
    ddf_snow > 0 || error("`ddf_snow` must be positive.")
    ddf_ice >= 0 || error("`ddf_ice` must be non-negative.")
    0.0 <= refreezing_fraction <= 1.0 || error("`refreezing_fraction` must be between 0 and 1.")
    temperature_sigma > 0 || error("`temperature_sigma` must be positive.")
    H_snow_max > 0 || error("`H_snow_max` must be positive.")

    return PDDModel(
        grid,
        SnowpackPhysicalConstants(Float64; kwargs...),
        Float64(ddf_snow),
        Float64(ddf_ice),
        Float64(refreezing_fraction),
        Float64(temperature_sigma),
        Float64(H_snow_max),
        _normalize_pdd_method(pdd_method),
    )
end
