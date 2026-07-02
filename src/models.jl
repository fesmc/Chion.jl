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

"""
    PDDModel(grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.6,
             temperature_sigma=5.0)

Bulk positive-degree-day surface mass-balance model configuration. Evolving
state is stored in `PDDState` and owned by `Simulation.now`. Monthly forcing
steps use a PISM-style expectation integral with normally-distributed
unresolved temperature variability.
"""
struct PDDModel{G <: SnowpackGrid}
    grid::G
    ddf_snow::Float64
    ddf_ice::Float64
    refreezing_fraction::Float64
    temperature_sigma::Float64
end

function PDDModel(
    grid::SnowpackGrid;
    ddf_snow::Real=3.0,
    ddf_ice::Real=8.0,
    refreezing_fraction::Real=0.6,
    temperature_sigma::Real=5.0,
)
    ddf_snow > 0 || error("`ddf_snow` must be positive.")
    ddf_ice >= 0 || error("`ddf_ice` must be non-negative.")
    0.0 <= refreezing_fraction <= 1.0 || error("`refreezing_fraction` must be between 0 and 1.")
    temperature_sigma > 0 || error("`temperature_sigma` must be positive.")

    return PDDModel(
        grid,
        Float64(ddf_snow),
        Float64(ddf_ice),
        Float64(refreezing_fraction),
        Float64(temperature_sigma),
    )
end
