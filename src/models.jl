"""
Model definitions exposed by Chion's simulation-first API.
"""

"""Dynamic surface albedo scheme for `BESSIModel`."""
struct DynamicAlbedo end
"""Constant surface albedo scheme for `BESSIModel`."""
struct ConstantAlbedo end
"""Use prescribed surface albedo from `SnowpackForcing`."""
struct PrescribedAlbedo end
"""BESSI low-density densification scheme."""
struct BESSIDensification end
"""HTESSEL-style low-density densification scheme."""
struct HTESSELDensification end
"""Use the constant fresh-snow density parameter."""
struct ConstantFreshSnowDensity end
"""Use the temperature/wind-dependent fresh-snow density parameterization."""
struct ParameterizedFreshSnowDensity end

_albedo_symbol(::DynamicAlbedo) = :dynamic
_albedo_symbol(::ConstantAlbedo) = :constant
_albedo_symbol(::PrescribedAlbedo) = :prescribed
_densification_symbol(::BESSIDensification) = :bessi
_densification_symbol(::HTESSELDensification) = :htessel
_fresh_snow_symbol(::ConstantFreshSnowDensity) = :constant
_fresh_snow_symbol(::ParameterizedFreshSnowDensity) = :parameterized

"""
    BESSIModel(grid; albedo=DynamicAlbedo(), densification=BESSIDensification(), ...)

Configuration for the layered BESSI snowpack model. Evolving state is stored in
`CurrentState` and owned by `Simulation.now`.
"""
struct BESSIModel{G <: AbstractSnowpackGrid, C <: SnowpackPhysicalConstants} <: AbstractSnowModel{G}
    grid::G
    c::C
    Ntot::Int
    mass_max::Float64
    mass_split::Float64
    mass_min::Float64
    rho_max::Float64
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
    grid::AbstractSnowpackGrid;
    albedo=DynamicAlbedo(),
    densification=BESSIDensification(),
    fresh_snow_density=ConstantFreshSnowDensity(),
    Ntot::Int=DEFAULT_NTOT,
    mass_max::Real=DEFAULT_MASS_MAX,
    mass_split::Real=DEFAULT_MASS_SPLIT,
    mass_min::Real=DEFAULT_MASS_MIN,
    rho_max::Real=DEFAULT_RHO_MAX,
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
        albedo_scheme=_albedo_symbol(albedo),
        low_density_densification=_densification_symbol(densification),
        fresh_snow_density_scheme=_fresh_snow_symbol(fresh_snow_density),
        kwargs...,
    )
    return BESSIModel(
        grid,
        c,
        Ntot,
        Float64(mass_max),
        Float64(mass_split),
        Float64(mass_min),
        Float64(rho_max),
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

function SnowpackDomain(model::BESSIModel)
    grid = model.grid
    return SnowpackDomain(;
        c=model.c,
        Ntot=model.Ntot,
        ncol=ncols(grid),
        mass_max=model.mass_max,
        mass_split=model.mass_split,
        mass_min=model.mass_min,
        rho_max=model.rho_max,
        x=grid.x,
        y=grid.y,
        js=grid.js,
        is=grid.is,
        mask=grid.mask,
    )
end

model_domain(model::BESSIModel) = SnowpackDomain(model)
model_domain(model::AbstractSnowModel) = model.grid

"""PISM-style expectation-integral monthly PDD parameterization."""
struct StochasticMonthlyPDD
    temperature_sigma::Float64
end

function StochasticMonthlyPDD(; temperature_sigma::Real=5.0)
    temperature_sigma > 0 || error("`temperature_sigma` must be positive.")
    return StochasticMonthlyPDD(Float64(temperature_sigma))
end

"""
    PDDModel(grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.6,
             monthly_method=StochasticMonthlyPDD())

Bulk positive-degree-day surface mass-balance model configuration. Evolving
state is stored in `PDDState` and owned by `Simulation.now`. Monthly forcing
steps use a PISM-style expectation integral with normally-distributed
unresolved temperature variability.
"""
struct PDDModel{G <: AbstractSnowpackGrid} <: AbstractSnowModel{G}
    grid::G
    ddf_snow::Float64
    ddf_ice::Float64
    refreezing_fraction::Float64
    monthly_method::StochasticMonthlyPDD
end

function PDDModel(
    grid::AbstractSnowpackGrid;
    ddf_snow::Real=3.0,
    ddf_ice::Real=8.0,
    refreezing_fraction::Real=0.6,
    monthly_method::StochasticMonthlyPDD=StochasticMonthlyPDD(),
)
    ddf_snow > 0 || error("`ddf_snow` must be positive.")
    ddf_ice >= 0 || error("`ddf_ice` must be non-negative.")
    0.0 <= refreezing_fraction <= 1.0 || error("`refreezing_fraction` must be between 0 and 1.")

    return PDDModel(
        grid,
        Float64(ddf_snow),
        Float64(ddf_ice),
        Float64(refreezing_fraction),
        monthly_method,
    )
end


"""
    ITMModel(grid; c_rad=0.513, c_temp=0.362, T_melt=0.0, refreezing_fraction=0.6)

Index-temperature model placeholder. Construction is supported; `run!` is not
implemented yet.
"""
struct ITMModel{G <: AbstractSnowpackGrid} <: AbstractSnowModel{G}
    grid::G
    c_rad::Float64
    c_temp::Float64
    T_melt::Float64
    refreezing_fraction::Float64
end

function ITMModel(
    grid::AbstractSnowpackGrid;
    c_rad::Real=0.513,
    c_temp::Real=0.362,
    T_melt::Real=0.0,
    refreezing_fraction::Real=0.6,
)
    return ITMModel(grid, Float64(c_rad), Float64(c_temp), Float64(T_melt), Float64(refreezing_fraction))
end
