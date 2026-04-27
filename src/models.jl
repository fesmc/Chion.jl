"""
Model definitions exposed by Chion's simulation-first API.
"""

"""Dynamic surface albedo scheme for `BESSIModel`."""
struct DynamicAlbedo end
"""Constant surface albedo scheme for `BESSIModel`."""
struct ConstantAlbedo end
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
_densification_symbol(::BESSIDensification) = :bessi
_densification_symbol(::HTESSELDensification) = :htessel
_fresh_snow_symbol(::ConstantFreshSnowDensity) = :constant
_fresh_snow_symbol(::ParameterizedFreshSnowDensity) = :parameterized

"""
    BESSIModel(grid; albedo=DynamicAlbedo(), densification=BESSIDensification(), ...)

Layered snowpack model backed by the BESSI process implementation.
"""
struct BESSIModel{G <: AbstractSnowpackGrid} <: AbstractSnowModel{G}
    grid::G
    domain::SnowpackDomain
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
    density_init::Real=DEFAULT_DENSITY_INIT,
    temperature_init::Real=DEFAULT_TEMPERATURE_INIT,
    kwargs...,
)
    c = SnowpackPhysicalConstants(
        Float64;
        albedo_scheme=_albedo_symbol(albedo),
        low_density_densification=_densification_symbol(densification),
        fresh_snow_density_scheme=_fresh_snow_symbol(fresh_snow_density),
        kwargs...,
    )
    domain = SnowpackDomain(;
        c=c,
        Ntot=Ntot,
        ncol=ncols(grid),
        mass_max=mass_max,
        mass_split=mass_split,
        mass_min=mass_min,
        density_init=density_init,
        temperature_init=temperature_init,
    )
    return BESSIModel(grid, domain)
end

"""
    PDDModel(grid; ddf_snow=3.0, ddf_ice=8.0, refreezing_fraction=0.6)

Bulk positive-degree-day surface mass-balance model. It stores per-column
snowpack water equivalent, cumulative ice-sheet SMB, and runoff diagnostics.
"""
mutable struct PDDModel{G <: AbstractSnowpackGrid} <: AbstractSnowModel{G}
    grid::G
    ddf_snow::Float64
    ddf_ice::Float64
    refreezing_fraction::Float64
    snowpack_swe::Vector{Float64}
    smb_ice::Vector{Float64}
    runoff::Vector{Float64}
    pdd_sum::Vector{Float64}
end

function PDDModel(
    grid::AbstractSnowpackGrid;
    ddf_snow::Real=3.0,
    ddf_ice::Real=8.0,
    refreezing_fraction::Real=0.6,
)
    ddf_snow > 0 || error("`ddf_snow` must be positive.")
    ddf_ice >= 0 || error("`ddf_ice` must be non-negative.")
    0.0 <= refreezing_fraction <= 1.0 || error("`refreezing_fraction` must be between 0 and 1.")
    ncol = ncols(grid)
    return PDDModel(
        grid,
        Float64(ddf_snow),
        Float64(ddf_ice),
        Float64(refreezing_fraction),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
        zeros(Float64, ncol),
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
