module Chion

using Printf
using Dates
using NCDatasets
using Adapt: Adapt, adapt, @adapt_structure
using KernelAbstractions
using KernelAbstractions: CPU
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using SpecialFunctions: erfc
using TOML

# ---------------------------------------------------------------------------
# Numerics, schemes, parameters and grid
# ---------------------------------------------------------------------------
include("backends.jl")
include("constants.jl")
include("numerics.jl")
include("schemes.jl")
include("boundary_conditions.jl")
include("params.jl")
include("grid.jl")

# ---------------------------------------------------------------------------
# Forcing
# ---------------------------------------------------------------------------
include("forcing.jl")
include("dataloaders.jl")
include("forcing_source.jl")

# ---------------------------------------------------------------------------
# State and column physics
# ---------------------------------------------------------------------------
include("column_state_utils.jl")
include("state.jl")
include("processes/albedo.jl")
include("processes/layer_structure.jl")
include("processes/accumulation.jl")
include("processes/melt.jl")
include("processes/diurnal_shortwave.jl")
include("processes/surface_fluxes.jl")
include("processes/energy_flux.jl")
include("processes/densification.jl")
include("processes/percolation.jl")
include("processes/refreezing.jl")
include("processes/pdd.jl")
include("processes/itm.jl")
include("step.jl")

# ---------------------------------------------------------------------------
# Model, outputs, simulation
# ---------------------------------------------------------------------------
include("model.jl")
include("diagnostics.jl")
include("io.jl")
include("outputs.jl")
include("run_output.jl")
include("simulation.jl")
include("config.jl")

# ---------------------------------------------------------------------------
# Exports — public API
# ---------------------------------------------------------------------------

# Grid and parameters
export Grid, Layering, SpatialLayout
export Params, PDDParams, ITMParams, ParamsGrISv1, ParamsAISv1

# Schemes
export ConstantFreshSnowDensity, ParameterizedFreshSnowDensity
export ConstantAlbedo, DynamicAlbedo, PrescribedAlbedo, AgingAlbedo, SemixAlbedo
export DangSnowAlbedo, WarrenWiscombeSnowAlbedo
export BESSIDensification, HTESSELDensification
export SemixSurfaceEnergyBalance, BESSISurfaceEnergyBalance
export SemixTurbulentFluxes, BESSITurbulentFluxes
export CloudProxyLongwave, GraybodyLongwave
export NoDiurnalCycle, DiurnalSubsteps, DiurnalTemperatureCycle, NoDiurnalTemperatureCycle
export SimplePDD, PISMPDD

# Boundary conditions
export BoundaryConditions, SurfaceEnergyBalance, PrescribedSurfaceFlux, PrescribedSurfaceTemperature
export Insulated, GeothermalFlux, PrescribedBottomTemperature

# Models
export BESSIModel, PDDModel, ITMModel
export BESSIState, PDDState, ITMState, initial_state
export CPU, to_backend, cpu_state, set_state!, set_active_mask!, set_input!, sync_forcing!, update_forcing!
export integrate!, update_diagnostics!, column_enthalpy, column_profile

# Forcing
export SnowpackForcing, Forcing
export forcing_grid, forcing_area, forcing_years, load_forcing_year, slice_forcing, relabel_year
export load_forcing_file, read_forcing_columns
export clear_sky_transmissivity
export shortwave_down_from_transmissivity, wind_speed_from_geometry
export temperature_variability_sigma, red_noise_monthly, update_air_pressure!

# Simulation and output
export Simulation, Clock, run!, time_step!, spin_up!, finished, yearly_step!, smbstats
export NoSpinup, CycledSpinup, ClimatologySpinup
export Output, OutputConfig, available_output_variables, run_config, write_restart
export Hour, Day, Week, Month, Year     # output frequencies, from Dates

end
