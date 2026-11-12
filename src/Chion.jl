module Chion

using Printf
using Dates
using NCDatasets
using Adapt: Adapt, adapt, @adapt_structure
using KernelAbstractions
using ProgressMeter: Progress, next!, update!
using SpecialFunctions: erfc

# ---------------------------------------------------------------------------
# Core domain, forcing, and physics
# ---------------------------------------------------------------------------
include("backends.jl")
include("constants.jl")
include("domain.jl")
include("forcing.jl")
include("timing.jl")
include("column_state_utils.jl")
include("models.jl")
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
include("diagnostics.jl")
include("dataloaders.jl")
include("io.jl")
include("monthly_output.jl")
include("run_types.jl")
include("reporting.jl")
include("runtime.jl")
include("integrators.jl")
include("simulation.jl")
include("outputs.jl")
include("forcing_source.jl")
include("experiment.jl")
include("config.jl")

# ---------------------------------------------------------------------------
# Exports — public API
# ---------------------------------------------------------------------------

# Grid
export SnowpackGrid

# Models
export BESSIModel, PDDModel, ITMModel
export BESSIState, PDDState, ITMState
export build_model, initial_state

# Forcing
export SnowpackForcing

# Simulation
export Simulation, SimulationResult, RunOptions
export SimulationIntegrator
export init_integrator, step!, yearly_step!, run!, finalize!, finished, set_active_mask!
export update_air_pressure!, sync_forcing!

# Experiments: forcing sources, spin-up, outputs
export Forcing, Output, Experiment, ExperimentResult
export NoSpinup, CycledSpinup, ClimatologySpinup
export forcing_grid, forcing_area, forcing_years, load_forcing_year, slice_forcing, relabel_year
export available_output_variables, run_config
export Hour, Day, Week, Month, Year     # output frequencies, from Dates

# Data loading
export load_forcing_file, read_forcing_columns

# Low-level / power-user exports
export update_diagnostics!
export cpu_state, gpu_state
export get_state, print_state
export add_timing!, timing_rows, print_timing_summary
export StepTimingStats

end
