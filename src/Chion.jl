module Chion

using Printf
using Dates
using Serialization
using NCDatasets
using Adapt: Adapt, adapt, @adapt_structure
using CUDA
using KernelAbstractions
using ProgressMeter: Progress, next!, update!

# ---------------------------------------------------------------------------
# Core domain, forcing, and physics
# ---------------------------------------------------------------------------
include("constants.jl")
include("domain.jl")
include("forcing.jl")
include("timing.jl")
include("step.jl")
include("column_state_utils.jl")
include("solvers/batched_tridiagonal_solver.jl")
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
include("dataloaders.jl")
include("io.jl")
include("convenience_types.jl")
include("reporting.jl")
include("models.jl")
include("state.jl")
include("processes/pdd.jl")
include("run_types.jl")
include("runtime.jl")
include("integrators.jl")
include("simulation.jl")

# ---------------------------------------------------------------------------
# Exports — public API
# ---------------------------------------------------------------------------

# Grid
export SnowpackGrid

# Models
export BESSIModel, PDDModel, StochasticMonthlyPDD
export AbstractState, AbstractSnowModelState, ReferenceState, CurrentState, PDDState
export build_model, initial_state
export DynamicAlbedo, ConstantAlbedo, PrescribedAlbedo
export BESSIDensification, HTESSELDensification
export ConstantFreshSnowDensity, ParameterizedFreshSnowDensity

# Forcing
export SnowpackForcing

# Simulation
export Simulation, SimulationResult, RunOptions
export SimulationIntegrator
export NetcdfOutput, MonthlyState, io_dict, write_nc!
export init_integrator, step!, yearly_step!, run!, finalize!, finished, set_active_mask!, state
export update_air_pressure!, sync_forcing!

# Data loading
export load_forcing_file

# Low-level / power-user exports
export step!
export step_interval_threads!
export ColumnarStepWorkspace
export update_diagnostics!
export cpu_state, gpu_state
export continuous_bottom_deplete!
export update_surface_albedo!
export get_state, print_state
export add_timing!, timing_rows, print_timing_summary
export StepTimingStats

end
