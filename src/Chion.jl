module Chion

using Printf
using Dates
using NCDatasets
using Base.Threads: @threads
using Adapt: Adapt, adapt, @adapt_structure
using CUDA
using KernelAbstractions
using ProgressMeter: Progress, next!

# ---------------------------------------------------------------------------
# Core domain, forcing, and physics
# ---------------------------------------------------------------------------
include("constants.jl")
include("domain.jl")
include("forcing.jl")
include("timing.jl")
include("step.jl")
include("column_state_utils.jl")
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
include("integrators.jl")
include("runtime.jl")
include("diagnostics.jl")
include("output_runtime.jl")
include("simulation.jl")

# ---------------------------------------------------------------------------
# Exports — public API
# ---------------------------------------------------------------------------

# Grid
export SnowpackGrid

# Models
export BESSIModel, PDDModel, ITMModel
export AbstractSnowModelState, BESSIState, PDDState, ITMState
export build_model, initial_state
export DynamicAlbedo, ConstantAlbedo
export BESSIDensification, HTESSELDensification
export ConstantFreshSnowDensity, ParameterizedFreshSnowDensity

# Forcing
export SnowpackForcing

# Simulation
export Simulation, SimulationResult, SimulationOptions, OutputOptions
export SimulationIntegrator
export init_integrator, step!, run!, finalize!, finished, set_forcing!, state

# Data loading
export load_forcing_file

# Low-level / power-user exports
export step!
export ColumnarStepWorkspace
export continuous_bottom_deplete!
export update_surface_albedo!
export get_state, print_state
export add_timing!, timing_rows, print_timing_summary
export StepTimingStats

end
