"""
Array-first snowpack model with Terrarium-style state containers.
"""

module SnowpackModel

using Printf
using Adapt: Adapt, adapt, @adapt_structure
using CUDA
using KernelAbstractions

include("model_constants.jl")

export SnowpackPhysicalConstants
export SnowpackStepForcing
export SnowpackStepFields
export AbstractSnowpackDomain
export SnowpackDomain
export StepWorkspace
export threaded_workspaces
export ColumnarStepWorkspace
export step!
export go_percolation!
export go_refreezing!
export continuous_bottom_deplete!
export get_state
export print_state
export go_densification!
export update_surface_albedo!
export StepTimingStats
export add_timing!
export timing_rows
export print_timing_summary
export column_count
export cpu_domain
export gpu_domain
export cuda_available
export kernelabstractions_available
export summarize_domain_state
export summarize_domain_state!
export variables
export compute_auxiliary!
export compute_tendencies!

include("abstractions.jl")
include("backend_utils.jl")
include("snowpack_types.jl")
include("step_forcing_types.jl")
include("step_workspaces.jl")
include("column_helpers.jl")
include("albedo.jl")
include("mass_balance.jl")
include("energy_flux.jl")
include("densification.jl")
include("percolation.jl")
include("refreezing.jl")
include("diurnal_shortwave.jl")
include("batch.jl")
include("timing.jl")
include("step_process_helpers.jl")
include("step.jl")
include("step_batch.jl")
include("state_io.jl")

end # module
