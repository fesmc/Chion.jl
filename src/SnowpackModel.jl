"""
Array-first snowpack model with Terrarium-style state containers.
"""

module SnowpackModel

using Printf
using Adapt: Adapt, adapt, @adapt_structure
using CUDA
using Enzyme
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

include("domain_interfaces.jl")
include("execution_backends.jl")
include("snowpack_domain.jl")
include("step_forcing.jl")
include("step_scratch.jl")
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
include("domain_summaries.jl")
include("step_timing.jl")
include("step.jl")
include("step_field_batches.jl")
include("state_access.jl")
include("enzyme_ad.jl")

end # module
