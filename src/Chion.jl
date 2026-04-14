module Chion

using Printf
using Dates
using HDF5
using Base.Threads: @threads
using Adapt: Adapt, adapt, @adapt_structure
using CUDA
using Enzyme
using KernelAbstractions

include("model_constants.jl")
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
include("case_runtime.jl")
include("cases/definitions.jl")
include("cases/forcing_file_definition.jl")
include("cases/api.jl")

export step!
export SnowpackDomain
export SnowpackStepForcing
export SnowpackStepFields
export SnowpackPhysicalConstants
export RunConfig, RunResult
export SnowpackCase
export CASE_NETCDF_VARIABLE_GROUPS, CASE_NETCDF_VARIABLES
export StepTimingStats, physics
export prescribed_case, synthetic_case, run_case
export StepWorkspace, threaded_workspaces
export ColumnarStepWorkspace
export continuous_bottom_deplete!
export update_surface_albedo!
export get_state, print_state
export add_timing!, timing_rows, print_timing_summary

end
