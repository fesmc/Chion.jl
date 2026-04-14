module Chion

using Printf

include("SnowpackModel.jl")
using .SnowpackModel
include("case_runtime.jl")
include("case_api.jl")

export step!
export SnowpackDomain
export SnowpackStepForcing
export SnowpackStepFields
export SnowpackPhysicalConstants
export SnowpackStateFields
export RunConfig, RunResult
export SnowpackCase
export CASE_NETCDF_VARIABLE_GROUPS, CASE_NETCDF_VARIABLES
export TimingStats, physics
export prescribed_case, synthetic_case, run_case
export StepWorkspace, threaded_workspaces
export ColumnarStepWorkspace
export continuous_bottom_deplete!
export update_surface_albedo!
export get_state, print_state
export StepTimingStats, add_timing!, timing_rows, print_timing_summary

end
