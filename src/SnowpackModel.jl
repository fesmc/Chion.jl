"""
Column-based snowpack model with dynamic layering.
Based on Born et al. (2019) algorithm.

This module keeps the public API in one place and delegates implementation to
focused source files by responsibility.
"""

module SnowpackModel

using Printf

include("model_constants.jl")

export SnowpackPhysicalConstants
export SnowpackColumn
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

include("snowpack_types.jl")
include("column_helpers.jl")
include("albedo.jl")
include("mass_balance.jl")
include("energy_flux.jl")
include("densification.jl")
include("percolation.jl")
include("refreezing.jl")
include("timing.jl")
include("step.jl")
include("state_io.jl")

end # module
