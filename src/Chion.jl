module Chion


import Printf; using Printf



export Prinf




# SnowpackModel playground
include("SnowpackModel.jl")
using .SnowpackModel    # Needed so we can export names from sub-modules at the top-level

export SnowpackColumn, step!
export SnowpackPhysicalConstants
export get_state, print_state

export calc_density_gradient_HL80
export calc_density_gradient_powerlaw_ref

end