module Chion

using Printf

include("SnowpackModel.jl")
using .SnowpackModel

export SnowpackColumn, step!
export SnowpackPhysicalConstants
export continuous_bottom_deplete!
export get_state, print_state

end
