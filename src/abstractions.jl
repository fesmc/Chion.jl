"""
Shared abstract domain and variable metadata types.
"""

abstract type AbstractSnowpackDomain{NF} end

struct PrognosticVariable{name}
    description::String
end

struct AuxiliaryVariable{name}
    description::String
end

variables(::Any) = ()
compute_auxiliary!(args...) = nothing
