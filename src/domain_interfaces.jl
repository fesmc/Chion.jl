"""
Shared abstract domain and variable metadata types.
"""

"""
    AbstractSnowpackDomain{NF}

Abstract supertype for array-backed snowpack state containers that expose the
core layer arrays and auxiliary diagnostics used by the process kernels.
"""
abstract type AbstractSnowpackDomain{NF} end

struct PrognosticVariable{name}
    description::String
end

struct AuxiliaryVariable{name}
    description::String
end

"""
    variables(x)

Return a tuple describing the prognostic and auxiliary variables exposed by
`x`. The generic fallback returns an empty tuple for objects that do not
publish variable metadata.
"""
variables(::Any) = ()

"""
    compute_auxiliary!(args...)

Generic fallback for updating auxiliary diagnostics in-place. Concrete domain
types override this when they expose derived state such as snow cover or
surface albedo.
"""
compute_auxiliary!(args...) = nothing
