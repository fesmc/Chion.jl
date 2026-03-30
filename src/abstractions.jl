"""
Terrarium-style abstract interfaces for snowpack processes, models, and states.
"""

abstract type AbstractProcess{NF} end
abstract type AbstractModel{NF} end
abstract type AbstractSnowModel{NF} <: AbstractModel{NF} end

abstract type AbstractSnowpackState{NF} end
abstract type AbstractSnowpackDomain{NF} <: AbstractSnowpackState{NF} end

struct PrognosticVariable{name}
    description::String
end

struct AuxiliaryVariable{name}
    description::String
end

@inline varname(::PrognosticVariable{name}) where {name} = name
@inline varname(::AuxiliaryVariable{name}) where {name} = name

variables(::Any) = ()
compute_auxiliary!(args...) = nothing
compute_tendencies!(args...) = nothing
