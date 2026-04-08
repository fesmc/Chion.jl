"""
Scratch storage for stepping and energy-flux solves.
"""

@inline _workspace_array(::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} = zeros(NF, dims...)
@inline _workspace_array(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} =
    similar(storage, NF, dims...)

struct EnergyWorkspace{LT,DT,UT,RT,IT,PT,TT,KT}
    lower::LT
    diag::DT
    upper::UT
    rhs::RT
    interface_conductance::IT
    previous_temperature::PT
    layer_thickness::TT
    thermal_conductivity::KT
end

function EnergyWorkspace(::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

function EnergyWorkspace(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(storage, NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

EnergyWorkspace(domain::AbstractSnowpackDomain) =
    EnergyWorkspace(domain.mass, number_type(domain.c), domain.Ntot)

struct StepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

function StepWorkspace(::Type{NF}, Ntot::Int) where {NF <: AbstractFloat}
    return StepWorkspace(
        _workspace_array(NF, Ntot),
        EnergyWorkspace(NF, Ntot),
    )
end

function StepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    return StepWorkspace(
        _workspace_array(domain.mass, NF, domain.Ntot),
        EnergyWorkspace(domain.mass, NF, domain.Ntot),
    )
end

threaded_workspaces(domain::AbstractSnowpackDomain) = [StepWorkspace(domain) for _ in 1:Threads.maxthreadid()]

struct ColumnarStepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

function ColumnarStepWorkspace(::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return ColumnarStepWorkspace(
        _workspace_array(NF, Ntot, ncol),
        EnergyWorkspace(NF, Ntot, ncol),
    )
end

function ColumnarStepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    return ColumnarStepWorkspace(
        _workspace_array(domain.mass, NF, domain.Ntot, column_count(domain)),
        EnergyWorkspace(domain.mass, NF, domain.Ntot, column_count(domain)),
    )
end

Adapt.@adapt_structure EnergyWorkspace
Adapt.@adapt_structure StepWorkspace
Adapt.@adapt_structure ColumnarStepWorkspace
