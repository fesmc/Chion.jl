"""
Scratch storage for stepping and energy-flux solves.
"""

"""
    _workspace_array(storage, ::Type{NF}, dims...)

Allocate scratch storage compatible with `storage`, preserving the active
backend while changing element type and shape.
"""
@inline _workspace_array(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N} =
    similar(storage, NF, dims...)

"""
    EnergyWorkspace

Scratch arrays reused by the implicit temperature solver in
[`go_energy_flux!`](@ref).
"""
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

"""
    EnergyWorkspace(storage, ::Type{NF}, dims...)

Allocate energy-solver scratch arrays on the same backend as `storage`.
Returns an `EnergyWorkspace` whose fields are mutated by the energy-flux
solver.
"""
function EnergyWorkspace(storage, ::Type{NF}, dims::Vararg{Int,N}) where {NF <: AbstractFloat, N}
    allocate() = _workspace_array(storage, NF, dims...)
    return EnergyWorkspace(allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate(), allocate())
end

"""
    EnergyWorkspace(domain)

Allocate energy-flux scratch storage sized for `domain`.
"""
EnergyWorkspace(domain::AbstractSnowpackDomain) =
    EnergyWorkspace(domain.mass, number_type(domain.c), domain.Ntot)

"""
    ColumnarStepWorkspace

Column-major scratch storage that holds temporary state for every column in a
batch run, regardless of whether the backing arrays live on CPU or GPU.
"""
struct ColumnarStepWorkspace{LWT,ET}
    liquid_water_before_energy::LWT
    energy::ET
end

"""
    ColumnarStepWorkspace(storage, ::Type{NF}, Ntot, ncol)

Allocate column-major scratch arrays that hold per-layer temporary state for
every column in a batch on the same backend as `storage`.
"""
function ColumnarStepWorkspace(storage, ::Type{NF}, Ntot::Int, ncol::Int) where {NF <: AbstractFloat}
    return ColumnarStepWorkspace(
        _workspace_array(storage, NF, Ntot, ncol),
        EnergyWorkspace(storage, NF, Ntot, ncol),
    )
end

"""
    ColumnarStepWorkspace(domain)

Allocate batch stepping scratch compatible with `domain`'s backend and sized
for all columns.
"""
function ColumnarStepWorkspace(domain::AbstractSnowpackDomain)
    NF = number_type(domain.c)
    return ColumnarStepWorkspace(domain.mass, NF, domain.Ntot, column_count(domain))
end

Adapt.@adapt_structure EnergyWorkspace
Adapt.@adapt_structure ColumnarStepWorkspace
