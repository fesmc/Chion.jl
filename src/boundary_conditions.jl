#=
Boundary conditions of the column heat equation (as Laddie's `BoundaryConditions` slot of the
model). The surface condition closes row 1 of the implicit heat solve, the bottom condition the
last row (the snow/firn or, with a substrate, the bottom of the ice substrate); both enter the
enthalpy budget through `boundary_enthalpy`.
=#

"""
Abstract supertype of the surface boundary condition of the heat solve.

Available subtypes: [`SurfaceEnergyBalance`](@ref), [`PrescribedSurfaceFlux`](@ref),
[`PrescribedSurfaceTemperature`](@ref).
"""
abstract type AbstractSurfaceBoundaryCondition end

"""$(TYPEDEF)

Robin condition from the linearised surface energy balance, Q(T_s) = G_s (T_s − T₁) with
G_s = 2K₁/Δz₁, switching to the Dirichlet condition T_s = T₀ when the free surface would exceed
the melting point; the excess energy melts snow or ice."""
struct SurfaceEnergyBalance <: AbstractSurfaceBoundaryCondition end

"""
$(TYPEDEF)

Neumann condition: a constant heat flux F into the surface replaces the energy balance,
F = G_s (T_s − T₁); above the melting point the surface is held at T₀ and the excess melts.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct PrescribedSurfaceFlux{NF <: AbstractFloat} <: AbstractSurfaceBoundaryCondition
    "heat flux into the surface (W m⁻², default `0`)"
    flux::NF = 0.0
end

"""
$(TYPEDEF)

Dirichlet condition: the surface temperature is fixed, T_s = T, and the heat flux into the column
is G_s (T − T₁). No melt is diagnosed; use temperatures at or below T₀.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct PrescribedSurfaceTemperature{NF <: AbstractFloat} <: AbstractSurfaceBoundaryCondition
    "surface temperature (K, default `263.15`)"
    temperature::NF = 263.15
end

"""
Abstract supertype of the bottom boundary condition of the heat solve.

Available subtypes: [`Insulated`](@ref), [`GeothermalFlux`](@ref),
[`PrescribedBottomTemperature`](@ref).
"""
abstract type AbstractBottomBoundaryCondition end

"""$(TYPEDEF)

No heat flux through the bottom of the column (the default)."""
struct Insulated <: AbstractBottomBoundaryCondition end

"""
$(TYPEDEF)

Neumann condition: a constant heat flux q enters the bottom layer from below,
cᵢ mₙ dTₙ/dt = … + q.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct GeothermalFlux{NF <: AbstractFloat} <: AbstractBottomBoundaryCondition
    "upward heat flux into the column (W m⁻², default `0.05`)"
    flux::NF = 0.05
end

"""
$(TYPEDEF)

Dirichlet condition half a layer below the bottom row: the bottom layer exchanges heat
Gₙ (T_b − Tₙ) with Gₙ = 2Kₙ/Δzₙ, e.g. with the near-surface ice temperature of an ice-sheet model.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct PrescribedBottomTemperature{NF <: AbstractFloat} <: AbstractBottomBoundaryCondition
    "temperature below the column (K, default `263.15`)"
    temperature::NF = 263.15
end

for T in (:PrescribedSurfaceFlux, :PrescribedSurfaceTemperature, :GeothermalFlux, :PrescribedBottomTemperature)
    @eval $T(value::Real) = $T(float(value))
end

"""
$(TYPEDEF)

The boundary conditions of the column heat solve.

    BoundaryConditions(; surface = SurfaceEnergyBalance(), bottom = Insulated())

```jldoctest
julia> BoundaryConditions(bottom = GeothermalFlux(0.06)).bottom
GeothermalFlux(flux = 0.06)
```

# Fields
$(TYPEDFIELDS)
"""
struct BoundaryConditions{S <: AbstractSurfaceBoundaryCondition, B <: AbstractBottomBoundaryCondition}
    "surface condition (default `SurfaceEnergyBalance()`)"
    surface::S
    "bottom condition (default `Insulated()`)"
    bottom::B
end

BoundaryConditions(; surface::AbstractSurfaceBoundaryCondition = SurfaceEnergyBalance(),
    bottom::AbstractBottomBoundaryCondition = Insulated()) = BoundaryConditions(surface, bottom)

validate(bc::PrescribedSurfaceTemperature) = _require(isfinite(bc.temperature) && bc.temperature > 0,
    "`temperature` must be a positive temperature in K.")
validate(bc::PrescribedBottomTemperature) = _require(isfinite(bc.temperature) && bc.temperature > 0,
    "`temperature` must be a positive temperature in K.")
validate(bc::Union{PrescribedSurfaceFlux, GeothermalFlux}) = _require(isfinite(bc.flux), "`flux` must be finite.")
validate(bc::BoundaryConditions) = (validate(bc.surface); validate(bc.bottom); nothing)

Adapt.adapt_structure(to, bc::BoundaryConditions) = bc

function Base.show(io::IO, bc::Union{PrescribedSurfaceFlux, PrescribedSurfaceTemperature, GeothermalFlux,
        PrescribedBottomTemperature})
    name = only(fieldnames(typeof(bc)))
    print(io, nameof(typeof(bc)), "(", name, " = ", getfield(bc, name), ")")
    return nothing
end

# ---------------------------------------------------------------------------
# Kernel-side closures (see `_go_energy_flux_resolved!`)
# ---------------------------------------------------------------------------

"""The linearised surface flux `(constant, linear)`, Q ≈ constant − linear T_s, under `bc`."""
@inline _surface_flux_terms(::AbstractSurfaceBoundaryCondition, flux_constant, flux_linear) = flux_constant, flux_linear
@inline _surface_flux_terms(bc::PrescribedSurfaceFlux, flux_constant, flux_linear) = bc.flux, zero(flux_linear)

"""Conductance between the bottom row and the bottom boundary (zero without a Dirichlet bottom)."""
@inline _bottom_conductance(::AbstractBottomBoundaryCondition, conductivity, thickness) = zero(conductivity)
@inline _bottom_conductance(::PrescribedBottomTemperature, conductivity, thickness) =
    oftype(conductivity, 2) * conductivity / _safe_positive(thickness)

"""Off-diagonal coupling `β G_b` of the bottom row to the boundary (`zero` without a Dirichlet bottom)."""
@inline _bottom_coupling(::AbstractBottomBoundaryCondition, β, conductance) = zero(β)
@inline _bottom_coupling(::PrescribedBottomTemperature, β, conductance) = β * conductance

"""Temperature of the virtual row below the column (the `p` that starts the elimination)."""
@inline _bottom_temperature(::AbstractBottomBoundaryCondition, reference) = zero(reference)
@inline _bottom_temperature(bc::PrescribedBottomTemperature, reference) = oftype(reference, bc.temperature)

"""Right-hand side of the bottom row: the old temperature plus the Neumann source −β q = Δt q / (cᵢ mₙ)."""
@inline _bottom_rhs(::AbstractBottomBoundaryCondition, old_temperature, β) = old_temperature
@inline _bottom_rhs(bc::GeothermalFlux, old_temperature, β) = old_temperature - β * bc.flux

"""Energy (J m⁻²) that entered through the bottom over the step, added to the surface `heating`."""
@inline _with_bottom_heating(::Insulated, heating, conductance, bottom_temperature, dt_seconds) = heating
@inline _with_bottom_heating(bc::GeothermalFlux, heating, conductance, bottom_temperature, dt_seconds) =
    heating + bc.flux * dt_seconds
@inline _with_bottom_heating(bc::PrescribedBottomTemperature, heating, conductance, bottom_temperature, dt_seconds) =
    heating + conductance * (bc.temperature - bottom_temperature) * dt_seconds
