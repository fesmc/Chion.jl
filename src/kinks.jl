#=
Treatment of the kinks of the column physics: the switches of `max`, `min` and `clamp` that the
state crosses. As Laddie's `AbstractKinks`: exact by default (bit-identical), or smoothed over a
band of half-width ε per physical quantity for smoother gradients.
=#

"""
Abstract supertype of how the model treats its kinks, the switches of `max`, `min` and `clamp`.
Pass an instance as `Params(kinks = …)`: [`ExactKinks`](@ref) (the default) or [`SmoothKinks`](@ref).

The kinks make the model piecewise smooth: a derivative exists almost everywhere but jumps where a
switch flips, and gradient-based calibration can stall there. Discrete events (layer split and
merge, snow ↔ bare ice) are not kinks in this sense and stay as they are.
"""
abstract type AbstractKinks end

"""$(TYPEDEF)

Every `max`, `min` and `clamp` as it is (the default): bit-identical to a model without a kink
treatment."""
struct ExactKinks <: AbstractKinks end

"""
$(TYPEDEF)

Smooth the kinks over a band of half-width ε, one ε per physical quantity: `max(x, 0)` becomes the
C¹ quadratic blend (x + ε)² / 4ε for |x| < ε and is exact outside the band; `max(a, b)`,
`min(a, b)` and `clamp` follow from it. An ε of zero keeps that group exact.

| ε | where |
|---|---|
| `flux` (W m⁻²) | melt energy of bare ice, max(Q_net, 0) |
| `wind` (m s⁻¹) | floor of the wind speed (0.1 m s⁻¹) of the SEMIX turbulent fluxes |
| `fraction` (–) | clamp of the shortwave deficit of the cloud-proxy longwave to [0, 1] |

Not smoothed are bounds that the state sits on rather than crosses, where a blend would bias every
step by about ε/4 (albedo at `alpha_wet`/`alpha_dry`, temperatures clamped to T₀), and the power of
the cloud fraction at zero.

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct SmoothKinks{NF <: Real} <: AbstractKinks
    "half-width for heat fluxes (W m⁻², default `1`)"
    flux::NF = 1.0
    "half-width for wind speeds (m s⁻¹, default `0.05`)"
    wind::NF = 0.05
    "half-width for fractions (–, default `0.01`)"
    fraction::NF = 0.01
end

SmoothKinks(values::Real...) = SmoothKinks(promote(map(float, values)...)...)

function validate(k::SmoothKinks)
    _require(all(name -> getfield(k, name) >= 0, fieldnames(SmoothKinks)), "SmoothKinks widths must be non-negative.")
    return nothing
end

Base.show(io::IO, k::SmoothKinks) =
    print(io, "SmoothKinks(", join(("$name = $(getfield(k, name))" for name in fieldnames(SmoothKinks)), ", "), ")")

"""`max(x, 0)` blended quadratically over |x| < ε (exact for ε = 0)."""
@inline function _smooth_max0(x, ε)
    blend = (x + ε)^2 / _safe_positive(4 * ε)
    return ifelse(abs(x) < ε, blend, max(x, zero(x)))
end

"""
    _kmax(kinks, a, b, group)

`max(a, b)` under the kink treatment `kinks`; `group` (`Val(:flux)`, …) names the ε of `SmoothKinks`.
"""
@inline _kmax(::ExactKinks, a, b, group) = max(a, b)
@inline _kmax(k::SmoothKinks, a, b, ::Val{group}) where {group} = b + _smooth_max0(a - b, oftype(a, getfield(k, group)))

"""`min(a, b)` under the kink treatment (see `_kmax`)."""
@inline _kmin(::ExactKinks, a, b, group) = min(a, b)
@inline _kmin(k::SmoothKinks, a, b, ::Val{group}) where {group} = a - _smooth_max0(a - b, oftype(a, getfield(k, group)))

"""`clamp(x, lo, hi)` under the kink treatment (see `_kmax`)."""
@inline _kclamp(::ExactKinks, x, lo, hi, group) = clamp(x, lo, hi)
@inline _kclamp(k::SmoothKinks, x, lo, hi, group) = _kmin(k, _kmax(k, x, lo, group), hi, group)
