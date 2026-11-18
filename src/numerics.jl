#=
Numerical guards used inside kernels.

Every guard is branch-free (`ifelse`, no `? :`) so that Reactant can trace it, and keeps the value
of the unguarded expression wherever that is defined, so that a guard never changes results. The
guards remove non-differentiable or undefined points instead of introducing new kinks.
=#

"""
    _safe_positive(x)

`x`, or `EPS_TINY` when `x ≤ EPS_TINY` (also for `NaN`): a positive denominator.
"""
@inline _safe_positive(x) = ifelse(x > EPS_TINY, x, oftype(x, EPS_TINY))

"""
    _safe_nonnegative(x)

Clamp `x` to zero from below (`NaN` gives zero).
"""
@inline _safe_nonnegative(x) = ifelse(x > zero(x), x, zero(x))

"""
    _safe_pow(x, p)

`xᵖ` for `x ≥ 0` and `p > 0`, with `0ᵖ = 0` reached without evaluating the power at zero. Plain
`x^p` is exact there too, but its derivative with respect to `p` is `xᵖ log x`, which is
`0 · (−∞) = NaN` at `x = 0` and breaks forward- and reverse-mode AD. The inner `ifelse` keeps the
power away from zero in both branches ("double where"); `max(x, floatmin)^p` would instead overflow
the derivative with respect to `x` for `p < 1`.
"""
@inline function _safe_pow(x, p)
    positive = x > zero(x)
    return ifelse(positive, ifelse(positive, x, one(x))^p, zero(x))
end
