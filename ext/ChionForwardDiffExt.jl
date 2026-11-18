"""Automatic differentiation of models built at a ForwardDiff `Dual` number type, loaded by `using ForwardDiff`."""
module ChionForwardDiffExt

using Chion
using ForwardDiff

# The value of a dual number (nested duals recursively), for host-side conversions.
Chion._primal(x::ForwardDiff.Dual) = Chion._primal(ForwardDiff.value(x))

end
