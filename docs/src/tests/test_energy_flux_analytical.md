# Energy Flux Analytical Validation Note

This page records the analytical checks used to verify the current
energy-balance documentation. These notes are retained as validation material;
they are not currently backed by a file in `test/`.

## Single-Layer Closed-Form Update

For one snow layer, the implemented implicit-Euler update is

```math
T^{n+1} = \frac{T^n + \lambda F_{\mathrm{const}}}{1 + \lambda F_{\mathrm{lin}}},
\quad
\lambda = \frac{\Delta t}{c_i m_s}.
```

This matches the single-layer branch in `src/processes/energy_flux.jl`.

## Uniform Profile Under Pure Diffusion

When all external surface-flux terms are disabled and the temperature profile is
uniform, the diffusion solve should preserve that uniform state up to
floating-point roundoff.

## Sensible-Energy Conservation

Under the same pure-diffusion setup, diffusion should redistribute heat
internally while conserving total sensible energy:

```math
E = \sum_i m_i c_i T_i.
```
