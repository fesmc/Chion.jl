# Energy Flux Analytical Tests

This page documents the three analytical checks in `test/test_energy_flux_analytical.jl`.

## Test 1: Single-layer closed-form update

For one snow layer, the solver should match the exact implicit-Euler solution:

```math
T^{n+1} = \frac{T^n + \lambda F_{\mathrm{const}}}{1 + \lambda F_{\mathrm{lin}}},
\quad
\lambda = \frac{\Delta t}{c_i m_s}.
```

The test compares model output to this formula with strict tolerance (`1e-12`).

## Test 2: Uniform profile is steady under pure diffusion

With all external fluxes disabled (``q_{sw}=q_{lw}=q_{sh}=q_{lh}=0`` and ``D_sh=\epsilon_{air}=\epsilon_{snow}=0``), a uniform temperature profile is an exact steady state.  
The test verifies that every active layer remains unchanged (up to floating-point roundoff).

## Test 3: Sensible-energy conservation under pure diffusion

Under the same zero-forcing setup, diffusion should only redistribute heat internally.  
The test checks conservation of total sensible energy:

```math
E = \sum_i m_i c_i T_i.
```

## Run command

```bash
julia --project=. test/test_energy_flux_analytical.jl
```

## Results snapshot

Recorded on **2026-03-12**:

| Test block | Passed assertions |
| --- | ---: |
| Single-layer closed-form update | 5/5 |
| Uniform profile steady state | 2/2 |
| Pure-diffusion energy conservation | 2/2 |
| **Total** | **9/9** |

Terminal summary:

```text
Activating project at `~/Documents/Chion.jl`
Test Summary:                | Pass  Total  Time
Energy Flux Analytical Cases |    9      9  0.6s
```
