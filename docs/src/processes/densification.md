```@meta
CurrentModule = Chion
```

# Densification

`go_densification!` updates layer density using the scheme encoded in
`SnowpackPhysicalConstants.low_density_densification`.

## Common Update Rule

For each active layer, the code computes a density tendency ``\dot{\rho}`` and
then applies

```math
\rho^{n+1} =
\min\left(
\rho_i,\;
\max\left(\rho^n,\; \rho^n + \dot{\rho}\Delta t\right)
\right).
```

The update is nondecreasing and capped at ice density ``\rho_i``.

The overburden pressure used by the densification laws is

```math
\sigma = g\left(M_{\mathrm{above}} + \frac{m}{2}\right),
```

with `g = 9.81 m s^-2`.

## Low-Density Branch (`rho < 550 kg m^-3`)

### BESSI Option (`:bessi`)

The BESSI-style tendency is

```math
\dot{\rho} =
0.011
\exp\!\left(-\frac{10160}{8.13\,T}\right)
(\rho_i - \rho)\max(A_t, 0),
```

where `A_t` is the accumulation proxy passed into `go_densification!`.

### HTESSEL Option (`:htessel`)

The HTESSEL-style tendency is

```math
\dot{\rho} = \rho \left(\frac{\sigma}{\eta} + \xi \right),
```

with

```math
\eta = 3.7\times 10^7
\exp\!\left(8.1\times 10^{-2}(T_0 - T) + 1.8\times 10^{-2}\rho\right),
```

and

```math
\xi = 2.8\times 10^{-6}
\exp\!\left(-4.2\times 10^{-2}(T_0 - T) - 460\max(0,\rho - 150)\right).
```

When HTESSEL is active, retained liquid water can also compact low-density snow
after the energy, melt, and percolation updates:

```math
\rho^{new} =
\min\left(
\rho_i,\;
\max\left(\rho,\; \rho + \rho\frac{\Delta m_w}{m}\right)
\right),
\qquad
\Delta m_w = \max\left(m_w^{after} - m_w^{before}, 0\right).
```

## Intermediate-Density Branch (`550 <= rho < 800 kg m^-3`)

The code uses

```math
\dot{\rho} =
25400
\exp\!\left(-\frac{60000}{8.13\,T}\right)
\rho\,f(\rho)\,(\Delta P)^3,
```

with

```math
f(\rho)=10^{-29.166(\rho/\rho_i)^3 + 84.422(\rho/\rho_i)^2 - 87.425(\rho/\rho_i) + 30.673}.
```

The effective pressure is

```math
\Delta P = P_{\mathrm{ice}} - P_{\mathrm{bubble}},
```

where

```math
P_{\mathrm{ice}} = \frac{9.81\,(M_{\mathrm{above}} + m/2)}{10^6}
```

and

```math
P_{\mathrm{bubble}}=
\frac{P_{\mathrm{atm}}}{10^6}
\left(
\frac{(1/\rho_e-1/\rho_i)}{(1/\rho-1/\rho_i)}-1
\right).
```

The current implementation uses `rho_e = 815 kg m^-3` and `P_atm = 101325 Pa`
internally.

## High-Density Branch (`rho >= 800 kg m^-3`)

The same Arrhenius prefactor is used with the alternative shape factor

```math
f(\rho)=\frac{3}{16}\,
\frac{1-\rho/\rho_i}{
\left(1-(1-\rho/\rho_i)^{1/3}\right)^3}.
```

## API

```@docs
go_densification!
```

```@docs; canonical=false
_htessel_low_density_rate
_apply_htessel_liquid_water_compaction!
```
