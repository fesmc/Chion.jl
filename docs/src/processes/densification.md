```@meta
CurrentModule = Chion.SnowpackModel
```

# Densification

`go_densification!` updates layer density ``\rho`` using a multi-regime parameterization. For low-density snow (``\rho < 550\ \mathrm{kg\,m^{-3}}``) the model can use either the BESSI law or the HTESSEL parameterisation, selected through `SnowpackPhysicalConstants(low_density_densification=...)`. Higher-density firn follows the BESSI-style parameterisation.

## Model Formulation

For each active layer ``m``, density ``\rho`` at time ``n+1`` is updated by:

```math
\rho_m^{n+1} = \min\left(\rho_i,\; \max\left(\rho_m^n,\; \rho_m^n + \dot{\rho}_m \Delta t\right)\right)
```

where ``\rho_i`` is the ice density and ``\Delta t`` is the timestep. The density is clamped to stay non-decreasing and below the density of ice ``\rho_i``.

The densification process is divided into low-density and higher-density branches.

### Regime 1: Low density (``\rho < 550``)

Two alternative parameterisations are available. We decided to implement an alternative to the BESSI choice low-density densification parameterisation. This is mostly due to the fact that BESSI's low-density law is underestimating compactification for seasonal snowpacks. In BESSI, the fresh-snow density was fixed at ``350\,\mathrm{kg\,m^{-3}}``, which is unrealisticly high for most fresh-snow regimes. 

#### Option A: BESSI (default, `low_density_densification = :bessi`)

```math
\dot{\rho} =
0.011 \mathrm{m^2kg^{-1}} \exp\!\left(-\frac{10160\mathrm{Jmol^{-1}}}{8.13\mathrm{Jmol^{-1}K^{-1}}\,T}\right)
(\rho_i-\rho)\,\max(A_t,0)
```

This is the original low-density branch. It depends on the accumulation proxy ``A_t`` and produces no compaction when ``A_t \le 0``.

#### Option B: HTESSEL (`low_density_densification = :htessel`)

For the HTESSEL option, the densification rate is given by

```math
\dot{\rho} = \rho \left(\frac{\sigma}{\eta} + \xi \right),
```

where ``\sigma`` is the overburden stress at the layer midpoint,

```math
\sigma = g\left(M_{\mathrm{overburden}} + \frac{m}{2}\right),
```

with ``g = 9.81\ \mathrm{m\,s^{-2}}``, overlying mass ``M_{\mathrm{overburden}}`` and layer mass ``m``. The effective snow viscosity is

```math
\eta = 3.7\times 10^7
\exp\!\left(8.1\times 10^{-2}(T_0-T) + 1.8\times 10^{-2}\rho\right),
```

and the thermal-metamorphism term is

```math
\xi = 2.8\times 10^{-6}
\exp\!\left(-4.2\times 10^{-2}(T_0-T) - 460\max(0,\rho-150)\right).
```

Here ``T_0`` is the freezing point of water. In this branch the low-density densification rate depends on temperature, density and overburden stress, not only on the accumulation proxy ``A_t``.

When HTESSEL is active, retained meltwater also produces an additional compaction step for low-density snow after the energy and melt update. If liquid water in a layer increases from ``m_w^{\mathrm{before}}`` to ``m_w^{\mathrm{after}}`` while the solid mass ``m`` stays fixed, the bulk density is updated as

```math
\rho^{\mathrm{new}} =
\min\left(\rho_i,\; \max\left(\rho,\; \rho + \rho\frac{\Delta m_w}{m}\right)\right),
\quad
\Delta m_w = \max\left(m_w^{\mathrm{after}} - m_w^{\mathrm{before}}, 0\right).
```

This extra compaction is only applied for ``\rho < 550\ \mathrm{kg\,m^{-3}}``.


### Regime 2: Intermediate density (``550 \le \rho < 800``, `hl=false`)
For densities between ``550 \le \rho < 800``, a semi-empirical model is used[^barnola]
```math
\frac{\dot{\rho}}{\rho} = A\,f\,(\Delta P)^n
= A_0 \exp\!\left(-\frac{Q}{RT}\right)\, f\,(\Delta P)^n \, ,
```
with the constant ``A_0=25400\mathrm{MPa^{-3}s^{-1}}``, the activation energy ``Q=60\mathrm{KJmol^{-1}}``, the gas constant ``R``, the temperature ``T``, the pressure difference between overburden pressure and inside gas pressure ``\Delta P`` and the empirical function ``f``, which is based on the spherical pore model of Wilkinson and Ashby (1975). The exponent is chosen to be ``n=3``. 
The function ``f`` is given by:

```math
f(\rho)=10^{-29.166(\rho/\rho_i)^3 + 84.422(\rho/\rho_i)^2 - 87.425(\rho/\rho_i) + 30.673}.
```

The pressure difference ``\Delta P`` (or effective pressure) is given by the difference between the vertical overburden pressure ``P_{\mathrm{ice}}`` and the inside gas pressure ``P_{\mathrm{bubble}}`` (only above a threshold density of ``\rho_e=815\mathrm{kgm^{-3}}``).
```math
\Delta P = P_{\mathrm{ice}} - P_{\mathrm{bubble}}
```

where:

```math
P_{\mathrm{ice}}=\frac{9.81\mathrm{ms^{-2}}\,(M_{\mathrm{overburden}}+m/2)}{10^6},
\quad
P_{\mathrm{bubble}}=
\frac{P_{\mathrm{atm}}}{10^6}\left(
\frac{(1/\rho_e-1/\rho_i)}{(1/\rho-1/\rho_i)}-1
\right)\;\;\;\mathrm{if}\;\;\;(\rho>\rho_e)
```

with the mass above ``M_{\mathrm{overburden}}``, the midpoint of the current layer ``m/2`` and the atmospheric pressure ``P_{\mathrm{atm}}``.

### Regime 3: High density (``\rho \ge 800``, `hl=false`)
At pressures above ``\rho \ge 800``, the same empirical law is used but with a different form of the function ``f``:

```math
f(\rho)=\frac{3}{16}\,
\frac{1-\rho/\rho_i}{
\left(1-(1-\rho/\rho_i)^{1/3}\right)^3}
```

## Inputs

- `accumulation_rate`: accumulation proxy used by the BESSI low-density branch and the optional HL branch. It is ignored by the HTESSEL low-density parameterisation.
- `hl`: switches to Herron-Langway behavior above ``550\ \mathrm{kg\,m^{-3}}``.
- `rho_e`, `P_atm`: parameters used in bubble-pressure correction.
- `column.c.low_density_densification`: selects `:bessi` or `:htessel` for the ``\rho < 550\ \mathrm{kg\,m^{-3}}`` branch.

## API

```@docs; canonical=false
go_densification!
```
[^barnola]: https://b.tellusjournals.se/articles/10.3402/tellusb.v43i2.15249 
