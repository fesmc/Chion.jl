```@meta
CurrentModule = Chion.SnowpackModel
```

# Densification

`go_densification!` updates layer density ``\rho`` using a multi-regime parameterization adapted from BESSI using an explicit Euler step.

## Model Formulation

For each active layer ``m``, density ``\rho`` at time ``n+1`` is updated by:

```math
\rho_m^{n+1} = \min\left(\rho_i,\; \max\left(\rho_m^n,\; \rho_m^n + \dot{\rho}_m \Delta t\right)\right)
```

where ``\rho_i`` is the ice density and ``\Delta t`` is the timestep. The density is clamped to stay non-decreasing and below the density of ice ``\rho_i``.

The densification process is divided in three distinct regimes. 

### Regime 1: Low density (``\rho < 550``)

```math
\dot{\rho} =
0.011 \mathrm{m^2kg^{-1}} \exp\!\left(-\frac{10160\mathrm{Jmol^{-1}}}{8.13\mathrm{Jmol^{-1}K^{-1}}\,T}\right)
(\rho_i-\rho)\,\max(A_t,0)
```

### Regime 2: HL option (``\rho \ge 550``, `hl=true`)
Not in use right now, not in BESSI description paper either. 
```math
\dot{\rho} =
0.575 \exp\!\left(-\frac{21400}{8.13\,T}\right)
(\rho_i-\rho)\,C_{yr}^{1/2}\,\max(A_t,0)^{1/2}
```

with ``C_{yr}=1000/(3600\cdot24\cdot365)``.

### Regime 3: Intermediate density (``550 \le \rho < 800``, `hl=false`)
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

### Regime 4: High density (``\rho \ge 800``, `hl=false`)
At pressures above ``\rho \ge 800``, the same empirical law is used but with a different form of the function ``f``:

```math
f(\rho)=\frac{3}{16}\,
\frac{1-\rho/\rho_i}{
\left(1-(1-\rho/\rho_i)^{1/3}\right)^3}
```

## Inputs

- `At`: accumulation proxy used by the low-density and HL regimes.
- `hl`: switches to Herron-Langway behavior above ``550\ \mathrm{kg\,m^{-3}}``.
- `rho_e`, `P_atm`: parameters used in bubble-pressure correction.

## API

```@docs
go_densification!(
    column::SnowpackColumn,
    At::Float64,
    dt_sec::Float64;
    hl::Bool=false,
    rho_e::Float64=815.0,
    P_atm::Float64=101325.0,
)
```
[^barnola]: https://b.tellusjournals.se/articles/10.3402/tellusb.v43i2.15249 