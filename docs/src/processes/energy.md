```@meta
CurrentModule = Chion
```

# Energy Balance

`go_energy_flux!` solves the snow temperature profile with an implicit 1D
conductive step plus linearized surface forcing. This page keeps the more
extensive derivation-style description, but the formulas below have been
checked against the current Julia implementation in `src/processes/energy_flux.jl`
and `src/processes/surface_fluxes.jl`.

## Surface Flux Parameterization

At the surface, the model combines absorbed shortwave radiation, longwave
radiation, sensible heat exchange, and precipitation / latent-heat terms into a
single surface forcing. In continuous form this can be written as

```math
c_i\, m_s \,\frac{\partial T_s}{\partial t}
= Q_{sw} + Q_{lw} + Q_{sh} + Q_{p} + Q_{lh},
```

with ice heat capacity ``c_i`` and surface-layer mass per area ``m_s``.

In the implementation, these terms are rewritten into a linearized form

```math
Q(T_s^{n+1}) \approx Q_{\mathrm{const}} - Q_{\mathrm{lin}} T_s^{n+1},
```

so that the single-layer update becomes

```math
T_s^{n+1}
=
\frac{T_s^n + \lambda Q_{\mathrm{const}}}
{1 + \lambda Q_{\mathrm{lin}}},
\qquad
\lambda = \frac{\Delta t}{c_i m_s}.
```

For multi-layer columns, the same surface forcing enters the top row of the
implicit tridiagonal diffusion system.

### Shortwave Radiation

If `q_sw_net` is not prescribed, absorbed shortwave is

```math
Q_{sw} = (1-\alpha)\,Q_{\mathrm{sw,down}},
```

where ``\alpha`` is the already diagnosed surface albedo stored in
`albedo_dynamic[idx]`.

### Longwave Radiation

If `q_lw_down` is not prescribed, the code uses

```math
Q_{\mathrm{lw}} =
\sigma\left(\epsilon_{\mathrm{air}} T_{\mathrm{air}}^{4}
- \epsilon_{\mathrm{snow}} T_{s}^{4}\right).
```

The emitted longwave term is nonlinear in surface temperature, so the solver
linearizes it about the previous-step temperature ``T_s^n``:

```math
(T_s^{n+1})^4 \approx 4(T_s^n)^3T_s^{n+1}-3(T_s^n)^4.
```

This gives the implemented constant and linear parts

```math
Q_{\mathrm{lw,const}} =
\sigma\left(
\epsilon_{\mathrm{air}}T_{\mathrm{air}}^4
+ 3\epsilon_{\mathrm{snow}}(T_s^n)^4
\right),
```

```math
Q_{\mathrm{lw,lin}} =
4\sigma\epsilon_{\mathrm{snow}}(T_s^n)^3.
```

If `q_lw_down` is prescribed, the incoming longwave term is replaced by that
value and only the outgoing ``\epsilon_{\mathrm{snow}} T_s^4`` term is
linearized.

### Sensible Heat

If `q_sh` is not prescribed, the code uses

```math
Q_{\mathrm{sh}}= D_{\mathrm{sh}}(T_{\mathrm{air}}-T_{s}),
```

with default ``D_{\mathrm{sh}} = 10\,\mathrm{W\,m^{-2}\,K^{-1}}``. In the
linearized form,

```math
Q_{\mathrm{sh,const}} = D_{\mathrm{sh}}T_{\mathrm{air}},
\qquad
Q_{\mathrm{sh,lin}} = D_{\mathrm{sh}}.
```

### Precipitation And Latent-Heat Terms

Snowfall and rain heat terms are directly diagnosed from the forcing
rates in `kg m^-2 s^-1`.

For snowfall, 

```math
H_{\mathrm{snow}} = P_{\mathrm{snow}} c_i,
\qquad
K_{\mathrm{snow}} = P_{\mathrm{snow}} c_i T_{\mathrm{air}},
```

so the snowfall contribution can be written as

```math
Q_{\mathrm{snow}}(T_s) =
K_{\mathrm{snow}} - H_{\mathrm{snow}}T_s
=
P_{\mathrm{snow}} c_i (T_{\mathrm{air}} - T_s).
```

For rainfall onto an existing snow surface, the current implementation uses a
temperature-independent term

```math
Q_{\mathrm{rain}} =
P_{\mathrm{rain}} c_w (T_{\mathrm{air}} - T_0).
```

If `q_lh` is prescribed, the code uses it directly and disables these internal
snowfall/rain heat-term diagnoses.

One implementation detail worth stating explicitly: when both snowfall and
rainfall are positive, the snowfall branch takes precedence in
`_diagnose_latent_heat_flux_coefficients`.

### Combined Linearized Surface Forcing

Collecting the terms used by the implementation gives

```math
Q_{\mathrm{const}} =
Q_{sw}
+ Q_{\mathrm{lw,const}}
+ Q_{\mathrm{sh,const}}
+ K_{\mathrm{snow/rain/latent}},
```

```math
Q_{\mathrm{lin}} =
Q_{\mathrm{lw,lin}}
+ Q_{\mathrm{sh,lin}}
+ H_{\mathrm{snow/latent}}.
```

In code, these are stored as `surface_flux_constant` and
`surface_flux_linear`.

## Diffusion

The vertical diffusion equation is

```math
c_i \rho_s \frac{\partial T}{\partial t}
= \frac{\partial}{\partial z}\!\left( K(\rho_s)\,\frac{\partial T}{\partial z} \right).
```

### Thermal Conductivity Options

`diffusion_model` selects one of three conductivity closures:

- `1`: Yen-style
- `2`: Sturm piecewise
- `3` or anything else: Van Dusen-style polynomial

The implemented formulas are:

```math
K(\rho) = K_i \left(\frac{\rho}{1000}\right)^{1.88}
```

for `diffusion_model = 1`,

```math
K(\rho)=
\begin{cases}
0.023 + 0.234\times 10^{-3}\rho, & \rho \le 156 \\
0.138 - 1.01\times 10^{-3}\rho + 3.233\times 10^{-6}\rho^2, & \rho > 156
\end{cases}
```

for `diffusion_model = 2`, and

```math
K(\rho)=2.1\times 10^{-2} + 4.2\times 10^{-4}\rho + 2.2\times 10^{-9}\rho^3
```

for the fallback branch.

### Interface Conductance

For each layer, the code first computes the layer thickness

```math
\Delta z_i = \frac{m_i}{\rho_i}.
```

For the interface between neighboring layers ``i`` and ``i+1``, the helper
`interface_conductance` computes

```math
G_{i+1/2}
=
\frac{K_i \Delta z_i + K_{i+1}\Delta z_{i+1}}
{(\Delta z_i + \Delta z_{i+1})^2}.
```

and we have

```math
\beta_i = -\frac{2\Delta t}{c_i m_i},
```

so the off-diagonal coefficients are proportional to ``\beta_i G_{i\pm1/2}``.

### Discrete System

For interior layers, the implemented stencil can still be understood in the
usual tridiagonal form

```math
a_i T_{i-1}^{n+1} + b_i T_i^{n+1} + c_i T_{i+1}^{n+1} = r_i,
```

with coefficients assembled from the layer masses, conductivities, and
thicknesses. The surface row additionally includes the linearized surface
forcing, while the bottom row only sees conductive exchange with the layer
above.

The multi-layer system is solved with the Thomas algorithm. 

## Melt-Point Constraint

If the solved surface temperature exceeds `T0`, the code:

1. marks the step as needing melt
2. diagnoses the energy required to bring the surface to `T0`
3. reruns the solve with the surface fixed at `T0`
4. clamps the final profile to `T <= T0`
5. returns any remaining positive surface energy as `melt_energy_available`

This two-pass handling avoids unphysical conductive fluxes from an
above-melting-point surface into deeper cold snow.

## API

```@docs
go_energy_flux!
```

```@docs; canonical=false
_go_energy_flux_resolved!
_debm_melt_window_fluxes
```
