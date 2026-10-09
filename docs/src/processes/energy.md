```@meta
CurrentModule = Chion
```

# Energy Balance

`go_energy_flux!` solves the temperature profile of the snow/firn column and
the optional ice substrate below it with one implicit conductive step and a
linearized surface energy balance. The surface temperature `Tsrf` is the
temperature of the snow–air interface; it carries no heat capacity of its own
(Robin boundary). The formulas below follow `src/processes/energy_flux.jl` and
`src/processes/surface_fluxes.jl`.

For daily forcing, the default [Diurnal Shortwave Cycle](diurnal_cycle.md)
subdivides a forcing step before this energy solve is applied.

## Surface Flux Parameterization

At the surface, the model combines absorbed shortwave radiation, longwave
radiation, sensible heat exchange, and precipitation / latent-heat terms into
the surface energy balance ``Q(T_s) = Q_{sw} + Q_{lw} + Q_{sh} + Q_{p} + Q_{lh}``.
The terms are linearized about the previous surface temperature,

```math
Q(T_s) \approx Q_{\mathrm{const}} - Q_{\mathrm{lin}} T_s,
```

and stored as `surface_flux_constant` and `surface_flux_linear`.

### Robin Surface Boundary

The interface temperature balances the surface energy balance against the
conductive flux into the first layer,

```math
Q(T_s) = G_s\,(T_s - T_1),
\qquad
G_s = \frac{2K_1}{\Delta z_1},
```

where ``T_1``, ``K_1`` and ``\Delta z_1`` are the temperature, conductivity and
thickness of the first layer (its centre lies ``\Delta z_1/2`` below the
surface). Solving for the interface temperature,

```math
T_s = \frac{Q_{\mathrm{const}} + G_s T_1}{Q_{\mathrm{lin}} + G_s},
```

and eliminating it leaves the first layer as a regular finite-volume cell,

```math
c_i m_1 \frac{\partial T_1}{\partial t}
= G_s\,(T_s - T_1) + G_{3/2}\,(T_2 - T_1),
```

so the surface forcing enters only the top row of the implicit tridiagonal
system. A thin top layer therefore responds quickly without an artificial
surface heat capacity.

`BESSIModel` uses `turbulent_flux_scheme=:semix` by default for sensible and
latent heat. Set `turbulent_flux_scheme=:bessi` to recover the BESSI bulk
formulation. This choice is independent of `seb_scheme`, which defaults to
`:semix` for consistent graybody longwave exchange, including incident absorption.
Set `seb_scheme=:bessi` for the legacy longwave treatment. Prescribed
`q_sh` or `q_lh` forcing always takes precedence over either parameterization.

### Shortwave Radiation

If `q_sw_net` is not prescribed, absorbed shortwave is

```math
Q_{sw} = (1-\alpha)\,Q_{\mathrm{sw,down}},
```

where ``\alpha`` is the already diagnosed surface albedo stored in
`albedo[idx]`.

### Longwave Radiation

If `q_lw_down` is not prescribed, the default `longwave_scheme=:cloud_proxy`
derives it once per daily forcing step as
``q_{\mathrm{lw,down}} = \epsilon_{\mathrm{eff}}\,\sigma T_{\mathrm{air}}^4`` with

```math
\epsilon_{\mathrm{eff}} = \epsilon_0 + \epsilon_T (T_{\mathrm{air}} - T_0) + \epsilon_n n,
\qquad
n = 1 - \frac{\tau}{\tau_{\mathrm{clear}}(z)},
\qquad
\tau = \frac{Q_{\mathrm{sw,down}}}{Q_{\mathrm{sw,TOA}}},
```

where the daily top-of-atmosphere shortwave follows from latitude and season
and ``\tau_{\mathrm{clear}} = 0.85 + 0.075\,z/\mathrm{km}``. The defaults
(``\epsilon_0 = 0.624``, ``\epsilon_T = 0.0032\,\mathrm{K^{-1}}``,
``\epsilon_n = 0.613``) were fitted to daily MAR longwave over Greenland;
``\epsilon_{\mathrm{eff}}`` is relative to the near-surface air temperature and
can exceed one under warm overcast skies. In the polar night, or for sub-daily
forcing whose shortwave is not a daily mean, ``n`` is the constant
`lw_night_cloud_fraction` (0.389). No forcing beyond the standard shortwave
is needed. `longwave_scheme=:graybody` uses ``\epsilon_{\mathrm{air}}`` instead:

```math
Q_{\mathrm{lw}} =
\epsilon_{\mathrm{snow}}\sigma\left(\epsilon_{\mathrm{air}} T_{\mathrm{air}}^{4}
- T_{s}^{4}\right).
```

The emitted longwave term is nonlinear in surface temperature, so the solver
linearizes it about the previous-step temperature ``T_s^n``:

```math
(T_s^{n+1})^4 \approx 4(T_s^n)^3T_s^{n+1}-3(T_s^n)^4.
```

This gives the implemented constant and linear parts

```math
Q_{\mathrm{lw,const}} =
\epsilon_{\mathrm{s}}\left(
q_{\mathrm{lw,down}}
+ 3\sigma(T_s^n)^4
\right),
\qquad
Q_{\mathrm{lw,lin}} =
4\sigma\epsilon_{\mathrm{s}}(T_s^n)^3,
```

where ``q_{\mathrm{lw,down}}`` is the prescribed, cloud-proxy or graybody
incident longwave and ``\epsilon_{\mathrm{s}}`` is the surface emissivity:
`ϵ_snow` over snow and `eps_ice` over bare ice. Only outgoing longwave is
linearized; the energy diagnostics record the nonlinear flux at the resolved
surface temperature. `seb_scheme=:bessi` uses unit incident absorption
instead.

### Sensible Heat

If `q_sh` is not prescribed, the default SEMIX scheme uses an aerodynamic
resistance corrected for atmospheric stability:

```math
Q_{\mathrm{sh}} = \frac{\rho_{\mathrm{air}}c_{p,\mathrm{air}}}{r_a}
(T_{\mathrm{air}}-T_s).
```

The resistance is based on wind speed, measurement height, snow or ice
roughness (`semix_z0m_snow`, or `semix_z0m_ice` over bare ice), and a
bulk-Richardson stability correction, ``1/(1 + b\,Ri_b)`` for
stable stratification with `semix_stable_coefficient` ``b`` (default 40). The
configured `semix_sensible_exchange_factor` (default 2.5) scales this exchange.

With `turbulent_flux_scheme=:bessi`, the code instead uses

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

If `q_lh` is prescribed, the code treats it as an additional turbulent latent
heat flux. It is added to these snowfall/rain heat-term diagnoses rather than
replacing them. If `q_lh` is not prescribed and relative humidity is
available, the default SEMIX parameterization diagnoses humidity exchange
using the same stability-corrected aerodynamic resistance as sensible heat:

```math
Q_{\mathrm{vap}} =
\frac{L\,\rho_{\mathrm{air}}}{r_a}(q_{\mathrm{air}}-q_s).
```

The latent heat ``L`` follows the surface phase, and
`semix_latent_exchange_factor` scales the exchange. With
`turbulent_flux_scheme=:bessi`, Chion instead uses the vapor-pressure
parameterization

```math
Q_{\mathrm{vap}} =
\frac{D_{\mathrm{lf}}}{p_{\mathrm{air}}}
\left(RH\,e_{\mathrm{sat,w}}(T_{\mathrm{air}}) - e_{\mathrm{sat,i}}(T_s)\right),
```

with

```math
D_{\mathrm{lf}} =
r\,\frac{D_{\mathrm{sh}}}{c_{p,\mathrm{air}}}\,0.622\,(L_v + L_m).
```

Here ``r`` is `latent_heat_flux_ratio`, defaulting to 1.0. The saturation
vapor pressure over ice is linearized about the previous surface temperature
for the implicit surface solve. `relative_humidity` may be supplied either as a
fraction (`0.0` to `1.0`) or as percent (`0.0` to `100.0`). If neither `q_lh`
nor relative humidity is
available, the turbulent latent heat term is zero.

When both snowfall and rainfall are positive, the snowfall branch takes precedence in `_diagnose_latent_heat_flux_coefficients`.

### Combined Linearized Surface Forcing

Collecting the terms used by the implementation gives

```math
Q_{\mathrm{const}} =
Q_{sw}
+ Q_{\mathrm{lw,const}}
+ Q_{\mathrm{sh,const}}
+ K_{\mathrm{snow/rain/latent}},
\qquad
Q_{\mathrm{lin}} =
Q_{\mathrm{lw,lin}}
+ Q_{\mathrm{sh,lin}}
+ H_{\mathrm{snow/latent}}.
```

## Diffusion

The vertical diffusion equation is

```math
c_i \rho_s \frac{\partial T}{\partial t}
= \frac{\partial}{\partial z}\!\left( K(\rho_s)\,\frac{\partial T}{\partial z} \right).
```

### Thermal Conductivity

The implementation uses the temperature-scaled, smooth snow/firn blend from
Calonne et al. (2019), Eq. (5); it transitions between the snow and firn
regressions around 450 kg m⁻³.

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
\left(
\frac{\Delta z_i}{2K_i}
+ \frac{\Delta z_{i+1}}{2K_{i+1}}
\right)^{-1}
=
\frac{2K_iK_{i+1}}
{K_{i+1}\Delta z_i + K_i\Delta z_{i+1}}.
```

and we have

```math
\beta_i = -\frac{\Delta t}{c_i m_i},
```

so the off-diagonal coefficients are proportional to ``\beta_i G_{i\pm1/2}``.

### Ice Substrate

With `ice_substrate_layers > 0` (default 5), the solve continues below the
snow/firn layers into a fixed-geometry layer of glacier ice. Substrate layer
``k`` has density ``\rho_i`` and thickness
``\Delta z_k = 2^{k-1}\,\Delta z_{\mathrm{top}}`` with
`ice_substrate_top_thickness_m` ``= \Delta z_{\mathrm{top}}`` (default 0.05 m,
1.55 m in total for five layers). Its temperatures are stored in
`state.ice_temperature` and are solved in the same tridiagonal system as the
snow, with the usual interface conductance between the deepest snow layer and
the top of the substrate. The base of the substrate is insulated.

When a column has no snow, the top substrate layer forms the surface and the
same Robin boundary applies, with ice emissivity and roughness. Bare ice
therefore cools when the surface energy balance is negative and has to be
re-warmed to the melting point before it melts. Melt and rain on bare ice run
off. With `ice_substrate_layers=0`, bare ice is instead held at the melting
point and only positive surface energy is used.

### Discrete System

For interior layers, the implemented stencil can still be understood in the
usual tridiagonal form

```math
a_i T_{i-1}^{n+1} + b_i T_i^{n+1} + c_i T_{i+1}^{n+1} = r_i,
```

with coefficients assembled from the layer masses, conductivities, and
thicknesses. The surface row additionally includes the eliminated Robin
boundary term, while the bottom row (the deepest snow layer or the base of the
ice substrate) only sees conductive exchange with the layer above.

The multi-layer system is solved with the Thomas algorithm. 

## Melt-Point Constraint

If the interface temperature from the first solve exceeds `T0`, the code:

1. marks the step as needing melt
2. fixes the interface at ``T_s = T_0`` and solves the system again, keeping
   the conduction from the surface into the first layer,
   ``G_s (T_0 - T_1)``
3. clamps the final profile to `T <= T0`
4. returns the remaining surface energy as melt energy,

```math
E_{\mathrm{melt}} =
\max\!\left(0,\;
\left[Q_{\mathrm{const}} - Q_{\mathrm{lin}} T_0 - G_s (T_0 - T_1)\right]\Delta t
\right).
```

The energy that warms the snow below a melting surface is thus conducted into
the column instead of being melted, and no conductive flux leaves an
above-melting-point surface.

## API

```@docs
go_energy_flux!
```

```@docs; canonical=false
_go_energy_flux_resolved!
```
