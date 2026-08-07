```@meta
CurrentModule = Chion
```

# Albedo

Surface albedo is stored in `state.albedo[idx]` and is updated before
the energy solve. The energy solver then treats that diagnosed value as the
surface albedo for absorbed shortwave radiation.

## Constant Scheme

When `albedo_scheme = :constant`, the code uses:

- `alpha_ice` when no snow is present
- `alpha_wet` when snow is present and the surface layer is at or above `T0`
- `alpha_dry` otherwise

## Dynamic Scheme

The dynamic scheme keeps the surface albedo between `alpha_wet` and
`alpha_dry`.

### Temperature Aging

For an existing snow surface, the code applies

```math
\alpha^\star =
\min\left(
\alpha^n,\;
\alpha^n - \Delta t_{\mathrm{days}}
\left(1.35\times 10^{-3}(T_s - T_0) + 0.0278\right)
\right),
```

followed by a clamp to `alpha_wet`. Scaling by the elapsed time
``\Delta t_{\mathrm{days}}`` keeps the aging rate consistent when a forcing
step is divided into diurnal substeps.

### Wetness Adjustment

The surface liquid-water content is diagnosed from the current pore volume:

```math
\mathrm{LWC}_{\mathrm{surf}} =
\frac{m_w / \rho_w}{m/\rho - m/\rho_i}.
```

When `LWC_surf > 0` and `max_lwc_albedo > 0`, the code linearly relaxes the
current albedo toward `alpha_wet` in proportion to
`LWC_surf / max_lwc_albedo`.

### Snowfall Refresh

Fresh snowfall brightens the surface according to

```math
\alpha^{n+} =
\min\left(
\alpha_{\mathrm{dry}},\;
\alpha^n + (\alpha_{\mathrm{dry}} - \alpha_{\mathrm{wet}})
\left(1 - e^{-\Delta m_{\mathrm{snow}}/3}\right)
\right),
```

where `Δm_snow` is the newly added snowfall mass in `kg m^-2`.

## Aging Scheme

Set `albedo=:aging` to use a snowfall-age albedo following
[Hoang et al. (2025)](https://doi.org/10.5194/cp-21-27-2025). Snowfall
resets `state.snow_age_days` to zero. During snow-free
forcing intervals the age advances by `dt_days`, and the snow albedo is

```math
\alpha_{snow} = \alpha_{firn} +
(\alpha_{fresh} - \alpha_{firn})
\exp\left(-N_{snowfall}/t^*\right).
```

The scheme reuses the existing albedo constants: `alpha_dry` is the fresh-snow
albedo, `alpha_wet` is the aged-snow/firn asymptote, and `alpha_ice` is used
without surface snow. The aging timescale defaults are
`aging_cold_timescale_days=20` and `aging_melting_timescale_days=5`.
The cold timescale is used below `T0`; the melting timescale is used at or
above `T0`. A column without surface snow uses `alpha_ice` and has zero snow
age. For time-varying surface temperature, the exponential decay is applied
incrementally over each forcing interval. This preserves the published law
when the timescale is constant and prevents snow from becoming brighter when
the surface returns below `T0`. Unlike the constant scheme, reaching the
melting point does not cause an instantaneous albedo jump.

The existing albedos and the two aging timescales can be overridden through
`BESSIModel` keywords.
For example:

```julia
model = BESSIModel(grid;
    albedo=:aging,
    alpha_dry=0.82,
    alpha_wet=0.60,
    alpha_ice=0.40,
    aging_cold_timescale_days=20.0,
    aging_melting_timescale_days=5.0,
)
```

## API

```@docs
update_surface_albedo!
```

```@docs; canonical=false
_constant_surface_albedo
_surface_liquid_water_content
_refresh_dynamic_albedo_from_snowfall!
_update_aging_surface_albedo_arrays!
_update_surface_albedo_arrays!
```
