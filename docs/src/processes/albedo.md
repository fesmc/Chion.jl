```@meta
CurrentModule = Chion
```

# Albedo

Surface albedo is stored in `domain.albedo_dynamic[idx]` and is updated before
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
\alpha^n - \left(1.35\times 10^{-3}(T_s - T_0) + 0.0278\right)
\right),
```

followed by a clamp to `alpha_wet`.

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

## API

```@docs
update_surface_albedo!
```

```@docs; canonical=false
_constant_surface_albedo
_surface_liquid_water_content
_refresh_dynamic_albedo_from_snowfall!
_update_surface_albedo_arrays!
```
