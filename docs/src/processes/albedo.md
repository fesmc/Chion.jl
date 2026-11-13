```@meta
CurrentModule = Chion
```

# Albedo

Surface albedo is stored in `state.albedo[idx]` and is updated before
the energy solve. The energy solver then uses it for the absorbed shortwave
radiation ``Q_{sw} = (1-\alpha)\,Q_{\mathrm{sw,down}}``.

The scheme is chosen with the `BESSIModel` keyword `albedo`:

| `albedo` | Behaviour |
| --- | --- |
| `:dynamic` (default) | prognostic BESSI albedo with temperature aging, wetness and gradual snowfall brightening |
| `:constant` | `alpha_dry`, `alpha_wet` or `alpha_ice` from the surface state |
| `:aging` | exponential decay with time since snowfall (Hoang et al., 2025) |
| `:semix` | two-band SEMIX snow albedo from grain size, dust, solar zenith angle and clouds |
| `:prescribed` | `prescribed_albedo` forcing for every column and step |

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `alpha_dry` | 0.81 | dry (fresh) snow albedo |
| `alpha_wet` | 0.70 | wet (aged) snow albedo |
| `alpha_ice` | 0.40 | bare-ice albedo |
| `max_lwc_albedo` | 0.1 | surface liquid-water content at which the albedo reaches `alpha_wet` |
| `aging_cold_timescale_days` | 20 | `:aging` timescale below ``T_0`` |
| `aging_melting_timescale_days` | 2 | `:aging` timescale at ``T_0`` |

The defaults are calibrated against MAR over Greenland. The Antarctic paper
runs use `alpha_dry=0.87`, `alpha_wet=0.82`.

In the formulas below, ``T_1`` is the temperature of the first snow layer
(`temperature[1]`), not the interface temperature `Tsrf`.

## Constant Scheme

When `albedo = :constant`, the code uses:

- `alpha_ice` when no snow is present
- `alpha_wet` when snow is present and the surface layer is at or above `T0`
- `alpha_dry` otherwise

## Dynamic Scheme

The dynamic scheme is the `BESSIModel` default. It keeps the surface albedo
between `alpha_wet` and `alpha_dry`; bare ice uses `alpha_ice`. The defaults
`alpha_dry = 0.81`, `alpha_wet = 0.70`, and `alpha_ice = 0.40` were calibrated
against MAR over Greenland in coupled runs.

The dynamic scheme follows the BESSI implementation of Zolles and Born
(2021), with temperature-dependent aging rates based on broadband albedo
observations (Aoki et al., 2003). The previous albedo is first clamped to
`[alpha_wet, alpha_dry]`. A column that was bare ice (``\alpha = 0.40``)
therefore starts from at least `alpha_wet` once any snow covers it.

### Temperature Aging

For an existing snow surface, the code applies

```math
\alpha^\star =
\min\left(
\alpha^n,\;
\alpha^n - \Delta t_{\mathrm{days}}
\left(1.35\times 10^{-3}(T_1 - T_0) + 0.0278\right)
\right),
```

followed by a clamp to `alpha_wet`, so warm snow darkens faster than cold snow. Scaling by the elapsed time
``\Delta t_{\mathrm{days}}`` keeps the aging rate consistent when a forcing
step is divided into diurnal substeps.

### Wetness Adjustment

The surface liquid-water content is diagnosed from the current pore volume:

```math
\mathrm{LWC}_{\mathrm{surf}} =
\frac{m_w / \rho_w}{m/\rho - m/\rho_i}.
```

When `LWC_surf > 0` and `max_lwc_albedo > 0`, define
`r = clamp(LWC_surf / max_lwc_albedo, 0, 1)`. The wetness adjustment is

```math
\alpha^{n+1} = \alpha_{\mathrm{wet}} +
(\alpha^\star - \alpha_{\mathrm{wet}})(1-r)^{\Delta t_{\mathrm{days}}}.
```

This preserves the original one-day adjustment and makes wetness relaxation
independent of subdivision at fixed wetness, in isolation from temperature
aging. Zero elapsed time applies no wetness adjustment; `r = 1` reaches
`alpha_wet` for any positive elapsed time. Combined aging, snowfall, and
evolving column conditions can still introduce timestep dependence.

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

where `Δm_snow` is the newly added snowfall mass in `kg m^-2`. Small events
brighten the surface gradually, whereas events of more than about
10 kg m⁻² restore an albedo close to that of fresh snow. The refresh happens
during accumulation. The `:constant` and `:aging` schemes instead reset the
albedo to `alpha_dry` after any snowfall.

## Aging Scheme

Set `albedo=:aging` to use a snowfall-age albedo following
[Hoang et al. (2025)](https://doi.org/10.5194/cp-21-27-2025). Snowfall
resets `state.snow_age_days` to zero and the albedo to `alpha_dry`. During
snow-free forcing intervals the age advances by `dt_days`, and the albedo
relaxes towards the aged-snow value step by step,

```math
\alpha^{n+1} = \alpha_{\mathrm{wet}} +
(\alpha^n - \alpha_{\mathrm{wet}})
\exp\left(-\Delta t_{\mathrm{days}}/\tau\right),
```

which for constant ``\tau`` is the published
``\alpha = \alpha_{\mathrm{firn}} + (\alpha_{\mathrm{fresh}} - \alpha_{\mathrm{firn}})\exp(-t/\tau)``
with ``t`` the time since snowfall.

The scheme reuses the existing albedo constants: `alpha_dry` is the fresh-snow
albedo, `alpha_wet` is the aged-snow/firn asymptote, and `alpha_ice` is used
without surface snow. The aging timescale defaults are
`aging_cold_timescale_days=20` and `aging_melting_timescale_days=2`.
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
    aging_melting_timescale_days=2.0,
)
```

## SEMIX Scheme

`albedo=:semix` uses the spectral snow albedo of SEMIX (Willeit et al., 2024).
It computes visible and near-infrared albedos for direct and diffuse light from

- the snow grain size, which grows with temperature and shrinks with snowfall
  (`semix_snow_grain_fresh`, `semix_snow_grain_old`, `semix_snow_0`,
  `semix_snow_1`, `semix_dT_age`, `semix_f_age_t`)
- the dust concentration in the snow, from optional `dust_deposition` forcing
  and amplified during melt of the seasonal snowpack
  (`semix_dust_con_scale`, `semix_w_snow_dust`, state `w_snow_max`)
- the cosine of the solar zenith angle (`coszm` forcing, or the daily mean
  from latitude and season)
- the subgrid orographic roughness (`z_sur_std` forcing; Dang scheme only)

following Dang et al. (2015) (`semix_snow_albedo=:dang`, default) or Warren
and Wiscombe (1980) (`:warren_wiscombe`). The bands are combined into a
broadband albedo with the cloud fraction (`cloud` forcing) and the visible
fraction `semix_frac_vu`. Over bare ice the scheme uses `alpha_ice`, or the
optional `prescribed_ice_albedo` forcing.

## Prescribed Scheme

`albedo=:prescribed` takes the surface albedo from the `prescribed_albedo`
forcing, which must be present for every column and step. It also applies over
bare ice.

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
