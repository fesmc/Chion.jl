```@meta
CurrentModule = Chion
```

# Layer Structure And Basal Transfer

These routines manage the adaptive layer stack and the basal export of excess
mass. They are called from accumulation and melt and define much of the model's
mass bookkeeping.

The column is discretised on a mass-following grid with at most `Ntot`
(default 15) layers. Each layer carries solid mass ``m_i``, liquid water
``m_{w,i}`` (kg m⁻²), density ``\rho_i`` and temperature ``T_i``, so the layer
thickness is ``\Delta z_i = m_i/\rho_i``. Mass leaving the bottom of the column
(bottom merges above ``\rho_i`` and the depth cap) is transferred to the
underlying ice and forms the snow-to-ice term of the ice-sheet SMB.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `Ntot` | 15 | maximum number of snow/firn layers |
| `mass_max` | 500 kg m⁻² | split threshold |
| `mass_split` | 300 kg m⁻² | mass moved into the new layer by a split |
| `mass_min` | 100 kg m⁻² | surface-merge threshold (inactive with near-surface layers) |
| `near_surface_layer_max_thicknesses_m` | (0.02, 0.05, 0.10, 0.30) m | target thicknesses of the top four layers |
| `ice_substrate_layers` | 5 | thermal ice layers below the snow |
| `ice_substrate_top_thickness_m` | 0.05 m | thickness of the top substrate layer (doubling downward) |

## Surface Split And Merge

The surface layer is split when `mass[1] > mass_max`.

- layer 1 keeps `surface_mass - mass_split`
- the new layer 2 receives `mass_split`
- liquid water is partitioned in proportion to the split mass
- density and temperature are copied into both layers

The surface layer is merged or rebalanced when `mass[1] < mass_min`. This
merge is switched off with the default near-surface layers (whenever the first
entry of `near_surface_layer_max_thicknesses_m` is finite), because the
remeshing controls the surface layers.

- if the top two layers together exceed `2 * mass_split`, the surface is
  topped up to exactly `mass_split`
- otherwise the two top layers are merged completely

## Bottom Merge And Slot Creation

When all `Ntot` layers are already active and the surface must split, the code
creates a free slot at the bottom:

- for `Ntot > 2`, it merges the deepest two active layers
- for very small columns, it depletes the bottom instead

If a bottom merge would imply density above `rho_i`, the excess mass is
exported to `mass_base` and `smb_ice`.

## Continuous Basal Depletion

`continuous_bottom_deplete!` removes a requested solid mass from the bottom
upward. Partial removal from a layer also removes liquid water in proportion to
the removed solid fraction:

```math
\Delta m_w = \Delta m \frac{m_w}{m}.
```

Removed solid mass is accumulated in `mass_base` and `smb_ice`; removed liquid
water is sent to runoff.

## Snow-Depth Cap

After accumulation, the code applies a snow-depth cap.

The reference depth is the fixed constant
``H_{\mathrm{ref}} = 22.5`` m (`BESSI_REFERENCE_SNOW_DEPTH_M`, the depth of
15 layers of 450 kg m⁻² at 300 kg m⁻³ in the original BESSI). It is independent
of `Ntot` and `mass_split`, so a reduced layer count can still represent the
same maximum physical depth.

If the active solid snow depth exceeds this reference depth, the code depletes
only enough basal mass to remove the excess depth (solid mass goes to the ice,
liquid water to runoff):

```math
\Delta H_{\mathrm{base}} = H_{\mathrm{solid}} - H_{\mathrm{ref}}.
```

## Near-Surface Layers

The coarse BESSI layers of several hundred kg m⁻² do not resolve the strong
near-surface temperature gradients that control the surface temperature and
the refreezing of meltwater. `near_surface_layer_max_thicknesses_m` (default
`(0.02, 0.05, 0.10, 0.30)` m) therefore sets **target** thicknesses for the
four uppermost layers. After accumulation and again at the end of each step
(after melt, sublimation, compaction and refreezing), material exceeding these
thicknesses is passed to the layer below, and layers thinner than their target
are refilled from below (`_fill_near_surface_layer_thicknesses!`). The
remeshing conserves solid mass, liquid water, volume and sensible heat and adds
no physical process.

An `Inf` entry disables that level. All four `Inf` recover the original BESSI
layering, in which the surface layer can hold up to `mass_max` and the
`mass_min` surface merge is active.

Layers below the refined near-surface layers (from layer 5 down) are split
again when they exceed `mass_max` (`_split_deep_layers!`), splitting off
`mass_split`. A full column first merges its two deepest layers and exports
any excess above ``\rho_i``.

## Ice Substrate

`ice_substrate_layers` (default 5) adds a fixed-geometry layer of glacier ice
below the snow/firn column. Its top layer is `ice_substrate_top_thickness_m`
(default 0.05 m) thick and thicknesses double downward (1.55 m in total for the
defaults). The substrate is thermal only: it is solved in the same implicit heat
equation as the snow layers with an insulated base, stores cold content and
forms the surface when no snow is left. It is initialised at
`temperature_init`. Bare ice must therefore be re-warmed to the melting point
before it melts. Without the substrate, bare ice would be held at the melting
point and negative energy fluxes discarded, which overestimates melt on
sub-daily time scales. The substrate is impermeable: melt and rain on bare ice
run off. Set `ice_substrate_layers=0` for the original bare-ice treatment at
the melting point.

## API

```@docs
continuous_bottom_deplete!
```

```@docs; canonical=false
_split_surface_layer!
_merge_surface_layer!
_merge_bottom_layer!
_continuous_bottom_deplete!
_enforce_snow_depth_cap!
_fill_near_surface_layer_thicknesses!
_split_deep_layers!
```
