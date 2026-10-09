```@meta
CurrentModule = Chion
```

# Layer Structure And Basal Transfer

These routines manage the adaptive layer stack and the basal export of excess
mass. They are called from accumulation and melt and define much of the model's
mass bookkeeping.

## Surface Split And Merge

The surface layer is split when `mass[1] > mass_max`.

- layer 1 keeps `surface_mass - mass_split`
- the new layer 2 receives `mass_split`
- liquid water is partitioned in proportion to the split mass
- density and temperature are copied into both layers

The surface layer is merged or rebalanced when `mass[1] < mass_min`.

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

The reference depth is the depth of a 15-layer column at
300 kg m^-3 density:

```math
H_{\mathrm{ref}} =
\frac{15 \times 1.5 \times m_{\mathrm{split}}}{300\,\mathrm{kg\,m^{-3}}}.
```

This cap is independent of the active `Ntot` in the current run, so a reduced
layer count can still represent the same maximum physical depth.

If the active solid snow depth exceeds this reference depth, the code depletes
only enough basal mass to remove the excess depth:

```math
\Delta H_{\mathrm{base}} = H_{\mathrm{solid}} - H_{\mathrm{ref}}.
```

## Near-Surface Layers

`near_surface_layer_max_thicknesses_m` (default `(0.02, 0.05, 0.10, 0.30)` m)
fixes the maximum thickness of the four uppermost layers. After accumulation
and again after melt, sublimation, compaction and refreezing, the column is
re-meshed conservatively so that the surface energy balance always sees a thin
top cell. Without these limits the surface layer can hold up to `mass_max`.

## Ice Substrate

`ice_substrate_layers` (default 5) adds a fixed-geometry layer of glacier ice
below the snow/firn column. Its top layer is `ice_substrate_top_thickness_m`
(default 0.05 m) thick and thicknesses double downward (1.55 m in total for the
defaults). The substrate is thermal only: it is solved in the same implicit heat
equation as the snow layers with an insulated base, stores cold content and
forms the surface when no snow is left. Bare ice must therefore be re-warmed to
the melting point before it melts; melt and rain on bare ice run off. Set
`ice_substrate_layers=0` for the original bare-ice treatment at the melting point.

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
```
