```@meta
CurrentModule = Chion.SnowpackModel
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

## Excess-Mass Cap

After accumulation, the code applies an excess-mass cap.

For very small columns (`Ntot <= 3`), the cap is

```math
M_{\max} = 1.5\,m_{\mathrm{split}}.
```

For the standard multi-layer case, the code uses the fixed BESSI reference cap

```math
M_{\mathrm{ref}} = 15 \times 1.5 \times m_{\mathrm{split}},
```

independent of the active `Ntot` in the current run.

If the active solid mass exceeds this reference cap, the code depletes the full excess basal mass

```math
\Delta m_{\mathrm{base}} = M_{\mathrm{solid}} - M_{\mathrm{ref}}.
```

## API

```@docs
continuous_bottom_deplete!
```

```@docs; canonical=false
_split_surface_layer!
_merge_surface_layer!
_merge_bottom_layer!
_continuous_bottom_deplete!
_enforce_mass_cap!
```
