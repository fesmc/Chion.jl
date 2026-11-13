```@meta
CurrentModule = Chion
```

# Percolation

`go_percolation!` moves liquid water downward with a bucket scheme. Each layer
retains water up to a volumetric content of `max_lwc` = 0.1 of its pore volume.
The pass runs once from the surface downward, before refreezing, so water in
excess of the holding capacity passes cold layers within the same step. Water
leaving the deepest layer runs off, because the underlying ice is impermeable.
Lateral runoff of water ponding on ice layers within the firn is not
represented.

## Model Formulation

For each active layer ``i``, the code first computes pore volume as

```math
\phi_i =
\frac{m_{s,i}}{\rho_i^{\mathrm{snow}}}
- \frac{m_{s,i}}{\rho_i^{\mathrm{ice}}}.
```

If ``\phi_i > 0``, liquid-water content is

```math
\mathrm{LWC}_i =
\frac{m_{w,i}}{\rho_w \phi_i}.
```

If ``\mathrm{LWC}_i > \mathrm{LWC}_{\max}``, the excess liquid water is

```math
\Delta m_{w,i} =
\left(\mathrm{LWC}_i - \mathrm{LWC}_{\max}\right)\rho_w\phi_i
```

and is routed to the next active layer or to runoff.

## Pore-Collapse Routing

If the pore volume is nonpositive, the current implementation routes all liquid
water in that layer downward immediately, or to runoff if there is no deeper
active snow layer.

Layers with zero solid mass cannot retain water: their liquid water goes
directly to runoff.

## Notes

- `max_lwc = 0.1` is fixed inside `step!`; it is a keyword only of the
  standalone `go_percolation!`
- returned runoff has units ``\mathrm{kg\,m^{-2}}``

## API

```@docs
go_percolation!
```

```@docs; canonical=false
_go_percolation!
```
