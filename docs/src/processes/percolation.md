```@meta
CurrentModule = Chion.SnowpackModel
```

# Percolation

`go_percolation!` redistributes retained liquid water downward until each layer
is at or below the configured retention threshold.

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

Layers with zero solid mass are treated the same way: their stored liquid water
is routed onward because they cannot retain it.

## Notes

- the current default is `max_lwc = 0.1`
- returned runoff has units ``\mathrm{kg\,m^{-2}}``
- `rho_i_tol` is still accepted by the function signatures but is not used by
  the current implementation

## API

```@docs
go_percolation!
```

```@docs; canonical=false
_go_percolation!
```
