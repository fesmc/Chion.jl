```@meta
CurrentModule = Chion.SnowpackModel
```

# Percolation

`go_percolation!` redistributes liquid water vertically through active snow layers and returns the runoff generated during that call.

## Model Formulation

For each layer ``i``, the model computes liquid water content (LWC) as pore-volume saturation:

```math
\mathrm{LWC}_i =
\frac{m_{w,i}}
{m_{s,i}\,\rho_w\left(\frac{1}{\rho_i^{\mathrm{snow}}}-\frac{1}{\rho_i^{\mathrm{ice}}}\right)}
```

where:

- ``m_{w,i}``: liquid water mass in layer ``i [\mathrm{kg\,m^{-2}}]``
- ``m_{s,i}``: snow mass in layer ``i`` ``[\mathrm{kg\,m^{-2}}]``
- ``\rho_w``: water density ``[\mathrm{kg\,m^{-3}}]``
- ``\rho_i^{\mathrm{snow}}``: bulk snow density in layer ``i`` ``[\mathrm{kg\,m^{-3}}]``
- ``\rho_i^{\mathrm{ice}}``: ice density ``[\mathrm{kg\,m^{-3}}]``

If ``\mathrm{LWC}_i > \mathrm{LWC}_{\max}``, excess liquid water percolates downward:

```math
\Delta m_{w,i} =
\left(\mathrm{LWC}_i - \mathrm{LWC}_{\max}\right)\,
\rho_w\,m_{s,i}\left(\frac{1}{\rho_i^{\mathrm{snow}}}-\frac{1}{\rho_i^{\mathrm{ice}}}\right)
```

and the layer is clipped to ``\mathrm{LWC}_{\max}``.

## Dense-Layer Rule

If a layer is nearly ice (``\rho_i^{\mathrm{snow}} > \rho_i^{\mathrm{ice}} - \rho_{\mathrm{tol}}``), it cannot retain liquid water in this scheme. All liquid water in that layer is transferred downward immediately (or to runoff if no receiving snow layer exists).

## Routing Logic

For each layer from top to bottom:

1. Compute percolating amount based on dense-layer or excess-LWC rule.
2. If the next layer exists and contains snow mass, add percolating water there.
3. Otherwise, add percolating water to runoff.

The loop follows the original BESSI-style control flow and stops at the first empty layer.

## API


```@docs; canonical=false
go_percolation!

```

## Notes

- `max_lwc` is dimensionless pore-space saturation (default `0.05`).
- Returned runoff has units ``[\mathrm{kg\,m^{-2}}]``.
- The `SnowpackColumn` method also accumulates returned runoff into `column.runoff`.
