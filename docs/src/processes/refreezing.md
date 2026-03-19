```@meta
CurrentModule = Chion.SnowpackModel
```

# Refreezing

`go_refreezing!` converts liquid water to ice using layer cold content and updates temperature, density, solid mass, and liquid water mass consistently.

## Model Formulation

For each layer ``i`` with snow and liquid water:

```math
Q_{cold,i} = (T_0 - T_i)\,c_i\,m_{s,i}
```

```math
Q_{lat,i} = m_{w,i}\,L_m
```

where:

- ``m_{s,i}``: snow/ice mass ``[\mathrm{kg\,m^{-2}}]``
- ``m_{w,i}``: liquid water mass ``[\mathrm{kg\,m^{-2}}]``
- ``T_i``: layer temperature ``[\mathrm{K}]``
- ``T_0``: melting point ``[\mathrm{K}]``
- ``c_i``: ice heat capacity ``[\mathrm{J\,kg^{-1}\,K^{-1}}]``
- ``L_m``: latent heat of fusion ``[\mathrm{J\,kg^{-1}}]``

### Case 1: Partial refreezing (``Q_{cold} < Q_{lat}``)

Only part of liquid water freezes:

```math
\Delta m_{ice} = \frac{Q_{cold}}{L_m}
```

Then:

- ``T_i \leftarrow T_0``
- ``m_{s,i} \leftarrow m_{s,i} + \Delta m_{ice}``
- ``m_{w,i} \leftarrow m_{w,i} - \Delta m_{ice}``
- density is rescaled with mass gain

### Case 2: Complete refreezing (``Q_{cold} \ge Q_{lat}``)

All liquid freezes:

- ``m_{s,i} \leftarrow m_{s,i} + m_{w,i}``
- ``m_{w,i} \leftarrow 0``
- layer temperature is recomputed by energy balance:

```math
T_i \leftarrow
\frac{m_{w,i}L_m/c_i + m_{w,i}T_0 + T_i m_{s,i}}
{m_{w,i}+m_{s,i}}
```

## Outputs

The function returns:

- `refrozen_mass`: total refrozen mass ``[\mathrm{kg\,m^{-2}}]``
- `released_latent_heat`: latent heat released ``[\mathrm{J\,m^{-2}}]``

## API

```@docs; canonical=false
go_refreezing!
```
