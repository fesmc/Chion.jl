```@meta
CurrentModule = Chion.SnowpackModel
```

# Accumulation and Ablation

## Accumulation

`apply_accumulation!` applies snowfall and rainfall forcing to the surface layer and enforces the dynamic layer-mass constraints.

### Surface Mass Update

With snowfall rate ``P_{snow}`` and rainfall rate ``P_{rain}`` in ``[\mathrm{kg\,m^{-2}\,s^{-1}}]``, and timestep ``\Delta t`` in seconds:

```math
\Delta m_{snow} = P_{snow}\,\Delta t,\qquad
\Delta m_{rain} = P_{rain}\,\Delta t.
```

If the column is empty (``N=0``), a new layer is created only when ``P_{snow}>0``. Rain alone does not create a snow layer.

Snowfall is added to layer 1 (`column.mass[1]`) and its density is mixed with fresh-snow density ``\rho_\mathrm{fresh}``.
The fresh-snow density is taken from the HTESSEL model:

```math
\rho_\mathrm{fresh} = a + b\,(T_\mathrm{air} - T_0) + c\,\sqrt{V},
```

with the constants ``a=109\,\mathrm{kgm^{-3}}``, ``b=6\,\mathrm{kgm^{-3}}``, ``c=26\,\mathrm{kgm^{-3.5}s^{0.5}}`` and the optional wind speed ``V`` in ``[\mathrm{m\,s^{-1}}]`` (fallback ``V=5\,\mathrm{m/s}``). 
The snow density is limited to a minimum of ``50\,\mathrm{kgm^{-3}}``.

![Fresh-snow density vs wind speed](../assets/fresh_snow_density_vs_wind.png)

The new snow density is given by
```math
\rho_1^{new} =
\frac{m_1^{old}+\Delta m_{snow}}
{m_1^{old}/\rho_1^{old} + \Delta m_{snow}/\rho_\mathrm{fresh}}.
```

Rainfall is added to liquid water mass in the surface layer (`column.mass_w[1]`).

### Layer-Structure Rules

After mass addition, the model enforces layer limits:

1. If surface mass exceeds `mass_max`, split the surface layer (`split_surface_layer!`).
2. If all layers are already active (`N == Ntot`), first merge the two bottom layers (`merge_bottom_layer!`) to make space.
3. If surface mass drops below `mass_min`, merge with layer 2 (`merge_surface_layer!`).
4. After the layer-structure updates, the model computes a total-column excess mass

```math
d m = \sum_{k=1}^{N} m_k - 15 \, m_{\mathrm{split}} \, 1.5.
```

If `dm > 0`, that excess is removed from the bottom upward with
`continuous_bottom_deplete!`.

### Basal Firn-to-Ice Transfer

The basal transfer now follows the logic of the original BESSI routine for the depletion step, and it is triggered
from a total-column excess-mass criterion.

Given a requested basal handoff `d_m_in` in ``\mathrm{kg\,m^{-2}}]``, the model removes that solid mass from the bottom upward:

1. Start at the deepest active layer.
2. If `d_m_in` exceeds the solid mass in that layer, remove the whole layer.
3. Otherwise remove only `d_m_in` from that layer.
4. On partial removal, remove liquid water proportionally:

```math
\Delta m_w = d_m \frac{m_w}{m}.
```

5. Add removed solid mass to `mass_base`.
6. Add removed liquid water to `runoff`.

The present trigger uses the 15-layer reference column as implemented in BESSI, even if the
active column uses a different total number of layers `Ntot`. In other words, the current code
uses the same excess-mass reference as BESSI:

```math
d m = M_{\mathrm{column}} - 15 \, m_{\mathrm{split}} \, 1.5.
```

The basal handoff is now mass-conserving and depth-aware:

## Ablation

`apply_melt!` removes melt mass from the top downward by converting snow mass into liquid water in the active layer(s).

Given a requested melt amount ``m_{melt}``:

1. Clamp to nonnegative melt.
2. Remove melt from layer 1 up to available mass.
3. Add the melted amount to layer liquid water mass.
4. Continue into deeper layers only when the current surface layer is fully depleted.

If a layer is depleted and a deeper snow layer exists, its liquid water is transferred to the next layer before removing the depleted surface layer. If the surface layer remains but falls below `mass_min` (and another layer exists), it is merged with the next layer.

If no receiving snow layer exists, remaining liquid water is routed to runoff:

```math
\mathrm{runoff} \leftarrow \mathrm{runoff} + m_{w,\mathrm{surface}}.
```

So `apply_melt!` is a phase-change/removal routine; most runoff is produced later by `go_percolation!` when liquid water exceeds retention limits.

Surface temperature state is synchronized after melt (`column.Tsrf = column.temperature[1]` when snow exists, else `T0`).

## Interaction in `step!`

In one timestep (`step!`), accumulation/ablation-related ordering is:

1. `apply_accumulation!`
2. densification (when sufficient active layers exist)
3. energy solve (`go_energy_flux!`)
4. if melt is diagnosed, `apply_melt!`
5. percolation (`go_percolation!`)
6. refreezing (`go_refreezing!`)

## API

```@docs; canonical=false
apply_melt!
```

```@docs; canonical=false
apply_accumulation!
```

```@docs; canonical=false
continuous_bottom_deplete!
```
