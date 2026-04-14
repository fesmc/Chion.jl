```@meta
CurrentModule = Chion
```

# Model State And Step Flow

This page summarizes the state layout and the main execution order implemented
in `src/step.jl`.

## State Model

Each column stores one active-layer count `N[idx]` together with layer-wise
arrays for:

| Field | Meaning | Units |
| --- | --- | --- |
| `mass` | solid snow / firn / ice mass per layer | `kg m^-2` |
| `mass_w` | retained liquid water mass per layer | `kg m^-2` |
| `density` | bulk snow density per layer | `kg m^-3` |
| `temperature` | layer temperature | `K` |
| `mass_base` | exported basal ice mass | `kg m^-2` |
| `smb_ice` | ice-sheet SMB contribution | `kg m^-2` |
| `runoff` | cumulative liquid-water export | `kg m^-2` |
| `Tsrf` | diagnosed surface temperature | `K` |
| `snow_cover` | diagnosed snow-cover fraction | `1` |
| `albedo_dynamic` | surface albedo used by the energy solver | `1` |

`SnowpackStepForcing` uses model-native units (`K`, `kg m^-2 s^-1`, `W m^-2`).
The higher-level case API accepts user-facing forcing in Celsius and
`mmWE day^-1`, then converts it into the native forcing arrays stored in
the case definition.

## Step Ordering

For a snow-covered column, `step!` currently executes the processes in this
order:

1. accumulation and rainfall input
2. albedo update
3. densification
4. energy solve
5. melt, when the energy solve diagnoses melt energy
6. percolation
7. HTESSEL liquid-water compaction, when that scheme is active and liquid water remains
8. refreezing
9. snow-cover refresh

If a column starts the step without snow, `step!` takes the bare-ice branch:
surface albedo is set to `alpha_ice`, positive bare-ice surface energy is
converted directly into SMB loss, and the snow-column process chain is skipped.

## Important Defaults

- `Ntot = 15`
- `mass_max = 500 kg m^-2`
- `mass_split = 300 kg m^-2`
- `mass_min = 100 kg m^-2`
- `rho_s = 315 kg m^-3` for the constant fresh-snow-density scheme
- `T0 = 273.15 K`

## Public Entry Points

```@docs
SnowpackPhysicalConstants
SnowpackDomain
SnowpackStepForcing
SnowpackStepFields
step!
```
