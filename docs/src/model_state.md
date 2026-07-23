```@meta
CurrentModule = Chion
```

# Model State And Step Flow

This page summarizes the state layout and the main execution order implemented
in `src/step.jl`.

## State Model

`Simulation` owns both state snapshots:

- `sim.ref` is a deep copy of the initial state.
- `sim.now` is the state advanced by `step!` and `run!`.

`BESSIModel` and `PDDModel` are configuration objects. Use
`initial_state(model)` when a workflow needs to seed or mutate state before
constructing a `Simulation`.

For `BESSIState`, each column stores one active-layer count `N[idx]` together
with layer-wise arrays for:

| Field | Meaning | Units |
| --- | --- | --- |
| `mass` | solid snow / firn / ice mass per layer | `kg m^-2` |
| `mass_w` | retained liquid water mass per layer | `kg m^-2` |
| `density` | bulk snow density per layer | `kg m^-3` |
| `temperature` | layer temperature | `K` |
| `mass_base` | exported basal ice mass | `kg m^-2` |
| `smb_ice` | cumulative net mass forcing to the ice sheet | `kg m^-2` |
| `runoff` | cumulative liquid-water export | `kg m^-2` |
| `melt` | cumulative melted solid mass | `kg m^-2` |
| `refreezing` | cumulative refrozen liquid mass | `kg m^-2` |
| `sublimation` | cumulative sublimated mass | `kg m^-2` |
| `latent_heat_flux_sum` | time-integrated latent heat-flux diagnostic | `W m^-2 day` |
| `Tsrf` | diagnosed surface temperature | `K` |
| `albedo` | surface albedo used by the energy solver | `1` |

`PDDState` stores `snowpack_swe`, `smb_ice`, `runoff`, and `pdd_sum` as
column vectors. Its snow reservoir is capped by `PDDModel.H_snow_max`.
Refrozen water leaves that reservoir as superimposed ice, and `smb_ice`
therefore contains only snow-to-ice conversion minus ice melt. For every PDD
step, precipitation is partitioned according to
`snowfall + rainfall = Δsnowpack_swe + Δsmb_ice + Δrunoff`.

`SnowpackForcing` stores one forcing matrix per field in model-native units
(`K`, `kg m^-2 s^-1`, `W m^-2`). Its constructor also accepts user-facing
temperature and precipitation fields in Celsius and `mmWE day^-1`.

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

If a column starts the step without snow, `step!` takes the bare-ice branch:
surface albedo is set to `alpha_ice`, positive bare-ice surface energy is
converted directly into SMB loss, and the snow-column process chain is skipped.

For BESSI, `smb_ice` is the cumulative mass transferred to the ice-sheet
reservoir through basal export and bare-ice surface mass changes. Monthly
NetCDF `smb_ice` is the change in that cumulative field during the month.

## Important Defaults

- `Ntot = 15`
- `mass_max = 500 kg m^-2`
- `mass_split = 300 kg m^-2`
- `mass_min = 100 kg m^-2`
- `rho_s = 315 kg m^-3` for the constant fresh-snow-density scheme
- `T0 = 273.15 K`

## Advanced Entry Points

Most users should run through `Simulation` and `run!`. Coupled workflows can
use `init_integrator`, update `integrator.sim.forcing`, call
`sync_forcing!`, step the initialized integrator, and finish with `finalize!`.
The lower-level state `step!` methods remain available for advanced workflows
that manage state, forcing slices, and workspaces directly.

```@docs
SnowpackStepForcing
BESSIState
PDDState
initial_state
init_integrator
sync_forcing!
finalize!
step!
```
