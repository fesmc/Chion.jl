```@meta
CurrentModule = Chion
```

# Simulation API And Outputs

The supported public workflow is:

1. build a `SnowpackGrid`
2. build `BESSIModel` or `PDDModel`
3. build `SnowpackForcing`
4. build `Simulation`
5. call `run!`

`Simulation` owns two model-state objects: `sim.ref` is the reference state and
`sim.now` is the evolving state. Model structs are configuration-only.

For coupling workflows, initialize an explicit stepper and update its forcing
between steps:

```julia
integrator = init_integrator(sim)
step!(integrator)                 # one scheduled forcing step
step!(integrator, 10)             # ten scheduled forcing steps
integrator.sim.forcing.air_temperature .= 265.15
sync_forcing!(integrator)
step!(integrator, 1.0, true)      # one externally supplied 1-day step
result = finalize!(integrator)
```

`run!(sim)` is a wrapper over `init_integrator`, `run!(integrator)`, and
`finalize!(integrator)`.

`BESSIModel` resolves the layered mass and energy balance. `PDDModel` uses a
bulk snow reservoir with configurable `ddf_snow`, `ddf_ice`,
`refreezing_fraction`, `temperature_sigma`, `H_snow_max`, and
`pdd_method`. The default `pdd_method=:simple` uses positive mean
temperature; `pdd_method=:pism` applies the Calov-Greve expectation integral
at every timestep.

## Input Units

`SnowpackForcing` accepts either model-native fields or user-facing fields:

- `air_temperature`, `snowfall_rate`, and `rainfall_rate` use native units
  (`K`, `kg m^-2 s^-1`, `kg m^-2 s^-1`).
- `air_temperature_c`, `snowfall_mm_day`, and `rainfall_mm_day` are converted
  to native units.
- `shortwave_down` is always `W m^-2`.
- scalar and time-vector forcing values are broadcast over columns.

`load_forcing_file` can select spatial columns directly with `mask_name` and
`mask_threshold`. Selection depends only on that mask variable.

## Outputs

`Simulation(...; write_netcdf=false)` skips NetCDF output. Use
`netcdf_variables` to request a state field such as `:thickness`, `:all`, or
`:monthly`.

For `BESSIModel`, `:monthly` writes monthly `smb_ice`, runoff, melt,
refreezing and sublimation changes, plus monthly mean latent heat flux and
albedo.

For `PDDModel`, `:all` writes `snowpack_swe`, `smb_ice`, `runoff`, and
`pdd_sum` every forcing step. `:monthly` writes month-end snowpack SWE and
monthly changes in the three cumulative fields. Monthly output for both models
is buffered and written to NetCDF in chunks after stepping.

PDD `smb_ice` is ice-facing, consistent with BESSI: refrozen water and snow
above `H_snow_max` are transferred to ice, while ice melt is negative SMB.
Seasonal snow retained in `snowpack_swe` is not credited to the ice sheet.

NetCDF output requires a `SnowpackGrid` with spatial coordinates.

`backend=:cpu` is an alias for `backend=:threads`; `backend=:gpu` requires a
functional CUDA environment.

## Reference

```@docs
SnowpackGrid
SnowpackForcing
BESSIModel
PDDModel
Simulation
RunOptions
SimulationResult
BESSIState
PDDState
initial_state
init_integrator
SimulationIntegrator
sync_forcing!
finished
finalize!
run!
load_forcing_file
```
