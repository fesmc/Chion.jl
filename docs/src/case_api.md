```@meta
CurrentModule = Chion
```

# Simulation API And Outputs

The supported public workflow is:

1. build a `SnowpackGrid`
2. build a model, usually `BESSIModel`
3. build `SnowpackForcing`
4. build `Simulation`
5. call `run!`

`PDDModel` and `ITMModel` can be constructed, but their `run!` methods are
placeholders until those physics paths are implemented.

## Input Units

`SnowpackForcing` accepts either model-native fields or user-facing fields:

- `air_temperature`, `snowfall_rate`, and `rainfall_rate` use native units
  (`K`, `kg m^-2 s^-1`, `kg m^-2 s^-1`).
- `air_temperature_c`, `snowfall_mm_day`, and `rainfall_mm_day` are converted
  to native units.
- `shortwave_down` is always `W m^-2`.
- scalar and time-vector forcing values are broadcast over columns.

## Outputs

`OutputOptions(save=...)` controls NetCDF output. Use `:none` or an empty
symbol vector to skip NetCDF, a variable symbol such as `:final_thickness`, a
group such as `:final`, or `"all"`.

Text summary and history CSV output are controlled separately with
`write_outputs=true`. NetCDF output requires a `SnowpackGrid` with spatial
coordinates.

## Reference

```@docs
SnowpackGrid
SnowpackForcing
BESSIModel
PDDModel
ITMModel
DynamicAlbedo
ConstantAlbedo
BESSIDensification
HTESSELDensification
ConstantFreshSnowDensity
ParameterizedFreshSnowDensity
Simulation
SimulationOptions
OutputOptions
SimulationResult
run!
load_forcing_file
```
