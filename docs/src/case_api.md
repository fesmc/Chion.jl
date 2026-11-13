```@meta
CurrentModule = Chion
```

# Simulation API And Outputs

For runs from forcing files with spin-up and outputs, use an `Experiment`
([Running Experiments](experiments.md)). This page describes the lower-level
`Simulation` API underneath it, for forcing held in memory, step-by-step
coupling to other models, and custom loops:

1. build a `SnowpackGrid`
2. build `BESSIModel`, `PDDModel`, or `ITMModel` (or `build_model(:bessi, grid; ...)`)
3. build `SnowpackForcing` (see its docstring below for the optional fields)
4. build `Simulation`
5. call `run!`

`Simulation` owns two model-state objects: `sim.ref` is the reference state and
`sim.now` is the evolving state. Model structs are configuration-only. By
default a `Simulation` runs the forcing once (`years=1`) on the CPU and writes
no NetCDF; pass `write_netcdf=true`, `netcdf_variables`, `output_dir` and
`backend=:gpu` to change this.

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
`finalize!(integrator)`. `set_active_mask!(integrator, mask)` restricts stepping
to a subset of columns, and `yearly_step!(integrator)` advances a BESSI
simulation by one forcing year and returns the annual coupling fields (SMB and
mean surface temperature).

`BESSIModel` resolves the layered mass and energy balance. Its defaults are
the setup calibrated against MAR over Greenland for daily forcing (see
[Model State And Step Flow](model_state.md#Important-Defaults)). The default
diurnal substeps require `latitude_deg` in the forcing; set
`diurnal_shortwave_substeps=false` when latitude is unavailable or the forcing
is already sub-daily. `PDDModel` uses a
bulk snow reservoir with configurable `ddf_snow`, `ddf_ice`,
`refreezing_fraction`, `temperature_sigma`, `H_snow_max`, and
`pdd_method`. The default `pdd_method=:simple` uses positive mean
temperature; `pdd_method=:pism` applies the Calov-Greve expectation integral
at every timestep.

`ITMModel` ports the Fortran insolation-temperature-melt scheme. It requires
explicit `latitude_deg`, `surface_height`, `ice_thickness` (m), and
`annual_pdd` (K day) in `SnowpackForcing`. `q_sw_net`, when supplied, is its
insolation driver; otherwise ITM uses `shortwave_down`. As in the original
ITM, this insolation is top-of-atmosphere radiation that the model attenuates
with its elevation-dependent transmissivity; supplying surface shortwave (for
example MAR `SWD`) attenuates it twice and strongly underestimates melt. It
preserves the
Fortran scheme's per-step `mmWE day^-1` fields (`smb`, `smbi`, `melt`,
`runoff`, `refreezing`) and cumulative `mmWE` fields.

## Input Units

`SnowpackForcing` accepts either model-native fields or user-facing fields:

- `air_temperature`, `snowfall_rate`, and `rainfall_rate` use native units
  (`K`, `kg m^-2 s^-1`, `kg m^-2 s^-1`).
- `air_temperature_c`, `snowfall_mm_day`, and `rainfall_mm_day` are converted
  to native units.
- `shortwave_down` is always `W m^-2`.
- `ice_thickness` is metres and `annual_pdd` is K day; both are required by
  `ITMModel`.
- scalar and time-vector forcing values are broadcast over columns.

`load_forcing_file` can select spatial columns directly with `mask_name` and
`mask_threshold`. Selection depends only on that mask variable. Each field is
read once straight into the selected columns; dimensions are recognised by the
names given as `x_name`/`y_name` (or `x`/`y`) and a time dimension, and any
other dimension (vertical level, sector) is read at its first index.
`read_forcing_columns` exposes the same reader for additional variables.

A host whose climate supplies only monthly air temperature and precipitation
can derive the rest with `shortwave_down_from_transmissivity` (shortwave from
a per-domain direct transmissivity; the longwave cloud proxy then diagnoses its
own cloudiness),
`wind_speed_from_geometry` (wind from elevation and slope),
`temperature_variability_sigma` and `red_noise_monthly` (daily temperature
variability about a smooth series). Chion does not call them itself.

## Outputs

`Simulation(...; write_netcdf=false)` skips NetCDF output. Use
`netcdf_variables` to request a state field such as `:thickness`, `:all`, or
`:monthly`.

For `BESSIModel`, `:monthly` writes monthly surface SMB (`surface_smb`:
precipitation − runoff + vapour exchange), `smb_ice`, runoff, melt,
refreezing and sublimation changes, plus albedo and the monthly mean surface
energy balance (`W m^-2`, positive into the surface): `shortwave_down`,
`shortwave_net`, `longwave_down`, `longwave_net`, `sensible_heat_flux`,
`latent_heat_flux` and `rain_heat_flux`.

For `PDDModel`, `:all` writes `snowpack_swe`, `smb_ice`, `runoff`, and
`pdd_sum` every forcing step. `:monthly` writes month-end snowpack SWE and
monthly changes in the three cumulative fields. Monthly output for both models
is buffered and written to NetCDF in chunks after stepping.

For `ITMModel`, `:all` writes its Fortran-compatible state and budget fields
at every forcing step. `:monthly` writes month-end snapshots of those fields.

NetCDF outputs include an unlimited, CF-style `t` record coordinate in days
since `1970-01-01 00:00:00`, using the `proleptic_gregorian` calendar. Full-state output
uses the forcing timestamps, while monthly output uses the timestamp of the
last forcing record in each month. For simulations that repeat the forcing for
multiple years, timestamps are advanced by one calendar year for each cycle.
PDD `smb_ice` is ice-facing, consistent with BESSI: refrozen water and snow
above `H_snow_max` are transferred to ice, while ice melt is negative SMB.
Seasonal snow retained in `snowpack_swe` is not credited to the ice sheet.

NetCDF output requires a `SnowpackGrid` with spatial coordinates.

Use `backend=:threads` for the CPU and `backend=:gpu` with a functional CUDA
environment. CUDA.jl is a weak dependency: `using CUDA` loads the GPU backend
(the `ChionCUDAExt` package extension). Without it, `backend=:gpu` stops with an
error that asks for `using CUDA`.

## Host Forcing Helpers

These functions build forcing fields for runs without a regional climate model
(for example, forcing from an ice-sheet model's climate). `clear_sky_transmissivity`
needs the physical constants of a model, `model.c`:

```julia
model = BESSIModel(grid)
tau_cs = clear_sky_transmissivity(model.c, 2000.0, 300.0)   # height (m), TOA SW (W m⁻²)
```

`PDDModel` and `ITMModel` accept `number_type` like `BESSIModel`
(`number_type=Float32` for single precision).

## Reference

```@docs
SnowpackGrid
SnowpackForcing
update_air_pressure!
build_model
set_active_mask!
yearly_step!
clear_sky_transmissivity
shortwave_down_from_transmissivity
wind_speed_from_geometry
temperature_variability_sigma
red_noise_monthly
BESSIModel
PDDModel
ITMModel
Simulation
RunOptions
SimulationResult
BESSIState
PDDState
ITMState
initial_state
init_integrator
SimulationIntegrator
sync_forcing!
finished
finalize!
run!
load_forcing_file
read_forcing_columns
```
