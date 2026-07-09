# Chion.jl

Chion is a fast, intermediate-complexity snowpack mass and energy balance model.
It supports single-column and gridded snowpack runs, CPU execution,
single-GPU execution through `CUDA.jl`, and optional NetCDF output.

The current public workflow is:

1. Build a `SnowpackGrid`.
2. Build a `BESSIModel` or `PDDModel`.
3. Build a `SnowpackForcing`, or load one with `load_forcing_file`.
4. Build a `Simulation`.
5. Execute it with `run!`.

This model is still work in progress. Not everything is validated yet, and the
API may still change.

## Features

- Snowpack mass and energy balance core with layer-based state
- Simulation-first public API for manual and file-backed forcing workflows
- CPU and single-GPU execution backends
- Per-year BESSI history and optional NetCDF output
- Pluto notebooks and script entry points for demos and larger runs

## Installation

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then load the package with:

```julia
using Chion
```

GPU runs require a working CUDA environment that `CUDA.jl` can use. NetCDF
loading and output use `NCDatasets.jl`.

## Quick Start

This is the smallest supported end-to-end workflow using forcing arrays.
High-level forcing inputs use Celsius and `mmWE/day`.

```julia
using Chion

grid = SnowpackGrid(1)

model = BESSIModel(grid;
    albedo=:dynamic,
    densification=:bessi,
    fresh_snow_density=:constant,
    Ntot=5,
)

forcing = SnowpackForcing(
    dt_days=[1.0, 1.0, 1.0],
    air_temperature_c=[-15.0, -12.0, -10.0],
    snowfall_mm_day=[2.0, 0.0, 0.0],
    rainfall_mm_day=0.0,
    shortwave_down=[80.0, 150.0, 220.0],
)

simulation = Simulation(model;
    forcing=forcing,
    years=1,
    backend=:threads,
    write_netcdf=false,
)

result = run!(simulation)
state = simulation.now
```

For coupled runs, use the initialized stepper API:

```julia
integrator = init_integrator(simulation)
integrator.sim.forcing.air_temperature .= 265.15
sync_forcing!(integrator)
step!(integrator, 1.0, true)
result = finalize!(integrator)
```

`simulation.ref` remains the reference state; `simulation.now` is the evolving
state.

## File-Backed Forcing

`load_forcing_file` reads a NetCDF forcing file and returns a named tuple with a
`SnowpackGrid` and `SnowpackForcing`.

```julia
loaded = load_forcing_file("forcing.nc")

simulation = Simulation(BESSIModel(loaded.grid);
    forcing=loaded.forcing,
    years=1,
    netcdf_variables=:all,
    write_netcdf=true,
)

result = run!(simulation)
```

The generic loader accepts common `(time, y, x)` layouts and MAR-style
`(x, y, TIME)` fields, including singleton layer dimensions such as
`(x, y, ATMLAY, TIME)`.

Spatial selection can be applied while loading:

```julia
loaded = load_forcing_file("forcing.nc";
    mask_name="MSK",
    mask_threshold=50.0,
)
```

## Running Scripts

There is an example script running the Greenland ice sheet on 10km resolution forced by MAR climatology in `examples/scripts/`.

Edit the `CONFIG` and `MODEL_OPTIONS` blocks in
`examples/scripts/run_gris_forcing_file_case.jl`, then run:

```bash
julia --project=. examples/scripts/run_gris_forcing_file_case.jl
```


## Documentation

Read the [Chion.jl documentation](https://fesmc.github.io/Chion.jl/dev/).
The documentation source lives in `docs/src/`. Build it locally with:

```bash
julia --project=docs docs/make.jl
```

Useful return values:

- `result.status`: run termination status
- `result.history`: per-year BESSI summary records
- `result.netcdf_path`: NetCDF output path when enabled

## Outputs And Backends

- `backend=:threads` runs on CPU.
- `backend=:cpu` is accepted as an alias for `:threads`.
- `backend=:gpu` runs on a CUDA device when available.
- `write_netcdf=false` skips NetCDF entirely.
- NetCDF output requires spatial grid coordinates.
- `history_year_stride` controls how often BESSI year metrics are recorded.
- `netcdf_variables=:all` selects the fields exposed by the chosen model.
- `netcdf_variables=:monthly` writes monthly BESSI or PDD diagnostics.
