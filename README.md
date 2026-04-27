# Chion.jl

Chion is a fast, intermediate-complexity snowpack mass and energy balance model.
It supports single-column and gridded snowpack runs, threaded CPU execution,
single-GPU execution through `CUDA.jl`, and optional summary, CSV, and NetCDF
outputs.

The current public workflow is:

1. Build a `SnowpackGrid`.
2. Build a `BESSIModel`.
3. Build a `SnowpackForcing`, or load one with `load_forcing_file`.
4. Build a `Simulation`.
5. Execute it with `run!`.

This model is still work in progress. Not everything is validated yet, and the
API may still change.

## Features

- Snowpack mass and energy balance core with layer-based state
- Simulation-first public API for manual and file-backed forcing workflows
- Threaded CPU and single-GPU execution backends
- Optional summary, history CSV, and NetCDF outputs
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

grid = SnowpackGrid(CPU(), 1)

model = BESSIModel(grid;
    albedo=DynamicAlbedo(),
    densification=BESSIDensification(),
    fresh_snow_density=ConstantFreshSnowDensity(),
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
    cycles=1,
    backend=:threads,
    save=:none,
    write_outputs=false,
)

result = run!(simulation)
state = SnowpackState(simulation.model.domain)
```

## File-Backed Forcing

`load_forcing_file` reads a NetCDF forcing file and returns a named tuple with a
`SnowpackGrid` and `SnowpackForcing`.

```julia
loaded = load_forcing_file("forcing.nc")

simulation = Simulation(BESSIModel(loaded.grid);
    forcing=loaded.forcing,
    cycles=1,
    save=:final,
)

result = run!(simulation)
```

The generic loader accepts common `(time, y, x)` layouts and MAR-style
`(x, y, TIME)` fields, including singleton layer dimensions such as
`(x, y, ATMLAY, TIME)`.

## Running Scripts

The script entry points in `examples/scripts/` are the fastest way to run the
packaged workflows.

```bash
julia --project=. examples/scripts/run_synthetic_simulation.jl \
  --backend=threads \
  --cycles=3 \
  --nx=2 \
  --ny=2 \
  --no-nc
```

```bash
julia --project=. examples/scripts/run_gris_forcing_file_case.jl \
  --forcing-file=/path/to/forcing.nc \
  --backend=gpu \
  --cycles=2 \
  --no-output \
  --no-nc
```

`examples/scripts/run_gris_equilibrium.jl` is still available as a
performance-oriented Greenland script while the public examples move to the
simulation API.

## Documentation

Detailed documentation lives in `docs/src/`. Build it locally with:

```bash
julia --project=docs docs/make.jl
```

Useful return values:

- `result.status`: run termination status
- `result.history`: per-cycle summary records
- `result.summary_path`, `result.history_csv_path`, `result.netcdf_path`: output paths when enabled

## Outputs And Backends

- `backend=:threads` runs on CPU with threading.
- `backend=:cpu` runs the same CPU path without threaded column stepping.
- `backend=:gpu` runs on a CUDA device when available.
- `save=:none` skips NetCDF entirely.
- NetCDF output requires spatial grid coordinates.
- `write_outputs=false` skips summary and history CSV files.
- `history_stride` controls how often cycle metrics are recorded.

## Tests

Run the test suite with:

```bash
julia --project=. test/runtests.jl
```

The tests cover the public simulation workflow, file-backed forcing, NetCDF
output behavior, and cleanup checks for removed legacy API names.
