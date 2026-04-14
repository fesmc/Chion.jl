# Chion.jl


## This model is still work in progress! Not everything is validated yet and many things still change. 

Chion is a fast, intermediate-complexity snowpack mass and energy balance model.
It supports single-column and gridded runs, threaded CPU execution, single-GPU execution through `CUDA.jl`, and high-level case builders for synthetic and prescribed forcings.

The current public workflow is:

1. Choose physics with `physics(...)`
2. Create a runnable case with `prescribed_case(...)` or `synthetic_case(...)`
3. Execute with `run_case(...)`

## Features

- Snowpack mass and energy balance core with layer-based state
- High-level case API for synthetic and prescribed forcing inputs
- Threaded CPU and single-GPU execution backends
- Optional summary, history CSV, and NetCDF outputs
- Pluto notebooks and script entry points for demos and larger runs

## Repository Layout

- `src/`: package code and public API
- `examples/scripts/`: command-line runners for runs
- `examples/pluto/`: Pluto notebooks for CPU, GPU scaffolding
- `examples/shared/`: helper utilities used by the notebooks and scripts
- `test/`: case API smoke tests
- `docs/`: documentation sources and generated site artifacts

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

### Optional Runtime Dependencies

- GPU runs require a working CUDA environment that `CUDA.jl` can use.
- File-backed prescribed workflows require prepared HDF5/NetCDF forcing files.
- If `libnetcdf` is not on the default library path, set `NETCDF_LIB` explicitly.

The existing launch scripts assume modules such as `julia`, `hdf5`, `netcdf-c`, and optionally `cuda`.

## Quick Start

This is the smallest supported end-to-end workflow using prescribed forcing arrays.
High-level forcing inputs use Celsius and `mmWE/day`.

```julia
using Pkg
Pkg.activate(".")
using Chion

case = Chion.prescribed_case(
    physics=Chion.physics(
        albedo=:dynamic,
        densification=:bessi,
        fresh_snow_density=:constant,
    ),
    ntot=5,
    nx=1,
    ny=1,
    dt_days=[1.0, 1.0, 1.0],
    air_temperature_c=[-15.0, -12.0, -10.0],
    snowfall_mm_day=[2.0, 0.0, 0.0],
    rainfall_mm_day=0.0,
    shortwave_down=[80.0, 150.0, 220.0],
    run=Chion.RunConfig(
        name="demo",
        backend=:threads,
        cycles=1,
        write_outputs=false,
        write_netcdf=false,
    ),
)

result = Chion.run_case(case)
state = Chion.get_state(result.domain, 1)
```

## Documentation

Detailed documentation lives in `docs/src/` and the rendered site groups the
material into:

- model state and stepping flow
- process documentation for albedo, accumulation/melt, layer structure, densification, energy, percolation, and refreezing
- the high-level case API and runtime outputs
- reference utilities and validation notes

Build the docs locally with:

```bash
julia --project=docs docs/make.jl
```

Useful return values:

- `result.status`: run termination status
- `result.history`: per-cycle summary records
- `result.domain`: final domain state
- `result.summary_path`, `result.history_csv_path`, `result.netcdf_path`: output paths when enabled

## Case Builders

### `prescribed_case(...)`

Build a runnable case directly from user-supplied forcing arrays, or from a
prepared external forcing file via `forcing_file=...`.
Use this for experiments, notebooks, and file-backed forcing pipelines.

### `synthetic_case(...)`

Generate a runnable case with built-in synthetic forcing for smoke tests and demos.
Supports `variant=:single_column` and `variant=:multi_column`, with `nx` and `ny` controlling the grid size.

## Running Scripts

The script entry points in `examples/scripts/` are the fastest way to run the packaged workflows.

### Synthetic Case

```bash
julia --project=. examples/scripts/run_synthetic_case.jl \
  --backend=threads \
  --cycles=3 \
  --nx=2 \
  --ny=2 \
  --no-nc
```

### Prepared Forcing File

```bash
julia --project=. examples/scripts/run_gris_forcing_file_case.jl \
  --forcing-file=/path/to/prepared_forcing.nc \
  --backend=gpu \
  --cycles=2 \
  --no-output \
  --no-nc
```

There is also a configuration-driven variant:

```bash
julia --project=. examples/scripts/run_gris_forcing_file_case_configured.jl
```

Its defaults can be overridden with environment variables such as `FORCING_PATH`, `BACKEND`, `CYCLES`, `WRITE_OUTPUTS`, `WRITE_NETCDF`, and `NETCDF_VARIABLES`.

## Pluto Notebooks

Interactive examples live in [`examples/pluto/README.md`](examples/pluto/README.md).
The current notebook set includes:

- `01_cpu_workflows.jl`: CPU workflows using the public case API
- `02_gpu_single_device.jl`: single-device GPU workflow
- `90_forcing_file_external_scaffold.jl`: optional external forcing-file scaffold

For cluster launches, `pluto.sh` starts a Pluto server with a dedicated runtime depot.

## Outputs And Backends

- `RunConfig(backend=:threads)` runs on CPU with threading
- `RunConfig(backend=:gpu)` runs on a CUDA device
- `backend=:cpu` is accepted and normalized to `:threads`
- `write_outputs=false` and `write_netcdf=false` are useful for smoke tests and timing runs
- `history_stride` controls how often cycle metrics are recorded
- `netcdf_variables` accepts `all`, `none`, group names, or explicit variable names

NetCDF output requires a grid layout, so it is available for gridded synthetic and prescribed cases.

## Tests

Run the test suite with:

```bash
julia --project=. test/runtests.jl
```

The tests cover the high-level case API, including direct prescribed, synthetic, and file-backed prescribed smoke cases.
