# Chion.jl

Chion is a fast, intermediate-complexity snowpack mass and energy balance model.
It supports single-column and gridded snowpack runs, CPU execution,
single-GPU execution through `CUDA.jl`, and optional NetCDF output.

The current public workflow is:

1. Build a `SnowpackGrid`.
2. Build a `BESSIModel`, `PDDModel`, or `ITMModel`.
3. Build a `SnowpackForcing`, or load one with `load_forcing_file`.
4. Build a `Simulation`.
5. Execute it with `run!`.

This model is still work in progress. Not everything is validated yet, and the
API may still change.

## Features

- Layered snowpack mass and energy balance (BESSI) with a SEMIX surface
  energy balance, a Robin surface boundary, fine near-surface layers, and a
  thermal ice substrate below the snow
- Defaults calibrated against MAR v3.14.3 over Greenland for daily forcing:
  dynamic albedo, cloud-proxy longwave (no extra forcing needed), SEMIX
  sensible heat, and diurnal substeps
- Bulk positive-degree-day (`PDDModel`) and insolation-temperature-melt
  (`ITMModel`) models
- Simulation-first public API for manual and file-backed forcing workflows
- CPU and single-GPU execution backends
- Per-year BESSI history and optional NetCDF output, including monthly fields

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

## Running An Experiment

An experiment is a model, its forcing, a spin-up and the outputs. The forcing is read one year at
a time, so long multi-file records fit in memory, and each output is written at its own frequency.

```julia
using Chion

forcing = Forcing("/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc";   # one file per year
    surface_fluxes=:parameterized)          # Chion computes the surface energy balance itself
model = BESSIModel(forcing)                 # the calibrated defaults on the forcing columns

experiment = Experiment(model, forcing;
    years=1940:2026,
    spinup=CycledSpinup(1940:1969; cycles=5),
    outputs=[
        Output("smb_monthly.nc"; variables=:monthly, frequency=:monthly),               # maps
        Output("budget_yearly.csv"; variables=:budget, frequency=:yearly, aggregate=:area_sum),  # Gt
        Output("firn.nc"; variables=:firn, frequency=Month(6), kind=:snapshot),       # profiles
    ],
    backend=:gpu)                           # :threads for the CPU; :gpu needs `using CUDA`

result = run!(experiment)                   # result.state is the final state
```

The same experiment as a TOML file (see `examples/configs/`):

```bash
julia --project -m Chion examples/configs/greenland_parameterized.toml --backend=gpu
```

- `Forcing(path; format, surface_fluxes, …)` reads MAR output (`format=:mar`: variable names, ice
  mask, constant wind, cell areas), single-site ALMA forcing (`format=:alma`, e.g. ESM-SnowMIP) or
  any NetCDF file with a `variables` name mapping;
  `Forcing(snowpack_forcing)` wraps forcing in memory and `Forcing(year -> …, years; grid)` takes a
  function for preprocessing the loaders do not cover. A `transform` applies changes such as a
  temperature offset to every year.
- `CycledSpinup(years; cycles)`, `ClimatologySpinup(forcing; years)` or `NoSpinup()`.
- `Output(path; variables, frequency, aggregate, kind, period, regions)`: `.nc` or `.csv`;
  frequencies `:step`, `:daily`, `:monthly`, `:yearly` or any period (`Day(5)`, `Month(6)`);
  `aggregate=:area_sum` gives ice-sheet or regional totals (Gt, km²). Fluxes are period totals,
  state variables period means or end-of-period snapshots. `available_output_variables(model)`
  lists the variables; groups are `:budget`, `:monthly`, `:snow`, `:firn` and `:default`.

## Quick Start (Low-Level API)

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
    latitude_deg=72.0,  # needed by the default diurnal substeps
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
`(x, y, TIME)` fields. Dimensions are recognised by the names passed as
`x_name`/`y_name` (or `x`/`y`) and a time dimension; extra dimensions such as
vertical levels or sectors are read at their first index. The default diurnal
substeps need latitude (`LAT` in MAR files); set
`diurnal_shortwave_substeps=false` for forcing without latitude or with
sub-daily time steps.

Spatial selection can be applied while loading:

```julia
loaded = load_forcing_file("forcing.nc";
    mask_name="MSK",
    mask_threshold=50.0,
)
```

## Examples

| Folder | Contents |
| --- | --- |
| `examples/configs/` | complete experiments as TOML files, run with `julia --project -m Chion <file>.toml` |
| `examples/runs/` | short scripts for experiments that need code: Northern Hemisphere snow from ERA5, the ten ESM-SnowMIP sites |
| `examples/figures/` | the paper figures, and `compare_with_mar.jl <run dir>` for the yearly budget against MAR |
| `examples/research/` | native 3-hourly forcing, longwave fit and benchmarks behind the paper (low-level API, settings from `CHION_*` environment variables) |

All runs behind the paper figures are in `run_paper_experiments.sh` (`sbatch run_paper_experiments.sh
greenland|native|antarctica|snow`); it shows how each output is made from a config and `--set` overrides.

The configs cover Greenland and Antarctica (parameterized, prescribed, firn
profiles for SUMup). Each run writes its output and a copy of its configuration to `output_dir`.

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
- `backend=:gpu` runs on a CUDA device; load CUDA.jl first (`using CUDA`), which
  activates the `ChionCUDAExt` package extension.
- `write_netcdf=false` skips NetCDF entirely.
- NetCDF output requires spatial grid coordinates.
- `history_year_stride` controls how often BESSI year metrics are recorded.
- `netcdf_variables=:all` selects the fields exposed by the chosen model.
- `netcdf_variables=:monthly` writes monthly BESSI, PDD, or ITM diagnostics.

## ITM

`ITMModel` is the Fortran-compatible insolation-temperature-melt bulk model.
It requires explicit `latitude_deg`, `surface_height`, `ice_thickness` (m),
and `annual_pdd` (K day) forcing. ITM uses `q_sw_net` as its insolation input
when supplied; otherwise it uses `shortwave_down`. Its internal budgets retain
the Fortran model's `mmWE/day` rates and cumulative `mmWE` diagnostics. The
offset varies with latitude as `itm_c + itm_b*(|lat| - |itm_lat0|)`, the same
in both hemispheres (`itm_lat0` [deg] is the reference latitude, its absolute
value; `|itm_lat0| >= 90` disables this).
