# Chion.jl

Chion is a fast, intermediate-complexity snowpack and firn model for ice sheets and seasonal snow,
written in Julia for CPUs and GPUs. It computes the surface mass and energy balance from
atmospheric forcing in independent vertical columns:

- **`BESSIModel`**, the multilayer model (after BESSI; Born et al., 2019): surface energy balance
  with bulk-aerodynamic turbulent fluxes and a parameterised longwave, heat conduction with a Robin
  surface boundary, densification, percolation and refreezing, a refined near-surface grid, a
  thermal ice substrate and diurnal substeps for daily forcing.
- **`PDDModel`**, a positive-degree-day model, and **`ITMModel`**, the insolation–temperature–melt
  model.

The defaults are calibrated against MAR v3.14.3 over Greenland with daily forcing
(`ParamsGrISv1()`, = `Params()`); `ParamsAISv1()` is the Antarctic setup.

This model is still work in progress: not everything is validated yet, and the API may change.

## Installation

```julia
using Pkg
Pkg.activate("path/to/Chion.jl")
Pkg.instantiate()
using Chion
```

NetCDF input and output use `NCDatasets.jl`. The GPU backend needs `CUDA.jl`, a weak dependency:
add it to your environment and run `using CUDA` before building a model with
`backend = CUDABackend()`.

## Quick start

One column with two days of hand-written forcing (temperatures in °C, precipitation in mm w.e. per
day):

```julia
using Chion

forcing = SnowpackForcing(
    dt_days = [1.0, 1.0],
    air_temperature_c = [-12.0, -10.0],
    snowfall_mm_day = [1.0, 0.0],
    rainfall_mm_day = 0.0,
    shortwave_down = [120.0, 160.0],
    latitude_deg = 72.0,            # used by the diurnal substeps and the longwave
)
model = BESSIModel(Grid(1); forcing)
sim = Simulation(model)
run!(sim)
model.temperature                   # layer temperatures (K); model.mass, model.density, …
```

## A gridded run

forcing → `Grid` → model (parameters, boundary conditions, backend) → `Simulation` (years, spin-up, outputs).

```julia
using Chion
using CUDA                                         # only for the GPU

forcing = Forcing("/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc";   # one file per year
    surface_fluxes = :parameterized,               # Chion computes the surface energy balance
    cache = "output/forcing_cache")                # decoded years: 0.1 s instead of ~3 s per MAR year
grid = Grid(forcing; layering = Layering(Ntot = 15))
model = BESSIModel(grid; forcing, params = Params(), backend = CUDABackend())   # CPU() by default

sim = Simulation(model;
    years = 1940:2026,
    spinup = CycledSpinup(1940:1969; cycles = 5),
    output = OutputConfig(dir = "output/greenland", outputs = [   # run directory
        Output("monthly.nc"; variables = :monthly, frequency = Monthly()),
        Output("budget_yearly.csv"; variables = :budget, frequency = Yearly(), aggregate = AreaSum()),
        Output("firn.nc"; variables = :firn, frequency = Every(Month(6)), kind = Snapshot()),
    ]))
run!(sim)                  # run!(sim; years = 1) runs one year and can be called again
```

The run directory holds the outputs, `log.txt`, the run metadata `chion_run.toml` and JLD2
restarts (`OutputConfig(restart_every = 10)`). The Greenland and Antarctic runs of the Chion paper
are named setups: `SimulationGreenlandMARv1(path; dir)` and `SimulationAntarcticaMARv1(…)`.

- **Forcing**: `Forcing(path; format, surface_fluxes, …)` reads MAR (`format = :mar`), single-site
  ALMA (`:alma`, e.g. ESM-SnowMIP) or any NetCDF file with a `variables` name mapping.
  `Forcing(snowpack_forcing)` wraps forcing in memory, `Forcing(year -> …, years; grid)` takes a
  function, and `transform` changes every year (e.g. a warming).
- **Spin-up**: `CycledSpinup(years; cycles)`, `ClimatologySpinup(forcing; years)` or `NoSpinup()`.
- **Outputs**: `Output(path; variables, frequency, aggregate, kind, period, regions, depths)`, as
  `.nc` or `.csv`. Frequencies are `EveryStep()`, `Daily()`, `Monthly()`, `Yearly()` or
  `Every(period)`. `AreaSum()` gives ice-sheet or regional totals (Gt, km²). Fluxes are period
  totals; state variables are period means (`Mean()`) or end-of-period snapshots (`Snapshot()`).
  `available_output_variables(model)` lists the variables.

## Parameterisations

Every exchangeable parameterisation is a type that carries its own coefficients, passed in a slot
of `Params`; the model code selects it by dispatch:

```julia
params = Params(
    albedo_scheme = AgingAlbedo(cold_timescale_days = 15),
    densification_scheme = HTESSELDensification(),
    turbulent_flux_scheme = SemixTurbulentFluxes(sensible_exchange_factor = 2.5),
    diurnal_cycle = NoDiurnalCycle(),
)
```

The surface and bottom boundary conditions of the heat solve are types too: `boundary =
BoundaryConditions(bottom = GeothermalFlux(0.05))`. Their values can be one per column and updated
during a run with `set_boundary!`. `subtypes(Chion.AbstractAlbedoScheme)` lists the
implementations of a scheme; the [parameter page](https://fesmc.github.io/Chion.jl/dev/parameters/)
lists every default.

## Configuration files

The same simulation as a TOML file: `[forcing]`, `[grid]` (= `Layering`), `[model]` (= `Params`
fields, schemes by name or as tables with `type`), `[spinup]`, `[[output]]`.

```bash
julia --project -m Chion my_run.toml --backend=gpu
julia --project -m Chion my_run.toml --set model.alpha_dry=0.84   # override any entry
```

In Julia: `run_config("my_run.toml"; backend, set)`.

## Coupling and stepping

- `time_step!(sim)` advances one forcing step; `yearly_step!(sim)` advances one year, the usual
  exchange interval with an ice-sheet model.
- `set_input!` and `set_active_mask!` change the forcing and the active columns between steps, and
  `set_boundary!` changes the boundary values (e.g. the ice temperature below the column).
- `integrate!(model, steps)` steps a model without a `Simulation`, as used for automatic
  differentiation.
- `to_backend(model, backend)` moves a model between devices; `cpu_state(model)` returns a host
  copy of its state.

## Backends and differentiability

- **CPU**: `CPU()` is the default, threaded over columns (`julia --threads N`).
- **GPU**: `CUDABackend()`, with one GPU thread per column; load CUDA.jl first.
- **AD**: Enzyme forward mode works through the column step on the CPU and the GPU. ForwardDiff
  works through a whole model built at `NF = Dual`. A Reactant backend (`ReactantBackend()`) is a
  prototype.
- **Number type**: `Grid(…; NF = Float32)` with `Params(; NF = Float32)` runs in single precision
  (about 2× faster on a GPU).

## ITM forcing

`ITMModel` is the Fortran-compatible insolation–temperature–melt bulk model. It requires
`latitude_deg`, `surface_height`, `ice_thickness` (m) and `annual_pdd` (K day) in the forcing.

- **Insolation input**: ITM uses `q_sw_net` when it is supplied, otherwise `shortwave_down`.
- **Top-of-atmosphere radiation**: this insolation must be top-of-atmosphere radiation, which ITM
  attenuates with its elevation-dependent transmissivity. Surface shortwave such as MAR `SWD` would
  be attenuated twice and strongly underestimate melt.
- **Units**: its budgets keep the Fortran model's rates in mm w.e. per day and cumulative
  diagnostics in mm w.e.

## Examples and the paper

The `examples/` folder is git-ignored; its configurations point to forcing data on the PIK cluster.

| Folder | Contents |
| --- | --- |
| `examples/quickstart.jl`, `examples/tutorial.jl` | the gridded Greenland run above, and a guided tour of the API |
| `examples/configs/` | complete simulations as TOML files |
| `examples/research/` | native 3-hourly forcing, calibration and benchmark scripts |

The runs and figures of the Chion model description paper are in a separate repository,
Chion-paper (a checkout next to this one). It has one short Julia file per model run and per
figure, Chion as a path dependency and a committed manifest.

## Documentation

Read the [Chion.jl documentation](https://fesmc.github.io/Chion.jl/dev/). Build it locally with:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```
