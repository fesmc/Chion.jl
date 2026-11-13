```@meta
CurrentModule = Chion
```

# Running Experiments

An `Experiment` combines a model, its forcing, a spin-up and a list of outputs. `run!` spins the
model up, runs every forcing year and writes each output at its own frequency. The forcing is
read one year at a time, so long multi-file records fit in memory, also on a GPU node.

```julia
using Chion

forcing = Forcing("/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc"; surface_fluxes=:parameterized)
model = BESSIModel(forcing)

experiment = Experiment(model, forcing;
    years=1940:2026,
    spinup=CycledSpinup(1940:1969; cycles=5),
    outputs=[Output("smb_monthly.nc"; variables=:monthly, frequency=:monthly),
             Output("budget_yearly.csv"; variables=:budget, frequency=:yearly, aggregate=:area_sum)],
    backend=:gpu)                 # needs `using CUDA`; :threads for the CPU
result = run!(experiment)         # an ExperimentResult; result.state is the final state
```

`backend=:gpu` needs `using CUDA` first (configuration files run with
`run_config` load it themselves). `Experiment(...; state=...)` starts from a
given state instead of `initial_state(model)`.

## Forcing

| Constructor | Use |
|---|---|
| `Forcing("dir/file_{year}.nc"; …)` | one NetCDF file per year (missing years are skipped) |
| `Forcing("file.nc"; …)` | one file of any length; years are sliced from its time axis |
| `Forcing(snowpack_forcing; grid, area)` | forcing already in memory |
| `Forcing(year -> snowpack_forcing, years; grid, area)` | a function for preprocessing the loaders do not cover |

Keywords of the file constructors:

- `format`: `:mar` for MAR output (variable names, ice mask `MSK` ≥ 50 %, a constant wind of
  5 m s⁻¹, the cell area `AREA`), `:alma` for single-site forcing in the ALMA convention of
  land-surface intercomparisons such as ESM-SnowMIP (`Tair`, `Qair`, `Psurf`, `Wind`, `Rainf`,
  `Snowf`, `SWdown`, `LWdown` in SI units; gaps are filled with the last valid value and the
  specific humidity becomes relative humidity over ice), `:generic` for any other file, `:auto`
  to detect.
- `surface_fluxes` (default `:prescribed`): `:prescribed` uses the surface energy fluxes of the file (MAR: net shortwave
  from `SWD − SWU`, `LWD`, `SHF`, `LHF` and the albedo `AL2`; use `albedo=:prescribed` in the
  model), `:parameterized` leaves them to Chion, `:albedo` keeps only the albedo (Chion computes
  the energy fluxes with it; use `albedo=:prescribed`), `:as_loaded` keeps whatever the loader
  provides. With the default `:prescribed` and a default `BESSIModel` (dynamic albedo), the MAR
  fluxes are combined with Chion's own albedo; pass `surface_fluxes=:parameterized` for a fully
  parameterised energy balance.
- `variables`: names of a generic file, e.g. `(air_temperature="t2m", snowfall="sf",
  rainfall="rf", shortwave_down="ssrd", latitude="lat")`; other `load_forcing_file` keywords
  (units, defaults) are passed through.
- `files = Dict(year => path)` replaces single years; `mask = (name, threshold)` selects columns;
  `wind` is a constant or a variable name; `area` the column areas (km²) as a vector or variable
  name, and `area_file` a separate file that holds it; `years` restricts the years used;
  `transform = (forcing, year) -> forcing` modifies every year (e.g. a temperature change).

`forcing_grid(forcing)`, `forcing_years(forcing)`, `forcing_area(forcing)` and
`load_forcing_year(forcing, year)` give access to the data, and `slice_forcing(f, r)` restricts
a `SnowpackForcing` to the time indices `r`; `BESSIModel(forcing; kwargs...)`
builds a model on the forcing columns. `relabel_year(snowpack_forcing, year)` moves the time
stamps of one year's forcing to another year (the solar geometry stays), for running the weather
of base years as future years, e.g. for a delta-method projection.

## Spin-up

- `CycledSpinup(years; cycles=5)` runs the real forcing `years` in order, `cycles` times.
- `ClimatologySpinup(forcing; years=200)` repeats another `Forcing` `years` times: the whole
  file of a single-file climatology (e.g. a day-by-day mean year), or the first year of a
  per-year source. In a configuration file, `[spinup]` with `type = "climatology"` and
  `path` takes the other settings (format, surface fluxes, …) from `[forcing]`.
- `NoSpinup()` starts from the initial state, or from `Experiment(...; state=...)`.

`years` of the `Experiment` may repeat years, e.g. `repeat(1940:1969, 5)`, to record every cycle
of an equilibrium run.

## Outputs

```julia
Output(path; variables=:default, frequency=:monthly, aggregate=:columns, kind=:mean,
       period=nothing, regions=nothing, area=nothing, depths=0:0.5:20, snow_cover_threshold=5.0,
       columns=nothing, grid=nothing)
```

- `path`: `.nc` (NetCDF, CF time axis with period bounds) or `.csv` (area aggregates or one column).
- `frequency`: `:step`, `:daily`, `:monthly`, `:yearly`, or a period: `Hour(3)`, `Day(5)`,
  `Week(1)`, `Month(6)`, `Year(10)`. Periods follow the calendar (months from January, days from
  1970-01-01). Each record covers one period; its time is the start of the period and
  `time_bnds` gives the first and last step.
- `aggregate`: `:columns` (one value per column, on the x–y grid when the grid has
  coordinates), `:area_sum` (fluxes and stores in Gt, snow cover in km², firn air content in
  km³, intensive quantities as area-weighted means) or `:area_mean`. `regions` takes a
  NamedTuple of column masks for several regions.
- `kind`: state variables as period means (`:mean`) or end-of-period values (`:snapshot`).
- `period = (1980, 2017)` writes only those years.
- `columns` records only a range of columns (e.g. one member of an ensemble run side by side),
  with `grid` their `SnowpackGrid`.

NetCDF files from an `Experiment` use a `time` dimension with `time_bnds`.
`available_output_variables(model)` lists the variables a model can write. The table below is
for `BESSIModel`; `PDDModel` and `ITMModel` write the forcing fluxes and their own state and
budget fields (`snowpack_swe`, `smb_ice`, `runoff`, … and `H_snow`, `melt`, `refreezing`, …).

| Variable | Kind | Units | |
|---|---|---|---|
| `surface_smb` | flux | kg m⁻² | precipitation − runoff + vapour exchange |
| `smb_ice` | flux | kg m⁻² | net mass forcing to the ice sheet |
| `melt`, `runoff`, `refreezing`, `sublimation` | flux | kg m⁻² | |
| `vapor_mass` | flux | kg m⁻² | net vapour exchange (deposition − sublimation − evaporation) |
| `precipitation`, `snowfall`, `rainfall` | flux | kg m⁻² | from the forcing |
| `latent_heat_flux`, `sensible_heat_flux` | rate | W m⁻² | period mean |
| `shortwave_down`, `shortwave_net`, `longwave_down`, `longwave_net`, `rain_heat_flux` | rate | W m⁻² | period mean; net = absorbed, positive into the surface |
| `swe`, `snow_mass`, `liquid_water` | state | kg m⁻² | solid + liquid, solid, liquid |
| `thickness`, `bulk_density`, `firn_air_content` | state | m, kg m⁻³, m | |
| `snow_cover` | state | 1 | SWE ≥ `snow_cover_threshold`; its mean is the snow-cover frequency |
| `albedo`, `Tsrf`, `snow_age_days`, `mass_base`, `N` | state | | |
| `mass`, `mass_w`, `density`, `temperature` | snapshot | | per snow layer |
| `ice_temperature` | snapshot | K | per ice-substrate layer |
| `density_profile`, `temperature_profile` | snapshot | kg m⁻³, K | on the depth bins / at the depths `depths` |

Groups: `:budget` (SMB, precipitation, melt, runoff, refreezing, vapour exchange), `:monthly`
(the monthly mass-budget and surface energy balance maps), `:snow`, `:firn` and `:default`.

Flux totals of a fine frequency add up to those of a coarser one to rounding precision; the tests
check daily against monthly and yearly totals.

## Configuration files

The same experiment as TOML, run with `julia --project -m Chion file.toml` (Julia ≥ 1.11) or
`run_config("file.toml")`. `--backend=gpu` and `--output-dir=DIR` override the file, and
`--set section.key=value` overrides any entry (e.g. `--set model.alpha_dry=0.84`; in Julia
`run_config(path; set=["model.alpha_dry" => 0.84])`). `run_config(path; dry_run=true)` returns
the `Experiment` without running it. The file is copied next to the results, with the
overrides applied.

```toml
name = "greenland_parameterized"
backend = "gpu"
years = [1940, 2026]          # first and last year; `repeat = 5` runs them five times
output_dir = "output/greenland"

[forcing]
path = "/data/MAR/MARv3.14.3-10km-daily-ERA5-{year}.nc"
format = "mar"
surface_fluxes = "parameterized"
[forcing.files]
2026 = "/data/MAR/MARv3.14.3-10km-daily-ERA5-2026_repaired.nc"

[model]                       # BESSIModel keywords
type = "bessi"

[spinup]
type = "cycled"               # none, cycled or climatology
years = [1940, 1969]
cycles = 5

[[output]]
file = "monthly.nc"
variables = "monthly"
frequency = "monthly"         # step, daily, monthly, yearly, "5 days", "6 months", …

[[output]]
file = "budget_yearly.csv"
variables = "budget"
frequency = "yearly"
aggregate = "area_sum"
```

Further keys:

| Section | Keys |
|---|---|
| top level | `name`, `backend`, `output_dir`, `years` (`[first, last]`, `"1940:2026"` or one year), `repeat` |
| `[forcing]` | `path`, `format`, `surface_fluxes`, `files`, `mask = ["MSK", 50]`, `wind`, `area`, `area_file`, `variables = {air_temperature = "t2m", …}`, `options = {…}` (other `Forcing` keywords) |
| `[model]` | `type = "bessi"`, `"pdd"` or `"itm"`; any model keyword, e.g. `albedo = "aging"`, `number_type = "Float32"`, `near_surface_layer_max_thicknesses_m = [0.02, 0.05, 0.1, 0.3]` |
| `[spinup]` | `type`, `years`, `cycles`; for `climatology` a `path` and an optional `[spinup.forcing]` table that overrides `[forcing]` |
| `[[output]]` | `file`, `variables`, `frequency`, `aggregate`, `kind`, `period = [1980, 2017]`, `depths = {start = 0, stop = 20, step = 0.5}`, `snow_cover_threshold` |

## Reference

```@docs
Forcing
Output
Experiment
ExperimentResult
available_output_variables
forcing_grid
forcing_years
forcing_area
load_forcing_year
slice_forcing
CycledSpinup
ClimatologySpinup
NoSpinup
relabel_year
run_config
```
