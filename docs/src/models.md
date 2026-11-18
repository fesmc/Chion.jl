```@meta
CurrentModule = Chion
```

# Models And Coupling

For runs from forcing files with spin-up and outputs, see
[Running Simulations](simulations.md). This page describes the pieces
underneath a `Simulation`, for forcing held in memory, step-by-step coupling
to other models, GPUs and automatic differentiation:

1. build a [`Grid`](@ref): the columns, their optional x–y layout and the
   vertical [`Layering`](@ref)
2. build the forcing: a file-backed [`Forcing`](@ref) or an in-memory
   [`SnowpackForcing`](@ref) (see its docstring below for the optional fields)
3. build `BESSIModel`, `PDDModel` or `ITMModel` with its parameters
   ([`Params`](@ref), [`PDDParams`](@ref), [`ITMParams`](@ref))
4. wrap it in a `Simulation` and call `run!`, or step the model directly

```jldoctest models
julia> using Chion

julia> grid = Grid(2; layering = Layering(Ntot = 20));

julia> forcing = SnowpackForcing(dt_days = ones(3), ncol = 2, air_temperature_c = -5.0,
           snowfall_mm_day = 2.0, rainfall_mm_day = 0.0, shortwave_down = 100.0, latitude_deg = 70.0);

julia> model = BESSIModel(grid; forcing, params = Params(albedo_scheme = AgingAlbedo()));

julia> sim = Simulation(model; verbose = false);

julia> run!(sim);

julia> finished(sim), model.Ntot, size(model.mass)
(true, 20, (20, 2))
```

A model holds its `grid`, the prognostic `state`, a `cache` (diagnostics,
scratch arrays and the active columns), the `params` and the `forcing` (the
source and the loaded calendar year). Their fields are forwarded, so
`model.mass` and `model.Tsrf` read the state, `model.thickness` the cache,
`model.alpha_dry` and `model.albedo_scheme` the parameters, `model.ncol` and
`model.Ntot` the grid, and `model.year` the loaded forcing year. A model knows
nothing about time beyond the loaded forcing year; a `Simulation` advances it.

`BESSIModel` resolves the layered mass and energy balance. Its default
parameters are the setup calibrated against MAR over Greenland for daily
forcing (see [Parameters](parameters.md)). The default diurnal substeps
require `latitude_deg` in the forcing; use
`Params(diurnal_cycle = NoDiurnalCycle())` when latitude is unavailable or the
forcing is already sub-daily. `PDDModel` uses a bulk snow reservoir with the
degree-day factors `ddf_snow` and `ddf_ice`, `refreezing_fraction`,
`temperature_sigma`, `H_snow_max` and the `method` of [`PDDParams`](@ref). The
default `method = SimplePDD()` uses the positive mean temperature;
`method = PISMPDD()` applies the Calov–Greve expectation integral at every
timestep.

`ITMModel` ports the Fortran insolation–temperature–melt scheme. It requires
explicit `latitude_deg`, `surface_height`, `ice_thickness` (m), and
`annual_pdd` (K day) in `SnowpackForcing`. `q_sw_net`, when supplied, is its
insolation driver; otherwise ITM uses `shortwave_down`. As in the original
ITM, this insolation is top-of-atmosphere radiation that the model attenuates
with its elevation-dependent transmissivity; supplying surface shortwave (for
example MAR `SWD`) attenuates it twice and strongly underestimates melt. It
preserves the Fortran scheme's per-step `mmWE day^-1` fields (`smb`, `smbi`,
`melt`, `runoff`, `refreezing`) and cumulative `mmWE` fields.

## Number Type

`NF` sets the floating-point type of the grid, the parameters, the state and
the forcing fields that the kernels read. The grid and the parameters must
agree: `Grid(n; NF = Float32)` with `Params(; NF = Float32)`, or
`BESSIModel(forcing; NF = Float32)`, which builds both. `PDDModel` and
`ITMModel` take `PDDParams(; NF)` and `ITMParams(; NF)` in the same way.

## Backends

A model lives on a KernelAbstractions backend: `CPU()` (multithreaded, the
default) or, after `using CUDA`, `CUDABackend()`. CUDA.jl is a weak
dependency: `using CUDA` loads the GPU backend (the `ChionCUDAExt` package
extension).

```julia
using CUDA
model = BESSIModel(grid; forcing, backend = CUDABackend())
model = to_backend(model, CUDABackend())     # or move an existing model
state = cpu_state(model)                     # host copy of the state
set_state!(model, state)                     # copy a host state (e.g. a restart) back
```

## Coupling

A `Simulation` can be advanced one forcing step or one year at a time, with
the model's inputs changed in between:

```julia
sim = Simulation(model; years = 2000:2100)
time_step!(sim)                              # one forcing step (spin-up first, if any)
set_input!(sim, :surface_height, z_s)        # e.g. from an ice-sheet model, per column
set_active_mask!(sim, ice_mask)              # step only these columns
fields = yearly_step!(sim)                   # the next forcing year (must start at its first step)
fields.smb_ice_yearly                        # kg m⁻² per column; fields.mean_T_srf_K in K
run!(sim; years = 10)                        # ten more years
```

`set_input!` replaces a forcing field of the loaded year (after editing `model.forcing.data` in
place, `sync_forcing!(model)` moves it to the kernels); `set_active_mask!`
takes a column vector or a y–x matrix on the grid, and resets columns that
become inactive to the initial state unless
`reset_newly_inactive = false`. `yearly_step!` returns the annual coupling
fields `smb_ice_yearly` and `mean_T_srf_K`, also as y–x matrices
(`*_grid`) when the grid has spatial coordinates.

Without a `Simulation`, `integrate!(model, steps)` advances the model over the
forcing steps `steps` of the loaded year, without output and without changing
the year; `update_forcing!(model, year)` loads another year. `integrate!` is
the entry point for automatic differentiation.

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

A model takes a `SnowpackForcing` directly (it becomes a one-year
`Forcing`). `load_forcing_file` can select spatial columns directly with
`mask_name` and `mask_threshold`. Selection depends only on that mask
variable. Each field is read once straight into the selected columns;
dimensions are recognised by the names given as `x_name`/`y_name` (or
`x`/`y`) and a time dimension, and any other dimension (vertical level,
sector) is read at its first index. `read_forcing_columns` exposes the same
reader for additional variables.

## Host Forcing Helpers

A host whose climate supplies only monthly air temperature and precipitation
can derive the rest with `shortwave_down_from_transmissivity` (shortwave from
a per-domain direct transmissivity; the longwave cloud proxy then diagnoses
its own cloudiness), `wind_speed_from_geometry` (wind from elevation and
slope), `temperature_variability_sigma` and `red_noise_monthly` (daily
temperature variability about a smooth series). Chion does not call them
itself. `clear_sky_transmissivity` needs the parameters of a model, whose
`CloudProxyLongwave` holds the coefficients:

```julia
tau_cs = clear_sky_transmissivity(model.params, 2000.0, 300.0)   # height (m), TOA SW (W m⁻²)
```

## Reference

```@docs
Grid
Layering
SpatialLayout
BESSIModel
PDDModel
ITMModel
to_backend
cpu_state
set_state!
set_active_mask!
set_input!
sync_forcing!
update_forcing!
integrate!
yearly_step!
SnowpackForcing
update_air_pressure!
load_forcing_file
read_forcing_columns
clear_sky_transmissivity
shortwave_down_from_transmissivity
wind_speed_from_geometry
temperature_variability_sigma
red_noise_monthly
```
