```@meta
CurrentModule = Chion
```

# Parameters

All physical constants and scheme choices of a `BESSIModel` are fields of one
immutable struct, [`Params`](@ref); the bulk models take [`PDDParams`](@ref)
and [`ITMParams`](@ref). The vertical discretisation of the columns is not a
parameter but part of the grid, [`Layering`](@ref) (see
[Models And Coupling](models.md)).

Every exchangeable process is a *scheme slot* of `Params` that takes an
instance of a scheme type, and the instance carries the coefficients of that
scheme. Constants shared by several schemes (the dry, wet and ice albedos, the
emissivities, the material constants) stay flat fields of `Params`:

```julia
params = Params(;
    NF = Float32,                                           # number type of the model
    albedo_scheme = AgingAlbedo(cold_timescale_days = 15),  # scheme with a coefficient
    turbulent_flux_scheme = SemixTurbulentFluxes(sensible_exchange_factor = 2.5),
    diurnal_cycle = DiurnalSubsteps(max_substeps = 4,
        temperature_cycle = NoDiurnalTemperatureCycle()),
    alpha_dry = 0.84,                                       # shared constant
)
model = BESSIModel(grid; forcing, params)
model.alpha_dry, model.albedo_scheme                        # forwarded to model.params
```

`Params(params; kwargs...)` replaces single fields of an existing `Params`.
Unknown keywords throw an error that lists the fields. Every number and the
floating-point coefficients of the schemes are converted to `NF`.

The defaults are the Greenland calibration against MAR v3.14.3 with daily
forcing (1980–1989 for the surface energy balance). They are a regression
contract: the paper runs depend on them.

## Named Setups

| Setup | Domain | Differences from `Params()` |
| --- | --- | --- |
| [`ParamsGrISv1`](@ref) | Greenland | none (the defaults) |
| [`ParamsAISv1`](@ref) | Antarctica | refitted `CloudProxyLongwave`, `SemixTurbulentFluxes` with weaker stable damping and a larger latent factor, brighter snow |

The Antarctic setup changes, relative to the Greenland calibration,

| Field | Value |
| --- | --- |
| `longwave_scheme.clear_sky_emissivity` | 0.7314 |
| `longwave_scheme.clear_sky_emissivity_temperature_slope` | 0.003378 K⁻¹ |
| `longwave_scheme.overcast_emissivity` | 0.9877 |
| `longwave_scheme.overcast_emissivity_temperature_slope` | 0.004171 K⁻¹ |
| `longwave_scheme.cloud_emissivity_exponent` | 2.0 |
| `longwave_scheme.night_cloud_fraction` | 0.9 |
| `turbulent_flux_scheme.sensible_exchange_factor` | 2.5 |
| `turbulent_flux_scheme.stable_coefficient` | 20 |
| `turbulent_flux_scheme.latent_exchange_factor` | 3.0 |
| `alpha_dry`, `alpha_wet` | 0.87, 0.82 |

Keywords of a named setup override single fields, e.g.
`ParamsAISv1(NF = Float32, alpha_dry = 0.86)`. In a configuration file the
setup is `setup = "GrISv1"` or `setup = "AISv1"` in `[model]` (see
[Running Simulations](simulations.md#Configuration-files)).

```@docs
Params
ParamsGrISv1
ParamsAISv1
```

## Schemes

| Slot of `Params` | Default | Alternatives | TOML names |
| --- | --- | --- | --- |
| `fresh_snow_density_scheme` | `ConstantFreshSnowDensity()` | `ParameterizedFreshSnowDensity()` | `constant`, `parameterized` |
| `albedo_scheme` | `DynamicAlbedo()` | `ConstantAlbedo()`, `AgingAlbedo()`, `SemixAlbedo()`, `PrescribedAlbedo()` | `dynamic`, `constant`, `aging`, `semix`, `prescribed` |
| `densification_scheme` | `BESSIDensification()` | `HTESSELDensification()` | `bessi`, `htessel` |
| `surface_energy_balance_scheme` | `SemixSurfaceEnergyBalance()` | `BESSISurfaceEnergyBalance()` | `semix`, `bessi` |
| `turbulent_flux_scheme` | `SemixTurbulentFluxes()` | `BESSITurbulentFluxes()` | `semix`, `bessi` |
| `longwave_scheme` | `CloudProxyLongwave()` | `GraybodyLongwave()` | `cloud_proxy`, `graybody` |
| `diurnal_cycle` | `DiurnalSubsteps()` | `NoDiurnalCycle()` | `substeps`, `none` |

The nested slots are `SemixAlbedo(snow_bands = …)` (`DangSnowAlbedo()`, the
default, or `WarrenWiscombeSnowAlbedo()`; TOML `dang`, `warren_wiscombe`) and
`DiurnalSubsteps(temperature_cycle = …)` (`DiurnalTemperatureCycle()`, the
default, or `NoDiurnalTemperatureCycle()`; TOML `sinusoidal`, `none`).
`using InteractiveUtils; subtypes(Chion.AbstractAlbedoScheme)` lists the
implementations of a slot.

### Fresh-snow density

See [Accumulation And Melt](processes/accumulation_ablation.md).

```@docs
AbstractFreshSnowDensityScheme
ConstantFreshSnowDensity
ParameterizedFreshSnowDensity
```

### Albedo

See [Albedo](processes/albedo.md). The albedos `alpha_dry`, `alpha_wet`,
`alpha_ice` and `max_lwc_albedo` are fields of `Params`.

```@docs
AbstractAlbedoScheme
DynamicAlbedo
ConstantAlbedo
AgingAlbedo
SemixAlbedo
AbstractSemixSnowAlbedo
DangSnowAlbedo
WarrenWiscombeSnowAlbedo
PrescribedAlbedo
```

### Densification

See [Densification](processes/densification.md).

```@docs
AbstractDensificationScheme
BESSIDensification
HTESSELDensification
```

### Surface energy balance, turbulent fluxes and longwave

See [Energy Balance](processes/energy.md). The emissivities `ϵ_snow`,
`eps_ice` and `ϵ_air` are fields of `Params`.

```@docs
AbstractSurfaceEnergyBalanceScheme
SemixSurfaceEnergyBalance
BESSISurfaceEnergyBalance
AbstractTurbulentFluxScheme
SemixTurbulentFluxes
BESSITurbulentFluxes
AbstractLongwaveScheme
CloudProxyLongwave
GraybodyLongwave
```

### Diurnal cycle

See [Diurnal Shortwave Cycle](processes/diurnal_cycle.md).

```@docs
AbstractDiurnalCycle
DiurnalSubsteps
NoDiurnalCycle
AbstractDiurnalTemperatureCycle
DiurnalTemperatureCycle
NoDiurnalTemperatureCycle
```

## Bulk Models

`PDDModel` takes [`PDDParams`](@ref), whose `method` slot selects the positive
degree days: `SimplePDD()` (the default) from the positive part of the mean
temperature, `PISMPDD()` from the Calov–Greve expectation integral at every
step. `ITMModel` takes [`ITMParams`](@ref), the coefficients of the Fortran
insolation–temperature–melt scheme.

```@docs
PDDParams
AbstractPDDMethod
SimplePDD
PISMPDD
ITMParams
```
