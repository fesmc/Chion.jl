# Chion.jl Documentation

```@meta
CurrentModule = Chion
```

Chion is an intermediate-complexity snowpack mass and energy balance model
with layered solid mass, liquid water, density, temperature, and surface
diagnostics. It also includes a bulk positive-degree-day model. The public
workflow is simulation-first: build a grid, choose `BESSIModel` or `PDDModel`,
provide forcing, construct a `Simulation`, and call `run!`.

This documentation treats the current Julia implementation as the source of
truth. The process pages below summarize the formulas and defaults that are
actually implemented in `src/`.
## Recommended Reading Order

1. [Model State And Step Flow](model_state.md)
2. [Process documentation](processes/albedo.md)
3. [Simulation API And Outputs](case_api.md)
4. [Reference Utilities](reference.md)
5. [Validation And Audit](validation.md)

## Quick Start

The smallest end-to-end workflow is:

```julia
using Chion

grid = SnowpackGrid(1)
model = BESSIModel(grid)
forcing = SnowpackForcing(
    dt_days=[1.0, 1.0],
    air_temperature_c=[-12.0, -10.0],
    snowfall_mm_day=[1.0, 0.0],
    rainfall_mm_day=[0.0, 0.0],
    shortwave_down=[120.0, 160.0],
)
simulation = Simulation(model; forcing=forcing, years=1, write_netcdf=false)
result = run!(simulation)
```

## Process Pages

- [Albedo](processes/albedo.md)
- [Accumulation And Melt](processes/accumulation_ablation.md)
- [Layer Structure And Basal Transfer](processes/layer_structure.md)
- [Densification](processes/densification.md)
- [Energy Balance](processes/energy.md)
- [Diurnal Shortwave Cycle](processes/diurnal_cycle.md)
- [Percolation](processes/percolation.md)
- [Refreezing](processes/refreezing.md)
