# Chion.jl Documentation

```@meta
CurrentModule = Chion
```

Chion is an intermediate-complexity snowpack mass and energy balance model
with layered solid mass, liquid water, density, temperature, and surface
diagnostics. The current implementation supports single-column and gridded
cases, threaded CPU execution, single-device GPU execution through `CUDA.jl`,
and a higher-level case API for prescribed and synthetic runs.

This documentation treats the current Julia implementation as the source of
truth. The process pages below summarize the formulas and defaults that are
actually implemented in `src/`.
## Recommended Reading Order

1. [Model State And Step Flow](model_state.md)
2. [Process documentation](processes/albedo.md)
3. [Case API And Outputs](case_api.md)
4. [Reference Utilities](reference.md)
5. [Validation And Audit](validation.md)

## Quick Start

The smallest end-to-end workflow is:

1. Choose physics with `physics(...)`
2. Create a runnable case with `prescribed_case(...)` or `synthetic_case(...)`
3. Execute with `run_case(...)`

## Process Pages

- [Albedo](processes/albedo.md)
- [Accumulation And Melt](processes/accumulation_ablation.md)
- [Layer Structure And Basal Transfer](processes/layer_structure.md)
- [Densification](processes/densification.md)
- [Energy Balance](processes/energy.md)
- [Percolation](processes/percolation.md)
- [Refreezing](processes/refreezing.md)
