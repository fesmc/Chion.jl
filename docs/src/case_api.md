```@meta
CurrentModule = Chion
```

# Case API And Outputs

The public case API is intentionally small:

1. create a runnable `SnowpackCase` with `prescribed_case` or `synthetic_case`
2. execute it with `run_case`
3. inspect the returned `RunResult`

## Input Units And Normalization

- `prescribed_case(...)` accepts air temperature in Celsius and snowfall /
  rainfall in `mmWE day^-1`
- `prescribed_case(...; forcing_file=...)` loads a prepared external forcing
  file and assumes any spatial masking has already been done outside Chion
- the internal forcing arrays and runtime kernels always operate in native
  units such as `K`, `kg m^-2 s^-1`, and `W m^-2`

## Outputs

`RunResult` stores the final domain, cycle history, run status, timing
diagnostics, and the output paths that were produced when file writing was
enabled. NetCDF variable selection is controlled with
`CASE_NETCDF_VARIABLE_GROUPS` and `CASE_NETCDF_VARIABLES`.

## Public Surface

- `physics(...)` selects the process parameterization bundle.
- `prescribed_case(...)` creates a runnable case from direct forcing arrays or a prepared forcing file.
- `synthetic_case(...)` creates a runnable case from built-in synthetic forcing.
- `run_case(case)` executes a case and returns a `RunResult`.
- `RunConfig(...)` controls backend, output writing, cycle count, and history stride.

## Reference

```@docs
prescribed_case
synthetic_case
run_case
```
