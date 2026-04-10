```@meta
CurrentModule = Chion
```

# Case API And Outputs

The case API separates input preparation from execution:

1. create or load reusable inputs as a [`CaseDefinition`](@ref)
2. turn those inputs into a runnable [`SnowpackCase`](@ref)
3. execute the case and inspect the returned [`RunResult`](@ref)

## Input Units And Normalization

- `prescribed_case(...)` accepts air temperature in Celsius and snowfall /
  rainfall in `mmWE day^-1`
- `ForcingData(...)` accepts either those user-facing units or model-native
  rates in `kg m^-2 s^-1`
- `SnowpackStepForcing` and the runtime kernels always operate in native units

## Outputs

`RunResult` stores the final domain, cycle history, run status, timing
diagnostics, and the output paths that were produced when file writing was
enabled. NetCDF variable selection is controlled with
[`CASE_NETCDF_VARIABLE_GROUPS`](@ref) and [`CASE_NETCDF_VARIABLES`](@ref).

## API Reference

```@docs
AbstractCaseSource
SyntheticCaseSource
MARCaseSource
CaseDefinition
SnowpackCase
physics
load_case
prescribed_case
synthetic_case
mar_case
build_case
run_case
ForcingData
GridLayout
SnowpackStateFields
RunConfig
RunResult
TimingStats
CASE_NETCDF_VARIABLE_GROUPS
CASE_NETCDF_VARIABLES
```
