```@meta
CurrentModule = Chion
```

# Reference Utilities

These helpers are useful when inspecting state, preparing execution scratch
space, moving between CPU/GPU storage, or summarizing a domain after a run.

## State Inspection

```@docs
AbstractSnowpackDomain
get_state
print_state
compute_auxiliary!
variables
column_count
```

## Domain Summaries And Backends

```@docs
summarize_domain_state
summarize_domain_state!
cpu_domain
gpu_domain
cuda_available
kernelabstractions_available
```

## Workspaces And Timing

```@docs
StepWorkspace
threaded_workspaces
ColumnarStepWorkspace
StepTimingStats
add_timing!
timing_rows
print_timing_summary
```
