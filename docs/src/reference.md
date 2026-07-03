```@meta
CurrentModule = Chion
```

# Advanced Reference

These helpers support advanced workflows that inspect state directly, manage
device storage, or inspect timing information. Some are intentionally
module-qualified rather than exported.

## State Inspection

```@docs
get_state
print_state
```

## Backends

```@docs
cuda_available
cpu_state
gpu_state
```

## Workspaces And Timing

```@docs
StepTimingStats
add_timing!
timing_rows
print_timing_summary
```
