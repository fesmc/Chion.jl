```@meta
CurrentModule = Chion
```

# Advanced Reference

These helpers are available for advanced workflows that inspect state directly,
manage device storage, or call low-level stepping APIs.

## State Inspection

```@docs
SnowpackState
get_state
print_state
```

## Backends

```@docs
cuda_available
```

## Workspaces And Timing

```@docs
ColumnarStepWorkspace
StepTimingStats
add_timing!
timing_rows
print_timing_summary
```
