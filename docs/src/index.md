# Chion.jl Documentation

```@meta
CurrentModule = Chion.SnowpackModel
```
Chion is a fast, intermediate complexity snowpack model. 
It is inspired by the [BESSI](https://tc.copernicus.org/articles/13/1529/2019/) model. 
## Function Reference

```@docs
step!
go_densification!
go_energy_flux!
go_percolation!
continuous_bottom_deplete!
go_refreezing!
apply_accumulation!
apply_melt!
split_surface_layer!
merge_surface_layer!
merge_bottom_layer!
get_state
print_state
```

## Undocumented Functions

- `step_density`
- `reset_column_at_index!`
- `_snow_thermal_conductivity`
- `_remove_surface_layer!`
