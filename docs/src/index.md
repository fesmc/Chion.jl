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
go_refreezing!
apply_accumulation!
apply_melt!
split_surface_layer!
merge_surface_layer!
merge_bottom_layer!
calc_density_gradient_HL80
calc_density_gradient_powerlaw_ref
get_state
print_state
```

## Undocumented Functions

- `step_density`
- `reset_column_at_index!`
- `_snow_thermal_conductivity`
- `_remove_surface_layer!`
