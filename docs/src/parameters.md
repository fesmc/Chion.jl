# Default Parameters

Default values of the multilayer model (`BESSIModel`) and the bulk models
(`PDDModel`, `ITMModel`). Scheme and layer settings are `BESSIModel`
keywords. The physical constants are forwarded to `SnowpackPhysicalConstants`
and can be overridden through the same model keywords, for example
`BESSIModel(grid; semix_sensible_exchange_factor=2.5, alpha_dry=0.87)`. The
current values of a model are in `model.c`.

The surface-energy-balance defaults were calibrated for the Greenland ice
sheet against MAR v3.14.3 (1980–1989).


## Model

| Parameter | Default | Description |
| --- | --- | --- |
| `albedo` | :dynamic | Surface-albedo scheme. |
| `seb_scheme` | :semix | Longwave formulation: `:semix` graybody with emissivity on incident and emitted longwave, `:bessi` legacy (incident fully absorbed). |
| `turbulent_flux_scheme` | :semix | Turbulent-flux scheme. |
| `longwave_scheme` | :cloud_proxy | Parameterisation of downwelling longwave radiation. |
| `densification` | :bessi | Low-density densification scheme. |
| `fresh_snow_density` | :constant | Fresh-snow density scheme. |
| `number_type` | Float64 | Floating-point precision. |

## Layers

| Parameter | Default | Description |
| --- | --- | --- |
| `Ntot` | 15 | Maximum number of snow/firn layers. |
| `mass_max` | 500 kg m⁻² | Mass above which a layer is split. |
| `mass_split` | 300 kg m⁻² | Mass of the lower part of a split layer. |
| `mass_min` | 100 kg m⁻² | Minimum surface-layer mass (only without near-surface thickness limits). |
| `near_surface_layer_max_thicknesses_m` | 0.02, 0.05, 0.10, 0.30 m | Target thicknesses of the uppermost four layers (`Inf` disables a level). |
| `ice_substrate_layers` | 5 | Number of thermal ice-substrate layers. |
| `ice_substrate_top_thickness_m` | 0.05 m | Thickness of the uppermost substrate layer (doubling with depth). |
| `density_init` | 300 kg m⁻³ | Initial snow density. |
| `temperature_init` | 273 K | Initial snow temperature. |
| `refreezing_correction` | 1.0 | Multiplier on the cold-content-limited refreezing capacity. |

## Diurnal

| Parameter | Default | Description |
| --- | --- | --- |
| `diurnal_shortwave_substeps` | true | Enable diurnal shortwave substeps. |
| `diurnal_shortwave_threshold` | 0 W m⁻² | Minimum excess of the noon peak over the daily mean. |
| `diurnal_shortwave_max_substeps` | 8 | Number of daily substeps. |
| `diurnal_shortwave_min_air_temperature_c` | −8 °C | Minimum air temperature for substepping. |
| `diurnal_temperature_cycle` | true | Enable diurnal air-temperature cycle. |
| `diurnal_temperature_amplitude_c` | 1 K | Base diurnal temperature amplitude. |
| `diurnal_temperature_amplitude_gradient_c_per_km` | 0 K km⁻¹ | Elevation gradient of the amplitude. |
| `diurnal_temperature_amplitude_reference_height_m` | 0 m | Reference elevation for the amplitude. |
| `diurnal_temperature_amplitude_max_c` | 1 K | Maximum diurnal temperature amplitude. |

## Material

| Parameter | Default | Description |
| --- | --- | --- |
| `rho_s`, `rho_i`, `rho_w` | 315, 917, 1000 kg m⁻³ | Fresh-snow, ice, and water densities. |
| `rho_s_a`, `rho_s_b`, `rho_s_c` | 109, 6, 26 | Coefficients of the parameterised fresh-snow density. |
| `ci`, `cw` | 2110, 4181 J kg⁻¹ K⁻¹ | Ice and water heat capacities. |
| `Lm`, `Lv` | ``3.34\times10^5``, ``2.501\times10^6`` J kg⁻¹ | Latent heats of melting and vaporisation. |
| `cp_air` | 1003 J kg⁻¹ K⁻¹ | Specific heat capacity of air. |
| `T0` | 273.15 K | Melting temperature. |
| `seconds_per_day` | 86400 s | Length of a model day. |

## Albedo

| Parameter | Default | Description |
| --- | --- | --- |
| `alpha_dry`, `alpha_wet`, `alpha_ice` | 0.81, 0.70, 0.40 | Dry-snow, wet-snow, and bare-ice albedos. |
| `max_lwc_albedo` | 0.1 | Liquid-water content at which the albedo relaxes fully to the wet-snow value. |
| `aging_cold_timescale_days` | 20 d | Cold-surface e-folding time (aging scheme). |
| `aging_melting_timescale_days` | 2 d | Melting-surface e-folding time (aging scheme). |

## Radiation

| Parameter | Default | Description |
| --- | --- | --- |
| `ϵ_snow`, `eps_ice` | 0.98, 0.98 | Snow and ice emissivities. |
| `lw_clear_sky_transmissivity` | 0.718 | Clear-sky transmissivity at sea level and 300 W m⁻² insolation. |
| `lw_clear_sky_transmissivity_per_km` | 0.0298 km⁻¹ | Elevation gradient of the clear-sky transmissivity. |
| `lw_clear_sky_transmissivity_insolation` | 0.0569 | Insolation dependence of the clear-sky transmissivity. |
| `lw_cloud_transmission_coefficient`, `lw_cloud_transmission_exponent` | 0.447, 2.75 | Cloud transmission ``1-a\,n^p e^{-bz}``. |
| `lw_cloud_transmission_decay_per_km` | 0.225 km⁻¹ | Elevation decay ``b`` of the cloud effect. |
| `lw_clear_sky_emissivity`, `lw_clear_sky_emissivity_temperature_slope` | 0.716, 0.00239 K⁻¹ | Clear-sky emissivity at ``T_0`` and its temperature dependence. |
| `lw_overcast_emissivity`, `lw_overcast_emissivity_temperature_slope` | 0.952, 0.00362 K⁻¹ | Overcast emissivity at ``T_0`` and its temperature dependence. |
| `lw_cloud_emissivity_exponent` | 3 | Exponent ``q`` of the cloud weighting. |
| `lw_night_cloud_fraction` | 0.85 | Cloud fraction in the polar night and for sub-daily forcing. |
| `lw_emissivity_max` | 1.0 | Upper bound of the effective emissivity. |
| `ϵ_air` | 0.80 | Constant atmospheric emissivity (graybody scheme). |
| `σ` | 5.670373e-8 W m⁻² K⁻⁴ | Stefan–Boltzmann constant. |

## Turbulence

| Parameter | Default | Description |
| --- | --- | --- |
| `semix_karman` | 0.4 | von Kármán constant. |
| `semix_surface_height` | 10 m | Reference height. |
| `semix_z0m_snow`, `semix_z0m_ice` | 0.001, 0.01 m | Momentum roughness lengths of snow and ice. |
| `semix_zm_to_zh` | 10 | Ratio of momentum to heat roughness length. |
| `semix_sensible_exchange_factor` | 3.0 | Sensible-heat exchange factor ``f_H``. |
| `semix_latent_exchange_factor` | 1.0 | Latent-heat exchange factor ``f_E``. |
| `semix_stable_coefficient` | 40 | Stability coefficient ``b`` for ``Ri_B\geq0``. |
| `D_sh` | 10 W m⁻² K⁻¹ | Sensible-heat transfer coefficient (BESSI scheme). |
| `latent_heat_flux_ratio` | 1.0 | Ratio ``r_{\mathrm{lh/sh}}`` (BESSI scheme). |

## PDD

| Parameter | Default | Description |
| --- | --- | --- |
| `ddf_snow`, `ddf_ice` | 3, 8 kg m⁻² K⁻¹ d⁻¹ | Degree-day factors for snow and ice. |
| `refreezing_fraction` | 0.6 | Maximum refreezing fraction of the remaining snow. |
| `temperature_sigma` | 5 K | Standard deviation of temperature (stochastic PDD). |

## PDD, ITM

| Parameter | Default | Description |
| --- | --- | --- |
| `H_snow_max` | 5000 kg m⁻² | Maximum snow reservoir. |

## ITM

| Parameter | Default | Description |
| --- | --- | --- |
| `trans_a`, `trans_b` | 0.46, 6e-5 m⁻¹ᐟ² | Atmospheric transmissivity ``\tau_a=a+b\sqrt{z_s}``. |
| `itm_c`, `itm_b`, `itm_lat0` | −45 W m⁻², −2 W m⁻² deg⁻¹, 65 °N | Melt offset and its latitude dependence. |
| `itm_t` | 10 W m⁻² K⁻¹ | Temperature sensitivity of melt. |
| `alb_snow_dry`, `alb_snow_wet`, `alb_ice` | 0.80, 0.65, 0.40 | Snow and ice albedos. |
| `alb_land`, `alb_forest`, `alb_ocean` | 0.20, 0.10, 0.10 | Background albedos of ice-free surfaces. |
| `melt_crit` | 0.5 kg m⁻² d⁻¹ | Melt rate above which snow is wet. |
| `H_snow_crit_desert`, `H_snow_crit_forest` | 10, 100 kg m⁻² | Snow depth for full snow albedo. |
| `Pmaxfrac` | 0.6 | Maximum refreezing fraction of rain and meltwater. |

## Antarctic Set

The Antarctic paper runs override the Greenland calibration with

| Parameter | Value |
| --- | --- |
| `lw_clear_sky_emissivity` | 0.7314 |
| `lw_clear_sky_emissivity_temperature_slope` | 0.003378 K⁻¹ |
| `lw_overcast_emissivity` | 0.9877 |
| `lw_overcast_emissivity_temperature_slope` | 0.004171 K⁻¹ |
| `lw_cloud_emissivity_exponent` | 2.0 |
| `lw_night_cloud_fraction` | 0.9 |
| `semix_sensible_exchange_factor` | 2.5 |
| `semix_stable_coefficient` | 20 |
| `semix_latent_exchange_factor` | 3.0 |
| `alpha_dry`, `alpha_wet` | 0.87, 0.82 |
