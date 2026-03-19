# Fresh-Snow Density Test

This page documents the checks in `test/test_fresh_snow_density.jl`.

## Parameterization

Fresh-snow density is prescribed as

```math
\rho_{\mathrm{fresh}} =
\operatorname{clamp}\!\left(
a + b\,(T_{\mathrm{air}} - T_0) + c\,\sqrt{V},
50,
\rho_i
\right),
```

with default coefficients

- ``a = 109\,\mathrm{kg m^{-3}}``
- ``b = 6\,\mathrm{kg m^{-3} K^{-1}}``
- ``c = 26\,\mathrm{kg m^{-3} (m s^{-1})^{-1/2}}``

where `V` is the absolute wind speed, `air_temperature` is air temperature, `T0` is the freezing point, and `rho_i` is ice density.

## What The Test Checks

## Test 1: Default wind speed behavior

`step!` defaults to `wind_speed = 5.0` ``\mathrm{m s^{-1}}``.  
The test confirms that calling `step!` without a wind-speed argument produces the same density as calling it explicitly with `wind_speed=5.0`.

## Test 2: Formula match during accumulation

The direct formula check is performed through `apply_accumulation!`, not `step!`.  

For the case

```math
T_{\mathrm{air}} = T_0 + 2\,\mathrm{ K},
\quad
V = 9\,\mathrm{m s}^{-1},
```

the expected fresh-snow density is

```math
\rho_{\mathrm{fresh}} = 109 + 6 \cdot 2 + 26 \cdot \sqrt{9} = 199\,\mathrm{ kg m}^{-3}.
```

The test verifies that the newly created surface layer matches this value to machine precision.

## Test 3: Diagnostic plot generation

The test also generates a plot of fresh-snow density versus wind speed for several subfreezing air temperatures and checks that the output file is created successfully.

![Fresh-snow density versus wind speed](../assets/fresh_snow_density_vs_wind.png)

## Run command

```bash
julia --project=. test/test_fresh_snow_density.jl
```

## Results snapshot

Recorded on **2026-03-18**:

| Test block | Passed assertions |
| --- | ---: |
| Default wind speed behavior | 5/5 |
| Formula match during accumulation | 2/2 |
| Plot generation | 2/2 |
| **Total** | **9/9** |

Terminal summary:

```text
Activating project at `~/Documents/Chion.jl`
Test Summary:                       | Pass  Total  Time
Fresh-snow density parameterization |    9      9  1.0s
```
