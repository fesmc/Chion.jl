# Fresh-Snow Density Validation Note

This page records the fresh-snow-density formula used by the current codebase
and the analytical spot checks that were used when documenting it. These notes
are retained as validation material; they are not currently backed by a file in
`test/`.

## Parameterization

Fresh-snow density is

```math
\rho_{\mathrm{fresh}} =
\operatorname{clamp}\!\left(
a + b\,(T_{\mathrm{air}} - T_0) + c\,\sqrt{V},
50,
\rho_i
\right),
```

with implementation defaults

- ``a = 109\,\mathrm{kg\,m^{-3}}``
- ``b = 6\,\mathrm{kg\,m^{-3}\,K^{-1}}``
- ``c = 26``

where `V` is the wind speed and the final result is clamped to
``[50, \rho_i]``.

## Analytical Checks

### Default wind-speed behavior

`SnowpackForcing` uses `wind_speed = 5.0 m s^-1` when no explicit wind speed is
supplied. The generic forcing-file loader uses the same default unless
`wind_default` is changed.

### Formula check

For

```math
T_{\mathrm{air}} = T_0 + 2\,\mathrm{K},
\qquad
V = 9\,\mathrm{m\,s^{-1}},
```

the expected fresh-snow density is

```math
\rho_{\mathrm{fresh}} = 109 + 6 \cdot 2 + 26 \cdot \sqrt{9} = 199\,\mathrm{kg\,m^{-3}}.
```

This is the value the current `_fresh_snow_density` implementation produces.

## Diagnostic Figure

![Fresh-snow density versus wind speed](../assets/fresh_snow_density_vs_wind.png)
