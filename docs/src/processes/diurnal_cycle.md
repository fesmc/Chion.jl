```@meta
CurrentModule = Chion
```

# Diurnal Shortwave Cycle

Chion can subdivide a daily BESSI forcing step to resolve the daytime
shortwave peak. The scheme reconstructs an idealized solar cycle whose
full-day mean equals the supplied shortwave forcing. An optional sinusoidal
air-temperature cycle can be applied over the same substeps.

The feature is disabled by default and applies only to `BESSIModel`.

## Configuration

Enable it when constructing the model:

```julia
model = BESSIModel(grid;
    diurnal_shortwave_substeps=true,
    diurnal_shortwave_threshold=0.0,
    diurnal_shortwave_max_substeps=3,
    diurnal_shortwave_min_air_temperature_c=-8.0,
    diurnal_temperature_cycle=false,
    diurnal_temperature_amplitude_c=5.0,
)
```

`diurnal_shortwave=true` is an alias that also enables
`diurnal_shortwave_substeps`.

The options mean:

| Option | Default | Meaning |
| --- | ---: | --- |
| `diurnal_shortwave_substeps` | `false` | enable adaptive subdivision |
| `diurnal_shortwave_threshold` | `0 W m^-2` | required peak-minus-mean shortwave excess |
| `diurnal_shortwave_max_substeps` | `3` | number of substeps when subdivision activates; valid range 1–24 |
| `diurnal_shortwave_min_air_temperature_c` | `-8 degC` | subdivision requires daily mean air temperature above this value |
| `diurnal_temperature_cycle` | `false` | reconstruct a diurnal air-temperature cycle |
| `diurnal_temperature_amplitude_c` | `5 degC` | half-amplitude of that temperature cycle |

Finite `latitude_deg` forcing is required when diurnal shortwave substeps are
enabled. Solar longitude is derived from each forcing timestamp.

## Solar Geometry

For latitude ``\phi`` and solar longitude ``\lambda``, the implementation
approximates solar declination as

```math
\delta =
\sin^{-1}\left(\sin\epsilon\,\sin\lambda\right),
\qquad
\epsilon = 23.439291^\circ.
```

The sunset hour angle is

```math
h_0 =
\cos^{-1}\left(-\tan\phi\,\tan\delta\right),
```

with explicit polar-night and polar-day limits of ``0`` and ``\pi``.

For hour angle ``h``, the unscaled solar shape is

```math
\mu(h) =
\sin\phi\sin\delta + \cos\phi\cos\delta\cos h.
```

Only the daylight interval ``[-h_0,h_0]`` contributes. Its integral is

```math
I_{\mathrm{day}} =
2\left[
h_0\sin\phi\sin\delta +
\cos\phi\cos\delta\sin h_0
\right].
```

## Energy-Conserving Shortwave Reconstruction

Let ``\overline{Q}_{sw}`` be the supplied daily mean shortwave flux. Chion
scales the solar shape by

```math
S = \overline{Q}_{sw}\frac{2\pi}{I_{\mathrm{day}}}.
```

For a substep spanning hour angles ``[h_a,h_b]``, define the daylight-clipped
bounds

```math
d_a = \max(h_a,-h_0),
\qquad
d_b = \min(h_b,h_0).
```

When ``d_b>d_a``, the interval-mean shortwave is

```math
\overline{Q}_{sw,[a,b]} =
\frac{S}{h_b-h_a}
\left[
(d_b-d_a)\sin\phi\sin\delta +
\cos\phi\cos\delta(\sin d_b-\sin d_a)
\right].
```

Nighttime intervals receive zero shortwave. Integrating all intervals over
``[-\pi,\pi]`` recovers the original daily mean, so subdivision does not alter
the prescribed daily shortwave energy.

If net shortwave `q_sw_net` is prescribed, the same reconstruction is applied
to it. Otherwise it is applied to `shortwave_down`, and albedo is handled by
the normal surface energy calculation.

## Activation Logic

For each column and forcing step, Chion first estimates the noon peak

```math
Q_{\mathrm{peak}} =
S\left(\sin\phi\sin\delta+\cos\phi\cos\delta\right).
```

The step is subdivided only when all of the following are true:

1. diurnal substeps are enabled
2. `diurnal_shortwave_max_substeps > 1`
3. the forcing duration is between `0.75` and `1.25` days
4. the selected daily mean shortwave is positive
5. daily mean air temperature is strictly above the configured minimum
6. air temperature, latitude, and solar longitude are finite
7. `Q_peak - Q_mean` is strictly greater than
   `diurnal_shortwave_threshold`
8. the column is not in polar night

If every condition passes, the code uses exactly
`diurnal_shortwave_max_substeps`; otherwise it uses one normal forcing step.
It does not select an intermediate number of substeps.

## Optional Temperature Cycle

When `diurnal_temperature_cycle=true`, the reconstructed temperature follows
a cosine cycle centered on the supplied daily mean:

```math
T(h) = \overline{T} + A_T\cos h.
```

The value passed to a substep is its interval average,

```math
\overline{T}_{[a,b]} =
\overline{T} +
A_T\frac{\sin h_b-\sin h_a}{h_b-h_a}.
```

This places the warmest point at solar noon and preserves the full-day mean
temperature. When the option is disabled, every substep uses the original air
temperature.

## Substep Forcing

The full day is divided uniformly in hour angle from ``-\pi`` to ``\pi``.
Each substep receives:

- a timestep scaled by its fraction of the day
- reconstructed interval-mean shortwave
- optionally reconstructed interval-mean air temperature
- unchanged precipitation rates and other atmospheric forcing

Because precipitation remains a rate while the timestep is shortened, total
daily snowfall and rainfall are preserved. Each substep then executes the
normal BESSI column process sequence.
