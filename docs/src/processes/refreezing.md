```@meta
CurrentModule = Chion
```

# Refreezing

`go_refreezing!` converts retained liquid water into solid mass using the cold
content of subfreezing layers.

## Model Formulation

For each active layer with solid mass, liquid water, and `T < T0`, the code
computes

```math
Q_{\mathrm{cold}} = (T_0 - T)c_i m_s,
\qquad
Q_{\mathrm{lat}} = m_w L_m.
```

### Partial Refreezing

If ``Q_{\mathrm{cold}} < Q_{\mathrm{lat}}``, only part of the water freezes:

```math
\Delta m_{\mathrm{refreeze}} = \frac{Q_{\mathrm{cold}}}{L_m}.
```

The code then sets:

- `T -> T0`
- `m_s -> m_s + Δm_refreeze`
- `m_w -> m_w - Δm_refreeze`
- `density -> min(density * (m_s + Δm_refreeze) / m_s, rho_i)`

### Complete Refreezing

If ``Q_{\mathrm{cold}} \ge Q_{\mathrm{lat}}``, all liquid water freezes. Using
the pre-update values of ``m_s``, ``m_w``, and ``T``, the code sets

```math
T^{new} =
\frac{m_w L_m / c_i + m_w T_0 + T m_s}{m_w + m_s}.
```

It then updates:

- `m_s -> m_s + m_w`
- `m_w -> 0`
- `density -> min(density * (m_s + m_w) / m_s, rho_i)` using the pre-update
  liquid mass

## API

```@docs
go_refreezing!
```

```@docs; canonical=false
_go_refreezing!
```
