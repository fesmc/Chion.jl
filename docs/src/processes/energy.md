```@meta
CurrentModule = Chion.SnowpackModel
```

# Energy Balance

`go_energy_flux!` solves layer temperatures with an implicit 1D conductive step plus surface flux forcing.

## Surface Flux Parameterization

The surface energy is given by the sum of net shortwave radiation ``Q_{sw}``, thermal longwave radiation ``Q_{lw}``, sensible heat exchange with the atmosphere ``Q_{sh}``, heat transport by precipitation ``Q_p`` and the exchange of latent heat due to refreezing or melting or rain and meltwater ``Q_{lh}``:

```math
c_i\, m_s \,\frac{\partial T_s}{\partial t}
= Q_{sw} + Q_{lw} + Q_{sh} + Q_p + Q_{lh}
```

with the heat capacity of ice ``c_i`` and the mass per area of the surface box ``m_{s,i}``.



### Shortwave Radiation 
The shortwave radiation is given by
```math
Q_{sw} = (1-\alpha)S_{boa}
```

where the albedo ``\alpha`` is temperature dependent. For ``T_s=0^\circ\mathrm{C}``, the albedo is given by ``\alpha=\alpha_{wet}`` and ``\alpha_{dry}`` otherwise. ``S_{boa}`` is the incoming solar radiation at the bottom of the atmosphere that is absorbed by the surface snow layer.

### Longwave radiation


```math
Q_{\mathrm{lw}} = \sigma\left(\epsilon_{\mathrm{air}} T_{\mathrm{air}}^{4} - \epsilon_{\mathrm{snow}} T_{s}^{4}\right)\,.
```
### Sensible Heat 
```math
Q_{\mathrm{sh}}= D_{\mathrm{sh}}(T_{\mathrm{air}}-T_{\mathrm{s}})
```

with 
```math
D_{\mathrm{sh}}= 1.29\cdot 10^{-2}\mathrm{K}^{-1}\cdot Apu.
```
Here ``D_{\mathrm{sh}}=10\,\mathrm{Wm^{-2}K^{-1}}``.

### Precipitation Heat 
```math
Q_{p,\mathrm{rain}} = P \,\rho_w \, c_w \,(T_{\mathrm{air}} - T_0),
\qquad
Q_{p,\mathrm{snow}} = P \,\rho_w \, c_i \,(T_{\mathrm{air}} - T_s)\,.
```

``Q_{p}`` can be rewritten as 

```math
Q_{p,\mathrm{snow}}(T_s) + K_{\mathrm{lh}} - H_{\mathrm{lh}}T_s
```

with 

```math
K_{\mathrm{lh}} = P \,\rho_w \, c_x
```

```math
H_{\mathrm{lh}} = P \,\rho_w \, c_x \, T_{\mathrm{air} }
```

with the heat capacity of snow or rain ``c_x``. 
In other words, ``H_{\mathrm{lh}}`` is the slope (derivative) of the precipitation heat flux w.r.t. surface temperature and ``K_{\mathrm{lh}}`` is the  ``T_s``-independent part.

In the case of rainfall (no snow), ``H_{\mathrm{lh}}=0`` and, therefore, rain contributes a temperature-independent heat source.

### Longwave Radiaton 




The emitted longwave radiation is nonlinear in the surface temperature ``\epsilon_{\mathrm{snow}} T_{s}^{4}`` making the system difficult to solve. 
To overcome this problem, the temperature at the surface at the time step ``n+1``, i.e., ``(T^{n+1}_s)`` is linearised around the temperature at the previous time step ``n`` (Taylor expansion): 

```math 
(T^{n+1}_s)^4 \approx 4(T^n)^3T^{n+1}-3(T^n)^4
```
Then the longwave radiation is linear in ``T^{n+1}_s`` and the surface energy can be split into a constant part and a linear part:

```math
Q(T^{n+1}_s) \approx Q_{\mathrm{const}}(T^n_s) - Q_{\mathrm{lin}}(T^n_s)\,T^{n+1}_s.
```

Collecting all the ``T^{n+1}_s`` (in)dependent terms gives

```math
Q_{\mathrm{const}} =
D_{sh}T_{2m}
+\sigma\left(\epsilon_{air}T_{2m}^4 + 3\epsilon_{snow}(T_s^n)^4\right)
+Q_{sw}
+K_{lh}
```

```math
Q_{\mathrm{lin}} =
D_{sh} + 4\sigma\epsilon_{snow}(T_s^n)^3 + H_{lh}.
```

##  Diffusion 

The diffusion equation is given by 

```math
c_i \rho_s \frac{\partial T_s}{\partial t}
= \frac{\partial}{\partial z}\!\left( K(\rho_s)\,\frac{\partial T_s}{\partial z} \right),
```

with the thermal conducitvity of snow (Yen, 1981):

```math
K(\rho_s) = K_i \left(\frac{\rho_s}{\rho_w}\right)^{1.88}.
```

Discretizing the equation gives for each interior layer ``i`` with thickness ``\Delta z_i``

```math
T_i^{n+1} - T_i^n = \frac{\Delta t}{\rho_i c \Delta z_i}\left(F^{n+1}_{i-1/2}  - F^{n+1}_{i+1/2} \right)
```

with the conductive heat flux across the interface between layers ``i`` and ``i+1`` (Fourier law):

```math
F^{n+1}_{i+1/2} = -K_{i+1/2} \frac{T^{n+1}_{i+1} - T^{n+1}_{i} }{\Delta z_{i+1/2} } = -G_{i+1/2}({T^{n+1}_{i+1} - T^{n+1}_{i} })
```

So we have a system of ``n`` (number of layers) equations in ``n`` unknowns ``T_i^{n+1}``.

The layer-center spacing is given by 

```math
\Delta z_{i+1/2} = \frac{\Delta z_i = \Delta z_{i+1}}{2}.
```

The interface conductivity ``K_{i+1/2}`` is given by a thickness-weighted arithmetic mean

```math
K_{i+1/2} = \frac{K_i \Delta z_i + K_{i+1}\Delta z_{i+1}}{\Delta z_{i} + \Delta z_{i+1}}
```

then 

```math
G_{i+1/2} = \frac{2(K_i \Delta z_i + K_{i+1}\Delta z_{i+1})}{(\Delta z_{i} + \Delta z_{i+1})^2}.
```

Now considering the contribution of the upper interface ``i+1/2`` to layer ``i``, we get 

```math
F^{n+1}_{i+1/2} = -G_{i+1/2}({T^{n+1}_{i+1} - T^{n+1}_{i} })
```

and putting this term into the layer equation gives

```math
T_i^{n+1} - T_i^n =  \frac{\Delta t}{\rho_i c \Delta z_i}G_{i+1/2}\left(\dots  - F^{n+1}_{i+1/2} \right) = \frac{\Delta t}{\rho_i c \Delta z_i} G_{i+1/2} (T_{i+1}^{n+1} - T_{i}^{n+1} )+\dots
```
So the coefficient multiplying ``T_i^{n+1}`` is 

```math
-\frac{\Delta t}{\rho_i c \Delta z_i}G_{i+1/2} = - \alpha_{i+1/2} = - \frac{\Delta t}{\rho_i c \Delta z_i} \frac{2(K_i \Delta z_i + K_{i+1}\Delta z_{i+1})}{(\Delta z_{i} + \Delta z_{i+1})^2} = c_i,
```

which is the superdiagonal entry in the matrix row ``i``. Similarly, the subdiagonal is given by 

```math
a_i = -\alpha_{i-1/2} = \frac{\Delta t}{\rho_i c \Delta z_i} \frac{2(K_i \Delta z_i + K_{i+1}\Delta z_{i+1})}{(\Delta z_{i} + \Delta z_{i+1})^2}.
```

The diagonal is set to conserve the stencel, i.e., ``b_i = 1 - a_i - c_i = 1 + \alpha_{i-1/2} + \alpha_{i+1/2}``.


So for every layer ``i``, we have

```math 
a_iT^{n+1}_{i-1} + b_i T^{n+1}_i +c_iT_{i+1}^{n+1} = r_i, 
```

with ``r_i = T^n_i + \mathrm{sources}``. In matrix form, we solve the tridiagonal matrix: 

```math
\begin{bmatrix}
b_1 & c_1 & 0   & \cdots & 0 \\
a_2 & b_2 & c_2 & \ddots & \vdots \\
0   & a_3 & b_3 & \ddots & 0 \\
\vdots & \ddots & \ddots & \ddots & c_{n-1} \\
0 & \cdots & 0 & a_n & b_n
\end{bmatrix}
\begin{bmatrix}
T_1^{n+1}\\
T_2^{n+1}\\
\vdots\\
T_n^{n+1}
\end{bmatrix}
=
\begin{bmatrix}
r_1\\
r_2\\
\vdots\\
r_n
\end{bmatrix}.
```


## Melt-Point Constraint

If the solved surface temperature exceeds ``T_0=0^\circ\mathrm{C}``, the surface temperature is clamped to ``T_0`` and the excessive heat is stored to calculate the the energy used to bring the surface to melt. The solver is then run an additional time to avoid unphysical temperature fluxes into deeper snow layers.  


## Thermal Conductivity Options

`diff_model` selects:

- `1`: Yen (1981)
- `2`: Sturm (1997) piecewise
- `3` (or other): Van Dusen (1929)

## API

```@docs; canonical=false
go_energy_flux!
```
