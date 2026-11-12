# Energy Flux Analytical Validation Note

This page summarizes the analytical checks of the energy solve with the Robin
surface boundary (see [Energy Balance](../processes/energy.md)). They are
executed in `test/test_longwave_consistency.jl` and `test/test_case_api.jl`.

## Insulated Equilibrium

With all surface fluxes and emission switched off and a uniform temperature,
the solve must leave interface and layer temperatures unchanged and diagnose
no melt. This is checked for a regular (100 kg m⁻²) and a vanishingly thin
(10⁻⁸ kg m⁻²) surface layer: eliminating the interface temperature must not
create a numerically singular top row, so the diagonal stays exactly one.

## Graybody Radiative Equilibrium

A surface at `T0` under incident longwave ``\sigma T_0^4`` from an isothermal
blackbody environment absorbs and emits the same graybody flux regardless of
emissivity, so neither snow nor bare ice may warm or melt. The original BESSI
treatment (`seb_scheme=:bessi`, unit incident absorption) is checked to retain
its spurious net gain ``(1-\epsilon)\,\sigma T_0^4``.

## Interface Conductance

For two layers under pure diffusion, the heat exchanged between them must
follow the harmonic interface conductance
``G_{i+1/2} = 2K_iK_{i+1}/(K_{i+1}\Delta z_i + K_i\Delta z_{i+1})``.

## Melting Boundary

When the surface energy balance would raise the interface above `T0`:

- melt is diagnosed even if the centre of the top cell is still below `T0`;
- the conduction ``G_s(T_0 - T_1)`` from the melting surface into the snow
  below is retained, and the melt energy is the remainder of the surface
  energy balance;
- the sensible-energy change of the column plus the melt energy equals the
  energy supplied at the surface.

## Ice Substrate

On a bare-ice column the public `go_energy_flux!` solves the ice substrate: a
negative surface energy balance cools both the interface and the top
substrate layer below `T0` instead of being discarded.
