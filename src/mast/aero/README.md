# MAST: aero/

Aerodynamic loads on rotors, towers and nacelles.

## Planned content

* airfoil polars: interpolation in angle of attack, Reynolds number and thickness
* blade element momentum: tip and hub losses, high-induction correction, skewed wake
* dynamic inflow (Øye) and dynamic stall (Beddoes–Leishman, Øye)
* tower shadow and tower drag, nacelle and hub drag
* later: free vortex wake
* the actuator line in REEF3D's resolved air phase uses these loads through `reef3d/`

## May include from

`core`, `structure`, `wind`
