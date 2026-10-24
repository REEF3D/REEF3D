# CORAL: aero/

Aerodynamic loads on rotors, towers and nacelles.

## Content

* `coral_airfoil.h/.cpp`: airfoil polars (cl, cd, cm against the angle of attack, linear, wrapped around for a -180 .. 180 deg table, held at the ends otherwise), several tables per airfoil blended in log(Re); readers of AeroDyn's airfoil files (AirfoilInfo v1: NumTabs, Re, the unsteady aerodynamics parameters kept by name for the dynamic stall model, NumAlf rows) and of plain tables
* `coral_bem.h/.cpp`: blade element momentum theory: one element solved with Ning's single-equation method (guaranteed bracket, Brent; as AeroDyn 15), Prandtl tip and hub loss, Buhl's high-induction correction, tangential induction, drag in the induction optional, the propeller brake region, Pitt-Peters skewed-wake factor; steady rotor in axial inflow (thrust, torque, power, Cp, Ct by the trapezoid rule over the blade stations); reader of AeroDyn 15 blade files

## Planned

* the aerodynamic load on the blades of the multibody model (rotating, flexible blades, platform motion, yaw and tilt; Phase 3.3), tower shadow and tower and nacelle drag
* dynamic inflow (Øye) and dynamic stall (Beddoes-Leishman as in AeroDyn, Øye)
* later: free vortex wake
* the actuator line in REEF3D's resolved air phase uses these loads through `reef3d/`

## May include from

`core`, `structure`, `wind`
