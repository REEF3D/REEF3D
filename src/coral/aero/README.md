# CORAL: aero/

Aerodynamic loads on rotors, towers and nacelles.

## Content

* `coral_airfoil.h/.cpp`: airfoil polars (cl, cd, cm against the angle of attack, linear, wrapped around for a -180 .. 180 deg table, held at the ends otherwise), several tables per airfoil blended in log(Re); readers of AeroDyn's airfoil files (AirfoilInfo v1: NumTabs, Re, the unsteady aerodynamics parameters kept by name for the dynamic stall model, NumAlf rows) and of plain tables
* `coral_bem.h/.cpp`: blade element momentum theory: one element solved with Ning's single-equation method (guaranteed bracket, Brent; as AeroDyn 15), Prandtl tip and hub loss, Buhl's high-induction correction, tangential induction, drag in the induction optional, the propeller brake region, Pitt-Peters skewed-wake factor; steady rotor in axial inflow (thrust, torque, power, Cp, Ct by the trapezoid rule over the blade stations); reader of AeroDyn 15 blade files

* `coral_rotor_aero.h/.cpp`: the aerodynamic load on the blades of a rotor in the multibody model: quasi-steady blade element momentum at each station from its actual motion (rotation, the nacelle's and the platform's motion), the inflow of the wind field, cone and tilt by the geometry, the skewed wake (Pitt-Peters with Coleman's wake angle), tower potential flow (AeroDyn's TwrPotent 1); channels of thrust, torque, power, Cp, Ct, tip speed ratio, yaw, pitch and the aerodynamic root moments of each blade
* `coral_rotor.h/.cpp`: the rotor builder (`rotor` keyword): hub and rigid blades as one body with the blades' distributed mass (ElastoDyn's blade file), the revolute shaft joint on the nacelle, the drivetrain and generator (`structure/coral_drivetrain.h`), the stations of AeroDyn's blade file, collective pitch constant or in time

## Planned

* flexible blades (beams from the blade's structural table) and pitch joints with actuators (Phase 3.3, second part); tower and nacelle drag
* dynamic inflow (Øye) and dynamic stall (Beddoes-Leishman as in AeroDyn, Øye)
* later: free vortex wake
* the actuator line in REEF3D's resolved air phase uses these loads through `reef3d/`

## May include from

`core`, `structure`, `wind`
