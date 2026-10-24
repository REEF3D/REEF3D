# CORAL: aero/

Aerodynamic loads on rotors, towers and nacelles.

## Content

* `coral_airfoil.h/.cpp`: airfoil polars (cl, cd, cm against the angle of attack, linear, wrapped around for a -180 .. 180 deg table, held at the ends otherwise), several tables per airfoil blended in log(Re); readers of AeroDyn's airfoil files (AirfoilInfo v1: NumTabs, Re, the unsteady aerodynamics parameters kept by name for the dynamic stall model, NumAlf rows) and of plain tables
* `coral_bem.h/.cpp`: blade element momentum theory: one element solved with Ning's single-equation method (guaranteed bracket, Brent; as AeroDyn 15), Prandtl tip and hub loss, Buhl's high-induction correction, tangential induction, drag in the induction optional, the propeller brake region, Pitt-Peters skewed-wake factor; steady rotor in axial inflow (thrust, torque, power, Cp, Ct by the trapezoid rule over the blade stations); reader of AeroDyn 15 blade files

* `coral_rotor_aero.h/.cpp`: the aerodynamic load on the blades of a rotor in the multibody model: quasi-steady blade element momentum at each station from its actual motion (rotation, the nacelle's and the platform's motion, blade deflection) and orientation (the chord's angle to the rotor plane: twist, pitch and elastic torsion), the inflow of the wind field, cone and tilt by the geometry, the skewed wake (Pitt-Peters with Coleman's wake angle), tower potential flow (AeroDyn's TwrPotent 1); channels of thrust, torque, power, Cp, Ct, tip speed ratio, yaw, pitch and the aerodynamic root moments of each blade; an approximate tangent (the aerodynamic damping of each station on its own body) keeps the rotor out of the finite-difference Jacobian
* `coral_rotor.h/.cpp`: the rotor builder (`rotor` keyword): the hub in a revolute shaft joint on the nacelle, the drivetrain and generator (`structure/coral_drivetrain.h`), the stations of AeroDyn's blade file; the blades rigid on the hub (their distributed mass from ElastoDyn's blade file as point masses, collective pitch in the aerodynamics), or on pitch bearings with second-order actuators (`structure/coral_pitch_actuator.h`), rigid or flexible: co-rotational beam elements from the structural table (flap and edge stiffness along the structural twist, mass), the aerodynamic stations on the nearest node

* `coral_rotor_tune.h/.cpp`: the steady operating curve of a rotor from `rotor_steady` (Cp_max and lambda_opt, the region-2 gain K, the rated wind speed, the pitch above rated that holds the rated power at the rated speed, the power sensitivity dP/dtheta and its linear fit for the gain scheduling): the input of the controller's tuning (`controller ... tune`)

## Planned

* tower and nacelle drag; full 6x6 blade sections (bend-twist coupling, offsets of the elastic axis and the centre of mass); modal blades for speed
* dynamic inflow (Øye) and dynamic stall (Beddoes-Leishman as in AeroDyn, Øye)
* later: free vortex wake
* the actuator line in REEF3D's resolved air phase uses these loads through `reef3d/`

## May include from

`core`, `structure`, `wind`
