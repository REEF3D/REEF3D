# CORAL: linear/

Linear analysis of a CORAL model.

## Content

* `coral_modes.h/.cpp`: linearisation of the model about its current configuration (M, C, K by finite differences, including the geometric stiffness of the constraint forces and the loads after a static solution), reduction to the null space of the constraint Jacobian, undamped modes (K φ = ω² M φ) or damped modes (first-order form, damping ratios)

## Planned

* linearisation about a moving operating point (rotating rotor), aero, hydro and control in the linear model
* Campbell diagrams over rotor speed
* multi-blade coordinate transform
* linear state-space models for controller design

## May include from

`core`, `structure`
