# CORAL: hydro/

Engineering hydrodynamic loads (no resolved fluid needed).

## Content

* `coral_waves.h/.cpp`: incident waves as linear components (regular, JONSWAP, Pierson-Moskowitz with seeded phases), dispersion, kinematics and dynamic pressure, start ramp
* `coral_hydro_db.h/.cpp`: hydrodynamic database (C, A(w), B(w), A_inf, X(w, beta), reference points, volume): CORAL text format (read/write), WAMIT reader (.1, .3, .hst)
* `coral_radiation.h/.cpp`: retardation function K(t) (exact for B piecewise linear), A_inf from Ogilvie's relation, state-space realization of K (eigensystem realization, diagonal form, exact step for a linear velocity)
* `coral_hydro_load.h/.cpp`: Cummins equation for the database bodies: buoyancy, hydrostatics, excitation, radiation (convolution or state-space), linear and quadratic damping
* `coral_morison.h/.cpp`: Morison members on rigid bodies (strip theory, relative velocity, current, added mass of the member)
* `coral_member_morison.h/.cpp`: Morison loads distributed along flexible members (frame nodes or superelement nodes through `point_set`), relative motion, current, buoyancy of the submerged strips, wet part of the strip cut by the surface

## Planned

* hydrostatics on the instantaneous wetted surface
* drift coefficients and QTFs: second-order drift (Newman, full QTF)
* MacCamy–Fuchs correction, marine growth, directional spreading
* slamming and breaking-wave impact models
* nonlinear Froude–Krylov loads from the incident wave field (REEF3D FNPF/NHFLOW or analytical)
* incident wave kinematics: REEF3D's wave theories when running without a fluid grid

## May include from

`core`
