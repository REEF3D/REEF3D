# CORAL: hydro/

Engineering hydrodynamic loads (no resolved fluid needed).

## Content

* `coral_waves.h/.cpp`: `wave_kinematics`, the interface of every source of wave kinematics (elevation, velocity, acceleration, dynamic pressure; the keyed calls `elevation_at` / `kinematics_at` name the load and the point, for sources that are sampled at fixed points); `wave_field`: incident waves as linear components (regular, JONSWAP, Pierson-Moskowitz with seeded phases), dispersion, kinematics and dynamic pressure, start ramp
* `coral_hydro_db.h/.cpp`: hydrodynamic database (C, A(w), B(w), A_inf, X(w, beta), reference points, volume): CORAL text format (read/write), WAMIT reader (.1, .3, .hst)
* `coral_radiation.h/.cpp`: retardation function K(t) (exact for B piecewise linear), A_inf from Ogilvie's relation, state-space realization of K (eigensystem realization, diagonal form, exact step for a linear velocity)
* `coral_hydro_load.h/.cpp`: Cummins equation for the database bodies: buoyancy, hydrostatics, excitation, radiation (convolution or state-space), linear and quadratic damping
* `coral_morison.h/.cpp`: Morison members on rigid bodies (strip theory, relative velocity, current, added mass of the member)
* `coral_body_drag.h/.cpp`: quadratic (viscous) damping of a body (`drag`): diagonal coefficients, or rho Cd A / 2 along an axis (the drag a potential-flow solver such as FNPF leaves out)
* `coral_fluid_load.h/.cpp`: loads of a resolved fluid (REEF3D, `fluid_body`) on rigid bodies: force and moment per body sampled at the fluid time levels, a polynomial in time through the last samples (a provisional sample of the next fluid stage turns the extrapolation into an interpolation, the predictor), the instantaneous added mass implicit on the accelerations of the current state; the artificial added mass with compensation for a fluid without its own (NHFLOW); one-way: the load acts at a point of the body (the REEF3D CoG, `set_attachment`) without added mass
* `coral_member_morison.h/.cpp`: Morison loads distributed along flexible members (frame nodes or superelement nodes through `point_set`), relative motion, current, buoyancy of the submerged strips, wet part of the strip cut by the surface
* `coral_sampled_kinematics.h/.cpp`: wave kinematics sampled by a resolved fluid (`fluid_kinematics`): the loads that use it (Morison on bodies and members, dynamic mooring lines) name their points (`load::kinematics_points`: strip centres, line nodes), the fluid fills velocity and elevation at these points once per fluid time level, CORAL interpolates in time (Lagrange through the last 2 or 3 levels; the acceleration is the time derivative of the same polynomial); points outside the fluid domain and the times before the first sample use the `wave_field` of the input
* `coral_member_fluid_load.h/.cpp`: the loads of a resolved body (REEF3D 6DOF body, `fluid_member`) distributed along a flexible member: the force on the hull in strips along the member axis, each strip force at its centre on the member element below it (linear weights to the element ends), a polynomial in time through the last samples like `fluid_load`; `excitation` subtracts the first sample (the static part)

## Planned

* hydrostatics on the instantaneous wetted surface
* drift coefficients and QTFs: second-order drift (Newman, full QTF)
* MacCamy–Fuchs correction, marine growth, directional spreading
* slamming and breaking-wave impact models
* nonlinear Froude–Krylov loads from the incident wave field (REEF3D FNPF/NHFLOW or analytical)
* incident wave kinematics: REEF3D's wave theories when running without a fluid grid
* sampled pressure on the wetted surface of CORAL bodies (Froude–Krylov from the resolved field)

## May include from

`core`
