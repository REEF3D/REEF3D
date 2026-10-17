# CORAL: core/

Core of CORAL: the data types and the multibody machinery every other folder builds on.

## Content

* `coral_types.h`: Eigen types and the index type used across CORAL
* `coral_version.h/.cpp`: name and version
* `coral_rotation.h/.cpp`: exponential and logarithm on SO(3), tangent operator, Euler angles
* `coral_body.h`: rigid body (mass, inertia, configuration x and R, velocities v global and w body frame)
* `coral_load.h`: the load interface (eval for trial states, commit per converged step) that gravity, structure, hydro, mooring, PTO, aero, soil and ice implement; optional analytic tangent (`has_tangent`, `tangent`: cK dF/dq + cC dF/dnu + cM dF/dnu' as triplets), such loads are left out of the finite-difference Jacobian
* `coral_constraint.h/.cpp`: constraints and joints: fixed, spherical, revolute, prismatic, cylindrical, planar (joint coordinate and rate, actuator force and torque in the joint)
* `coral_model.h/.cpp`: bodies, generalized coordinates (scalar degrees of freedom such as modal amplitudes, after the 6 per body), point masses at offsets on bodies, constraints, loads; generalized forces, inertia forces, constraint Jacobian; configuration save/set for the solvers
* `coral_points.h/.cpp`: material points of the structure (`point_set`: position, velocity, acceleration, forces back to the generalized coordinates) for distributed loads; `body_points` for points on bodies and frame nodes
* `coral_force.h/.cpp`: constant point force and moment on a body (optional ramp)
* `coral_linsolve.h/.cpp`: the Newton systems as triplets: dense LU up to 400 unknowns, sparse LU (COLAMD) above, symbolic analysis kept while the pattern stays
* `coral_timeint.h/.cpp`: Lie group generalized-α on R3 x SO(3) with index-3 constraints and Lagrange multipliers, Newton with the analytic tangents of the loads that have one plus a coloured finite-difference Jacobian of the rest (modified Newton: the Jacobian is kept while Newton converges fast); the generalized coordinates with the same scheme on R^n
* `coral_gravity.h/.cpp`: gravity on all bodies and point masses
* `coral_track.h/.cpp`: motion of bodies recorded at the ends of the (sub)steps and interpolated in between (cubic Hermite in position, velocity, rotation and angular velocity), for a coupled solver that needs the bodies at its own times (the stages of the fluid)
* `coral_static.h/.cpp`: static equilibrium (all loads at zero velocity, constraints by Lagrange multipliers): full Newton, falling back to Newton damped by the natural monotonicity test (Deuflhard) from the start configuration
* `coral_fdjac.h`: finite-difference Jacobian by columns grouped with a distance-2 colouring of the body coupling graph (`model::coloring`): a load declares which bodies it couples (`load::coupling()`: global, local to its bodies, or diagonal), so a frame of n nodes needs a few dozen residual evaluations instead of 6n

## Planned

* nodes with 3 DOF (lumped-mass lines) and flexible bodies in the same model
* sparse assembly of the constraint Jacobian (many joints)
* sub-cycling inside the fluid time step, predictor for the fluid coupling

## May include from

nothing outside `core`
