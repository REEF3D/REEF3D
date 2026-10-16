# CORAL: core/

Core of CORAL: the data types and the multibody machinery every other folder builds on.

## Content

* `coral_types.h`: Eigen types and the index type used across CORAL
* `coral_version.h/.cpp`: name and version
* `coral_rotation.h/.cpp`: exponential and logarithm on SO(3), tangent operator, Euler angles
* `coral_body.h`: rigid body (mass, inertia, configuration x and R, velocities v global and w body frame)
* `coral_load.h`: the load interface (eval for trial states, commit per converged step) that gravity, structure, hydro, mooring, PTO, aero, soil and ice implement
* `coral_constraint.h/.cpp`: constraints and joints: fixed, spherical, revolute, prismatic, cylindrical, planar (joint coordinate and rate, actuator force and torque in the joint)
* `coral_model.h/.cpp`: bodies, constraints, loads; generalized forces, inertia forces, constraint Jacobian
* `coral_timeint.h/.cpp`: Lie group generalized-α on R3 x SO(3) with index-3 constraints and Lagrange multipliers, Newton with a finite-difference Jacobian (modified Newton: the Jacobian is kept while Newton converges fast)
* `coral_gravity.h/.cpp`: gravity on all bodies

## Planned

* nodes with 3 DOF (lumped-mass lines) and flexible bodies in the same model
* analytic tangents of the elements (sparse Jacobian for large models)
* sub-cycling inside the fluid time step, predictor for the fluid coupling

## May include from

nothing outside `core`
