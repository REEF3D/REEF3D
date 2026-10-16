# MAST: core/

Core of MAST: the data types and the multibody machinery every other folder builds on.

## Planned content

* `mast_types.h`: Eigen types and the index type used across MAST
* `mast_version.h/.cpp`: name and version
* rotations: quaternions, rotation vectors, small-angle forms, Euler angles for output
* model: bodies, nodes, degrees of freedom and their numbering, body and node connectivity
* constraints and joints: fixed, revolute, prismatic, universal, spherical; Lagrange multipliers or penalty
* assembly of mass, damping and stiffness contributions and of the load vectors
* time integration: generalized-α (Chung–Hulbert) with Newton iterations, sub-cycling inside the fluid time step, an explicit predictor for the fluid coupling
* linear solvers: dense and sparse (Eigen) for the small structural systems
* the load interface that structure, aero, hydro, mooring, PTO, soil and ice implement

## May include from

nothing outside `core`
