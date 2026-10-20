# CORAL: mooring/

Mooring lines and cables.

## Content

* `coral_catenary.h/.cpp`: elastic catenary with a seabed portion and friction (fairlead forces from the fairlead position, profile, anchor forces)
* the seabed of a line (`default_seabed`): at its ground end (the lower of two); lines between two bodies (shared lines, bridles) have none unless given
* `coral_mooring.h/.cpp`: quasi-static line (`mooring_qs`) and dynamic lumped-mass line (`mooring_line`: axial stiffness and damping, normal and tangential drag and added mass, waves and current (kinematics at the converged nodes at both ends of a step, linear in between; `exact_kinematics` at every stage), smooth free-surface crossing (submerged fraction per node), seabed contact and friction (across and along the line, elliptic Coulomb limit), explicit RK4 sub-steps inside the multibody step with the ends on Hermite curves, the end nodes' mass and added mass carried by the bodies, start from the catenary, settling in the water at rest, stop with the line's name when unstable)
* `coral_mooring_static.cpp`: static equilibrium of the dynamic line's segments (`mooring_line::equilibrium`: damped Newton on the inner nodes, block-tridiagonal Jacobian by coloured finite differences), the start of every dynamic line (`no_static` in the input: from the catenary)

## Planned

* buoys and clump weights on the lines, multi-segment lines with different properties, line-to-line connections
* winches and line length changes, line breaking
* cables and risers with bending stiffness: geometrically exact rod in `structure/`

## May include from

`core`, `structure`, `hydro` (wave kinematics)
