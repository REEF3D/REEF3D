# CORAL: mooring/

Mooring lines and cables.

## Content

* `coral_catenary.h/.cpp`: elastic catenary with a seabed portion and friction (fairlead forces from the fairlead position, profile, anchor forces)
* lines of several sections (chain-wire-chain, lazy-wave buoyancy sections): junction points without mass (`mooring_junctions`: their start on the catenary of the whole line); a floating line's catenary is the catenary of its net buoyancy upside down
* bending stiffness of cables and risers (`EI`: energy EI/l (1 - cos phi) at each inner node, clamped ends on bodies with their moment, pinned otherwise)
* bathymetry (`coral_seabed.h/.cpp`: plane, grid of heights; or the resolved fluid's bed at the nodes)
* line failure (`mooring_break`: at a time or a tension of an end segment; the end lets go of its body or point and moves free, checkpoints replay the failures)
* the seabed of a line (`default_seabed`): at its ground end (the lower of two); lines between two bodies (shared lines, bridles) have none unless given
* `coral_mooring.h/.cpp`: quasi-static line (`mooring_qs`) and the mooring system (`mooring_system`: dynamic lumped-mass lines whose ends sit on bodies or on free points — bridle nodes, clump weights, buoys — advanced together; `mooring_line`: a system of one line; dynamics in `coral_mooring.cpp`, set-up, outputs and checkpoints in `coral_mooring_setup.cpp`) with the dynamic lines (axial stiffness and damping, normal and tangential drag and added mass, waves and current (kinematics at the converged nodes at both ends of a step, linear in between; `exact_kinematics` at every stage), smooth free-surface crossing (submerged fraction per node), seabed contact and friction (across and along the line, elliptic Coulomb limit), explicit RK4 sub-steps inside the multibody step with the ends on Hermite curves, the end nodes' mass and added mass carried by the bodies, start from the catenary, settling in the water at rest, stop with the line's name when unstable)
* `coral_mooring_static.cpp`: static equilibrium of the dynamic line's segments (`mooring_system::equilibrium`: damped Newton on the inner nodes, block-tridiagonal Jacobian by coloured finite differences), the start of every dynamic line (`no_static` in the input: from the catenary); for systems with points first the points on the catenaries of their lines (quasi-static), then Newton on the sparse Jacobian of all nodes and points, with dynamic relaxation (kinetic damping) before it if Newton does not converge

## Planned

* winches and line length changes
* geometrically exact rods (torsion, shear) for risers and umbilicals in `structure/`

## May include from

`core`, `structure`, `hydro` (wave kinematics)
