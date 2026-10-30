# CORAL: mooring/

Mooring lines and cables.

## Content

* `coral_catenary.h/.cpp`: elastic catenary with a seabed portion and friction (fairlead forces from the fairlead position, profile, anchor forces)
* lines of several sections (chain-wire-chain, lazy-wave buoyancy sections): junction points without mass (`mooring_junctions`: their start on the catenary of the whole line); a floating line's catenary is the catenary of its net buoyancy upside down
* bending stiffness of cables and risers (`EI`: energy EI/l (1 - cos phi) at each inner node, clamped ends on bodies with their moment, pinned otherwise)
* bathymetry (`coral_seabed.h/.cpp`: plane, grid of heights; or the resolved fluid's bed at the nodes)
* line failure (`mooring_break`: at a time or a tension of an end segment; the end lets go of its body or point and moves free, checkpoints replay the failures)
* the seabed of a line (`default_seabed`): at its ground end (the lower of two); lines between two bodies (shared lines, bridles) have none unless given
* `coral_mooring.h/.cpp`: quasi-static line (`mooring_qs`) and the mooring system (`mooring_system`: dynamic lumped-mass lines whose ends sit on bodies or on free points — bridle nodes, clump weights, buoys — advanced together; `mooring_line`: a system of one line; dynamics in `coral_mooring.cpp`, set-up, outputs and checkpoints in `coral_mooring_setup.cpp`) with the dynamic lines (axial stiffness and damping: only while a segment is taut, or with `BA_slack` in every segment as MoorDyn, normal and tangential drag and added mass, waves and current (kinematics at the converged nodes at both ends of a step, linear in between; `exact_kinematics` at every stage), smooth free-surface crossing (submerged fraction per node), seabed contact and friction (across and along the line, elliptic limit; stick-slip: anchors on shear springs up to the static limit, sliding at the kinetic one, a line on a slope stays put; or the regularized Coulomb law with `v_eps`; the anchors in `coral_mooring_friction.cpp`), explicit RK4 sub-steps inside the multibody step with the ends on Hermite curves, the end nodes' mass and added mass carried by the bodies, start from the catenary, settling in the water at rest, stop with the line's name when unstable)
* `coral_mooring_friction.cpp`: seabed friction of the dynamic lines and points (`seabed_friction`, `stick_friction`: anchors on shear springs, stick up to `mu_static` times the kinetic limit, slide on the closest point of the kinetic ellipse; the regularized law with `v_eps` > 0), in its own file for the compiler's inlining of the sub-step loop
* MoorDyn v2 options of the dynamic lines: `BA` < 0 as the damping ratio of a segment's axial mode (`BA/-zeta`),
  nonlinear elasticity from a stress-strain table (`mooring_props::EA_strain`, `EA_tension`; input `EA_curve`), a
  current profile over the height (`mooring_system::set_current_profile`; input `current_profile`), outputs at the
  nodes (`set_line_outputs`: node positions, segment tensions; input `outputs p|t`)
* `coral_mooring_static.cpp`: static equilibrium of the dynamic line's segments (`mooring_system::equilibrium`: damped Newton on the inner nodes, block-tridiagonal Jacobian by coloured finite differences), the start of every dynamic line (`no_static` in the input: from the catenary); for systems with points first the points on the catenaries of their lines (quasi-static), then Newton on the sparse Jacobian of all nodes and points, with dynamic relaxation (kinetic damping) before it if Newton does not converge

## Validation

`Regression/moordyn/`: the dynamic lines, line networks (bridle, clump weight, buoy, chain-wire-chain) and a
moored floating body against MoorDyn v2.7.2 (MoorDyn-C, prescribed fairlead motion and coupled planar body), in
still water, currents and waves (CORAL's wave field at MoorDyn's nodes); against published results: the flexible
jumper of Lowe and Langley (2006, top tension in motion and waves) and the chain of Azcona et al. (2017, measured
static tension, the ratios of the dynamic peaks); `Regression/unit/coral_moordyn_test.cpp` checks CORAL against
the stored MoorDyn time series and the published values. Known differences:
MoorDyn's seabed damping also pulls a node that rises from the seabed (CORAL's contact only pushes: snap loads of
a chain lifting off differ a little), MoorDyn's points take their whole volume above the free surface (CORAL: the
submerged part), MoorDyn's initial condition leaves the current out (CORAL's static start takes it), MoorDyn holds
the last tension of a stress-strain table beyond it (CORAL: the last slope).

## Planned

* winches and line length changes
* geometrically exact rods (torsion, shear) for risers and umbilicals in `structure/`

## May include from

`core`, `structure`, `hydro` (wave kinematics)
