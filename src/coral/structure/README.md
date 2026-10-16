# CORAL: structure/

Structural elements.

## Content

* `coral_spring.h/.cpp`: spring-damper between two points (tether: tension-only, pretension); the rigid body itself is in `core/`
* `coral_linear6.h/.cpp`: linear 6x6 stiffness and damping against the ground about a reference point (linearized mooring, foundation matrix)
* `coral_beam.h/.cpp`: co-rotational 3D beam element between two bodies (frame nodes or rigid bodies, with offsets): element frame from the chord and the mean of the end triads, local Timoshenko stiffness (EA, EI about both axes, GJ, shear GA; GA = 0 gives Euler-Bernoulli), stiffness-proportional damping β; `beam_section` with `tube()` (thin-walled κ = 0.5) and `blend()` for tapered members
* `coral_frame.h/.cpp`: members of n elements between two points, tapered or not; nodes are 6-DOF bodies named `<member>:<i>` and shared between members that meet; ends on the ground or on a rigid body; lumped mass (translation and rotary inertia per node)

## Planned

* nonlinear springs and dampers
* full 6x6 cross-section stiffness and mass (composite blades, bend–twist coupling), consistent mass
* analytic element tangents (Phase 2.3)
* modal body: assumed modes in a floating frame of reference (fast blades and towers)
* Craig-Bampton superelements
* drivetrain: shaft torsion, gearbox ratio, generator inertia
* rotor assembly: hub, blades, pitch bearings, nacelle, yaw bearing

## May include from

`core`
