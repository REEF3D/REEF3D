# CORAL: structure/

Structural elements.

## Content

* `coral_spring.h/.cpp`: spring-damper between two points (tether: tension-only, pretension); the rigid body itself is in `core/`
* `coral_linear6.h/.cpp`: linear 6x6 stiffness and damping against the ground about a reference point (linearized mooring, foundation matrix)
* `coral_beam.h/.cpp`: co-rotational 3D beam element between two bodies (frame nodes or rigid bodies, with offsets): element frame from the chord and the mean of the end triads, local Timoshenko stiffness (EA, EI about both axes, GJ, shear GA; GA = 0 gives Euler-Bernoulli), stiffness-proportional damping β; strain matrix B analytic (spin of the element frame, inverse left Jacobian of SO(3)), analytic tangent (material and geometric stiffness by differences of -Bᵀf, damping -β BᵀKB); `beam_section` with `tube()` (thin-walled κ = 0.5) and `blend()` for tapered members
* `coral_frame.h/.cpp`: members of n elements between two points, tapered or not; nodes are 6-DOF bodies named `<member>:<i>` and shared between members that meet; ends on the ground or on a rigid body (its half element becomes a point mass of the body); lumped mass (translation and rotary inertia per node)
* `coral_superelement.h/.cpp`: Craig-Bampton superelement of linear Timoshenko members (consistent mass, rotary inertia): ends on the ground clamped, ends on rigid bodies are the interface (rigid transition piece), Guyan shapes plus n fixed-interface modes as generalized coordinates, consistent gravity, modal and stiffness-proportional damping, analytic tangent; the nodes are material points (recovered motion for Morison loads, VTP output); fixed reference frame (fixed-bottom substructures)

## Planned

* nonlinear springs and dampers
* full 6x6 cross-section stiffness and mass (composite blades, bend–twist coupling), consistent mass
* modal body and superelements in a floating frame of reference (floating substructures, fast blades and towers)
* drivetrain: shaft torsion, gearbox ratio, generator inertia
* rotor assembly: hub, blades, pitch bearings, nacelle, yaw bearing

## May include from

`core`
