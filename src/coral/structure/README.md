# CORAL: structure/

Structural elements.

## Content

* `coral_spring.h/.cpp`: spring-damper between two points (tether: tension-only, pretension); the rigid body itself is in `core/`
* `coral_linear6.h/.cpp`: linear 6x6 stiffness and damping against the ground about a reference point (linearized mooring, foundation matrix)

## Planned

* point mass and point inertia, nonlinear springs and dampers
* nonlinear beam: geometrically exact, full 6x6 cross-section stiffness and mass (composite blades, bend–twist coupling); first through the adapter to REEF3D's `beam.h`
* linear 3D Timoshenko frame elements for jackets, monopiles and towers
* modal body: assumed modes in a floating frame of reference (fast blades and towers)
* Craig-Bampton superelements
* drivetrain: shaft torsion, gearbox ratio, generator inertia
* rotor assembly: hub, blades, pitch bearings, nacelle, yaw bearing

## May include from

`core`
