# MAST: structure/

Structural elements.

## Planned content

* rigid body (6 DOF, full inertia tensor)
* point mass and point inertia, linear and nonlinear springs and dampers between nodes
* nonlinear beam: geometrically exact, full 6x6 cross-section stiffness and mass (composite blades, bend–twist coupling); first through the adapter to REEF3D's `beam.h`
* linear 3D Timoshenko frame elements for jackets, monopiles and towers
* modal body: assumed modes in a floating frame of reference (fast blades and towers)
* Craig-Bampton superelements
* drivetrain: shaft torsion, gearbox ratio, generator inertia
* rotor assembly: hub, blades, pitch bearings, nacelle, yaw bearing

## May include from

`core`
