# CORAL: pto/

Power take-offs of wave energy converters.

## Content

* `coral_pto_element.h`: element interface (state q, qd, t; force and its Jacobians dF/dq, dF/dqd;
  absorbed power) and the composite of elements = one PTO, with the generator gain (declutching),
  force and power limits
* `coral_pto_elements.h/.cpp`: linear damper, spring-damper with preload (also negative stiffness,
  reactive control), Coulomb friction (regularised with tanh), end stops (penalty spring-damper that only pushes), latch (stiff spring-damper at the
  latched position)
* `coral_pto_joint.h/.cpp`: prismatic joint of one rigid body to the ground (axis fixed in space,
  attachment point on the body): q, qd and the generalised direction J
* `coral_pto_controller.h/.cpp`: latching and declutching controllers (zero crossings of qd with the
  crossing time interpolated), hydrodynamic table A(w), B(w) along the joint and the optimal passive
  damping B_opt
* the joint, elements and controllers serve the REEF3D 6DOF bodies of FNPF and NHFLOW (ctrl.txt
  X 500 - X 508, `reef3d/coral_sixdof_pto.cpp`; formerly `6DOF_pto_*` in `src/`)
* `coral_pto.h/.cpp`: PTO in a joint of the CORAL model (axial force or torque about the axis), built
  from the elements and controllers above: spring-damper with preload, Coulomb friction, force and
  power limits on the generator, end stops, latching (latch element with damper and optional spring)
  and declutching; absorbed power and energy of the generator

* `coral_joint_drive.h/.cpp`: prescribed motion of a joint coordinate (`drive`), a sine with a start ramp imposed by a stiff spring-damper; forced oscillations of a body in a resolved fluid (added mass and damping for a hydrodynamic database)

## Planned

* hydraulic PTO: cylinders, accumulators, valves, motor
* electrical generator model (efficiency, delivered power)
* PTO between two points (linear actuator not in a joint)

## May include from

`core`, `structure`
