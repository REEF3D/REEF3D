# CORAL: pto/

Power take-offs of wave energy converters.

## Content

* `coral_pto_element.h`: element interface (state q, qd, t; force and its Jacobians dF/dq, dF/dqd;
  absorbed power) and the composite of elements = one PTO, with the generator gain (declutching),
  force and power limits
* `coral_pto_elements.h/.cpp`: linear damper, spring-damper with preload (also negative stiffness,
  reactive control), Coulomb friction (regularised with tanh), end stops (penalty spring-damper that only pushes), latch (stiff spring-damper at the
  latched position); the damper and the Coulomb friction can also be drive-train losses
  (`loss()`: they act on the body in every state, are not scaled or limited with the generator,
  and their power goes to `P_loss`, not to the absorbed power)
* `coral_pto_joint.h/.cpp`: prismatic joint of one rigid body to the ground (axis fixed in space,
  attachment point on the body): q, qd and the generalised direction J
* `coral_pto_controller.h/.cpp`: latching and declutching controllers (zero crossings of qd with the
  crossing time interpolated), hydrodynamic table A(w), B(w) along the joint and the optimal passive
  damping B_opt
* the joint, elements and controllers serve the REEF3D 6DOF bodies of FNPF and NHFLOW (ctrl.txt
  X 500 - X 508, X 510 losses, `reef3d/coral_sixdof_pto.cpp`; formerly `6DOF_pto_*` in `src/`)
* `coral_pto.h/.cpp`: PTO in a joint of the CORAL model (axial force or torque about the axis), built
  from the elements and controllers above: spring-damper with preload, Coulomb friction, force and
  power limits on the generator, end stops, latching (latch element with damper and optional spring)
  and declutching; absorbed power and energy of the generator; drive-train losses `c_loss`,
  `coulomb_loss` (dissipated power and energy in the channels `_Ploss`, `_Eloss`)

* `coral_joint_drive.h/.cpp`: prescribed motion of a joint coordinate (`drive`), a sine with a start ramp imposed by a stiff spring-damper; forced oscillations of a body in a resolved fluid (added mass and damping for a hydrodynamic database)
* `coral_joint_force.h/.cpp`: prescribed force along a joint axis from a table `t F` (`jforce`), with scale, start ramp and periodic repetition; system identification of a WEC (a multisine force of the PTO motor gives the intrinsic admittance of the relative motion)
* `coral_joint_drive` also from a table `t q` (`drive ... file`): a C1 cubic Hermite curve through the points, scale, ramp, periodic (measured or irregular forced motions, WEC-Sim's PTO actuation motion)
* `coral_pto.h/.cpp` between two points (link `pto <name> a <body> <x y z> b <body> <x y z>`): q the change of the distance, the force along the line; the geometry converts the motion exactly (a cylinder between a flap and the ground, a rod, a tether)
* `coral_pto_hydraulic.h/.cpp`: hydraulic PTO as PTO-Sim: double-acting cylinder, four check valves (smooth opening), HP and LP gas accumulators (polytropic), motor with leakage and torque efficiency, shaft inertia and friction, generator torque `b_gen w`, relief valve; stiff, integrated implicitly (two-stage SDIRK, sub-steps) from the accepted state along the piston path of each force call, so the force is a function of the step's state; energies of the piston, the generator shaft and the dissipation (balance checked in the tests); its own element (not scaled or limited with the generator, absorbed power)
* `coral_generator.h/.cpp`: electrical side (wave-to-wire): gear or motion conversion (ratio, efficiency), ideal machine, efficiency map, permanent-magnet machine with field-oriented control (copper and iron losses, a held force costs power, motoring pays the losses), power electronics; for the direct drive and the hydraulic motor's shaft; channels `_Pe`, `_Ee`
* `coral_pto_wec_control.h/.cpp`: further controllers that set the gains and the force of the generator once per step: gain schedule (`schedule`, table `t k c [F0]`), a user library (`user`, C interface as DISCON, optional save and load of its state), model predictive control (`mpc`: single-DOF model with the radiation memory as states, the excitation predicted from the hydro database or a table, force limit and a soft stroke limit, accelerated projected gradient for the box); tuned passive and reactive gains at input (`tune`, `pto_reactive` in `coral_pto_controller.h`)

## Planned

* the hydraulic PTO, the generator and the further controllers for the REEF3D 6DOF bodies (X 5xx); for now in CORAL's PTO (also when CORAL drives the 6DOF bodies, `fluid_body`)
* multi-DOF MPC (several bodies, coupled modes), Kalman estimation of the excitation from measured motions

## May include from

`core`, `structure`
