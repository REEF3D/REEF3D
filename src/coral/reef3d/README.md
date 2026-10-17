# CORAL: reef3d/

The link between CORAL and REEF3D. The only folder that includes REEF3D headers and MPI.

## Content

* `coral_reef3d.h/.cpp`: REEF3D runs the CORAL model of `coral.txt` alongside the fluid solver (ctrl.txt `Z 40 1`; CFD, NHFLOW, FNPF): rank 0 reads and broadcasts the input, every rank advances the same model in step with the fluid time (CORAL steps of at most the dt of `coral.txt` per fluid step), rank 0 writes the time series to `REEF3D_CORAL/`; the static equilibrium and the linear analysis (`modes.dat`, `mode_shapes.dat`) run before the first fluid step when `coral.txt` asks for them
* fluid bodies (`fluid_body` in `coral.txt`): `coral_reef3d` is the external dynamics of the REEF3D 6DOF
  bodies (`6DOF_external.h`, FNPF so far, `fnpf_6DOF_external.cpp`): the fluid hands over the loads of
  each body (force and moment about the CoG without the weight, the instantaneous added mass with the
  cross terms of several bodies) at the first stage of its step, with the predictor also the loads of
  the next stage at a predicted pose; CORAL steps the whole model (joints, PTO, mooring, structure) and
  returns the poses and velocities at the stage times of the fluid from the recorded motion of the step
  (`coral/core/coral_track.h`); mass and inertia come from REEF3D unless `own_mass`; the time series of
  `REEF3D_CORAL/` gets the fluid loads (`fluid_<body>_Fx` .. `Mz`, added-mass part included)
* `coral_sixdof_pto.cpp`: the PTO of the REEF3D 6DOF body of FNPF and NHFLOW (ctrl.txt X 500 - X 508;
  members of `sixdof_obj`, formerly `6DOF_obj_pto.cpp`): joint, elements and controllers from
  `pto/`, the stage times and energy weights of the body's Runge-Kutta scheme, the linearly implicit
  end stops (X 500 2), the stability check and `REEF3D_6DOF_pto_<n>.dat`

## Planned

* fluid bodies in NHFLOW and CFD (load lag of the reforcing, artificial added mass), joints between fluid bodies
* resolved fluid kinematics at CORAL points (Morison members, mooring lines)
* strong coupling with the fluid where needed: sub-iterations (Aitken, IQN-ILS)
* actuator line and actuator disk in the air phase
* adapters to existing REEF3D code until it moves into CORAL: `beam.h`, `mooring_*`, `6DOF_obj_*`
* MPI: the structural system is solved on every rank after one reduction of the loads

## May include from

everything in CORAL, REEF3D headers, MPI
