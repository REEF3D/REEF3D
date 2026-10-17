# CORAL: reef3d/

The link between CORAL and REEF3D. The only folder that includes REEF3D headers and MPI.

## Content

* `coral_reef3d.h/.cpp`: (a thin wrapper around `coral/io/coral_coupling`, which does all the coupling work and can record it) REEF3D runs the CORAL model of `coral.txt` alongside the fluid solver (ctrl.txt `Z 40 1`; CFD, NHFLOW, FNPF): rank 0 reads and broadcasts the input, every rank advances the same model in step with the fluid time (CORAL steps of at most the dt of `coral.txt` per fluid step), rank 0 writes the time series to `REEF3D_CORAL/`; the static equilibrium and the linear analysis (`modes.dat`, `mode_shapes.dat`) run before the first fluid step when `coral.txt` asks for them
* fluid bodies (`fluid_body` in `coral.txt`): `coral_reef3d` is the external dynamics of the REEF3D 6DOF
  bodies (`6DOF_external.h`, FNPF so far, `fnpf_6DOF_external.cpp`): the fluid hands over the loads of
  each body (force and moment about the CoG without the weight, the instantaneous added mass with the
  cross terms of several bodies) at the first stage of its step, with the predictor also the loads of
  the next stage at a predicted pose; CORAL steps the whole model (joints, PTO, mooring, structure) and
  returns the poses and velocities at the stage times of the fluid from the recorded motion of the step
  (`coral/core/coral_track.h`); mass and inertia come from REEF3D unless `own_mass`; the time series of
  `REEF3D_CORAL/` gets the fluid loads (`fluid_<body>_Fx` .. `Mz`, added-mass part included)
* one-way loads (FNPF and NHFLOW): `fluid_body ... one_way` puts the loads of a 6DOF body that REEF3D
  holds fixed (X 11 0) or moves itself on any CORAL body, at the REEF3D CoG; `fluid_member` distributes
  them along a flexible member: the 6DOF body bins its hull force in strips along the member axis
  (`sixdof_obj::external_strips_*`, filled by `forces_fnpf` and `force_calc_stl`: the pressure of each
  triangle of the load integration, linear through its quadrature values, is integrated exactly over
  the parts of the triangle in each strip, so the strips do not depend on the size of the STL
  triangles), reduced over the ranks. FNPF delivers the loads of the first stage (t^n), so CORAL runs one fluid step behind; NHFLOW
  at the end of the step (`sixdof_nhflow::ext_loads`)
* resolved kinematics (`fluid_kinematics`, FNPF and NHFLOW): `coral_fluid_sampler.h/.cpp` interpolate
  velocity and free surface of the fluid grid at the CORAL points (FNPF: sigma grid, the top cell
  sampled just below the surface; NHFLOW: cell centres), once per fluid step, reduced over the ranks
* record and replay: `fluid_coupling record <file>` writes every fluid input into
  `REEF3D_CORAL/<file>`; the standalone `CORAL` replays it without the fluid (a changed structure
  with the same fluid points and strips: a tower, an RNA, a controller study on the same loads)
* `coral_sixdof_pto.cpp`: the PTO of the REEF3D 6DOF body of FNPF and NHFLOW (ctrl.txt X 500 - X 508;
  members of `sixdof_obj`, formerly `6DOF_obj_pto.cpp`): joint, elements and controllers from
  `pto/`, the stage times and energy weights of the body's Runge-Kutta scheme, the linearly implicit
  end stops (X 500 2), the stability check and `REEF3D_6DOF_pto_<n>.dat`

## Planned

* two-way fluid bodies in NHFLOW and CFD (load lag of the reforcing, artificial added mass), one-way CFD
* strong coupling with the fluid where needed: sub-iterations (Aitken, IQN-ILS)
* actuator line and actuator disk in the air phase
* adapters to existing REEF3D code until it moves into CORAL: `beam.h`, `mooring_*`, `6DOF_obj_*`
* MPI: the structural system is solved on every rank after one reduction of the loads

## May include from

everything in CORAL, REEF3D headers, MPI
