# CORAL: reef3d/

The link between CORAL and REEF3D. The only folder that includes REEF3D headers and MPI.

## Content

* `coral_reef3d.h/.cpp`: REEF3D runs the CORAL model of `coral.txt` alongside the fluid solver (ctrl.txt `Z 40 1`; CFD, NHFLOW, FNPF): rank 0 reads and broadcasts the input, every rank advances the same model in step with the fluid time (CORAL steps of at most the dt of `coral.txt` per fluid step), rank 0 writes the time series to `REEF3D_CORAL/`; the static equilibrium and the linear analysis (`modes.dat`, `mode_shapes.dat`) run before the first fluid step when `coral.txt` asks for them
* `coral_sixdof_pto.cpp`: the PTO of the REEF3D 6DOF body of FNPF and NHFLOW (ctrl.txt X 500 - X 508;
  members of `sixdof_obj`, formerly `6DOF_obj_pto.cpp`): joint, elements and controllers from
  `pto/`, the stage times and energy weights of the body's Runge-Kutta scheme, the linearly implicit
  end stops (X 500 2), the stability check and `REEF3D_6DOF_pto_<n>.dat`

## Planned

* fluid interfaces for CFD, NHFLOW and FNPF: kinematics and pressure at points, loads on the immersed body surface
* body motion back to the fluid (6DOF geometry, level set, immersed boundary)
* strong coupling with the fluid: sub-iterations (Aitken, IQN-ILS), added-mass estimates
* actuator line and actuator disk in the air phase
* adapters to existing REEF3D code until it moves into CORAL: `beam.h`, `mooring_*`, `6DOF_obj_*`
* MPI: the structural system is solved on every rank after one reduction of the loads

## May include from

everything in CORAL, REEF3D headers, MPI
