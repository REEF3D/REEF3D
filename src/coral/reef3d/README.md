# CORAL: reef3d/

The link between CORAL and REEF3D. The only folder that includes REEF3D headers and MPI.

## Content

* `coral_reef3d.h/.cpp`: REEF3D runs the CORAL model of `coral.txt` alongside the fluid solver (ctrl.txt `Z 40 1`; CFD, NHFLOW, FNPF): rank 0 reads and broadcasts the input, every rank advances the same model in step with the fluid time (CORAL steps of at most the dt of `coral.txt` per fluid step), rank 0 writes the time series to `REEF3D_CORAL/`

## Planned

* fluid interfaces for CFD, NHFLOW and FNPF: kinematics and pressure at points, loads on the immersed body surface
* body motion back to the fluid (6DOF geometry, level set, immersed boundary)
* strong coupling with the fluid: sub-iterations (Aitken, IQN-ILS), added-mass estimates
* actuator line and actuator disk in the air phase
* adapters to existing REEF3D code until it moves into CORAL: `beam.h`, `mooring_*`, `6DOF_pto_*`, `6DOF_obj_*`
* MPI: the structural system is solved on every rank after one reduction of the loads

## May include from

everything in CORAL, REEF3D headers, MPI
