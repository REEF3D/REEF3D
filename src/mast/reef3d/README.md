# MAST: reef3d/

The link between MAST and REEF3D. The only folder that includes REEF3D headers and MPI.

## Planned content

* `mast_...` classes in REEF3D style: the driver hook (`start`, `ini`) and a void version when MAST is off
* fluid interfaces for CFD, NHFLOW and FNPF: kinematics and pressure at points, loads on the immersed body surface
* body motion back to the fluid (6DOF geometry, level set, immersed boundary)
* strong coupling with the fluid: sub-iterations (Aitken, IQN-ILS), added-mass estimates
* actuator line and actuator disk in the air phase
* adapters to existing REEF3D code until it moves into MAST: `beam.h`, `mooring_*`, `6DOF_pto_*`, `6DOF_obj_*`
* MPI: the structural system is solved on every rank after one reduction of the loads

## May include from

everything in MAST, REEF3D headers, MPI
