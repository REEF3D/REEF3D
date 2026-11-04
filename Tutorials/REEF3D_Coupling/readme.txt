REEF3D runtime coupling (B 92 62): FNPF drives CFD, NHFLOW or FNPF models that run in the same MPI job.
Each case has one folder per model and a run.sh that runs DIVEMesh in every folder and starts all
models with one mpirun command (one app context per model, -wdir sets its folder). The explanation is
in the REEF3D User Guide, Chapter "Tutorial | Coupled Models".

1_FNPF_CFD_2D             2D FNPF wave tank -> 2D CFD box (one-way)
2_FNPF2D_CFD3D_oblique    2D FNPF -> 3D CFD box, waves at 30 deg (B 105), relaxation zones on all edges (B 520)
3_FNPF_FNPF_nesting       coarse FNPF (dx 0.2 m) -> fine FNPF (dx 0.05 m)
4_FNPF_CFD_floating_body  floating box: the CFD body (6DOF) drives the FNPF body (I 243)
5_FNPF_CFD_two-way        two-way: CFD with a wall at x+, the reflection enters FNPF (I 244)
6_FNPF_FNPF_CFD_chain     regional FNPF -> local FNPF -> CFD (I 246)
