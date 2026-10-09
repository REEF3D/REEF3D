#!/bin/bash
# Sphere Re = 300, FW-H validation, D/40: 512 ranks. Run directory: control.txt, ctrl.txt,
# 6DOF_motion.dat (from this folder), floating.stl (Regression/cases/cfd_3d_heave_sphere_6dof),
# REEF3D and DiveMESH built on Betzy. sbatch job.slurm
#SBATCH --job-name=sphere_d40
#SBATCH --account=nnXXXXk
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --time=06:00:00
#SBATCH --output=run_%j.out

set -e
for f in control.txt ctrl.txt floating.stl 6DOF_motion.dat REEF3D DiveMESH; do
    [ -e "$f" ] || { echo "missing $f in $(pwd)"; exit 1; }
done

module purge
module load foss/2023b          # GCC + OpenMPI; the toolchain REEF3D and DIVEMesh were built with

./DiveMESH > divemesh.log 2>&1
srun ./REEF3D > reef3d.log 2>&1
