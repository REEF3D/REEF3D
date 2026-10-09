#!/bin/bash
# KP505 resolved open water, Betzy: 300-step test in ./test (grid generation, stability of the
# direct forcing when the blades first sweep cells, time per step). sbatch job_test.slurm
#SBATCH --job-name=kp505_test
#SBATCH --account=nnXXXXk
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --time=01:00:00
#SBATCH --output=test_%j.out

set -e
module purge
module load foss/2023b          # GCC + OpenMPI; the toolchain REEF3D and DIVEMesh were built with

mkdir -p test
cp control.txt floating.stl REEF3D DiveMESH test/
sed '/^N 45/d' ctrl.txt > test/ctrl.txt
echo "N 45 300" >> test/ctrl.txt
cd test

./DiveMESH > divemesh.log 2>&1
srun ./REEF3D > reef3d.log 2>&1 || true

grep -a "knox" divemesh.log
grep -a "FW-H:" reef3d.log
grep -a "timer per step" reef3d.log | tail -1
grep -a "EMERGENCY" reef3d.log | head -3 || true
