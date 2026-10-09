#!/bin/bash
# KP505 resolved open water, Betzy: production run, 8 revolutions (N 41 0.8 s), state files every
# 0.1 s for restarts (P 40 1, P 42 0.1). sbatch job_run.slurm
#SBATCH --job-name=kp505
#SBATCH --account=nnXXXXk
#SBATCH --partition=normal
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:00:00
#SBATCH --output=run_%j.out

set -e
module purge
module load foss/2023b

./DiveMESH > divemesh.log 2>&1
srun ./REEF3D > reef3d.log 2>&1
