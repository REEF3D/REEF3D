#!/bin/bash
# REEF3D coupled run (B 92 62): the models run side by side in one MPI job, app numbers in this order: coarse fine
set -e
cd "$(dirname "$0")"
(cd coarse && DiveMESH > divemesh.log)
(cd fine && DiveMESH > divemesh.log)
mpirun -np 1 -wdir coarse REEF3D : -np 1 -wdir fine REEF3D > run.log 2>&1
