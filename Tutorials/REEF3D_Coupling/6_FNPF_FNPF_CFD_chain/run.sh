#!/bin/bash
# REEF3D coupled run (B 92 62): the models run side by side in one MPI job, app numbers in this order: regional local cfd
set -e
cd "$(dirname "$0")"
(cd regional && DiveMESH > divemesh.log)
(cd local && DiveMESH > divemesh.log)
(cd cfd && DiveMESH > divemesh.log)
mpirun -np 1 -wdir regional REEF3D : -np 1 -wdir local REEF3D : -np 1 -wdir cfd REEF3D > run.log 2>&1
