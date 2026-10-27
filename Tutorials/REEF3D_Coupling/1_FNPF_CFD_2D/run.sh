#!/bin/bash
# REEF3D coupled run (B 92 62): the models run side by side in one MPI job, app numbers in this order: fnpf cfd
set -e
cd "$(dirname "$0")"
(cd fnpf && DiveMESH > divemesh.log)
(cd cfd && DiveMESH > divemesh.log)
mpirun -np 1 -wdir fnpf REEF3D : -np 1 -wdir cfd REEF3D > run.log 2>&1
