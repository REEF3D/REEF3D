#!/bin/bash
# REEF3D with rebalancing (G 46): runs the given command, e.g.
#   tools/reef3d_rebalance.sh mpirun -np 64 ./REEF3D
# and repeats it as long as REEF3D stopped to rebalance: checkpoint, load map and the marker
# REEF3D_Checkpoint/rebalance (exit code 3; some mpirun versions pass 0 on, so the marker decides).
# The next start reads the checkpoint (I 45) and decomposes the global grid of DIVEMesh (M 50 1, 2)
# with the weighted partition of the load map; it removes the marker.
# REEF3D_REBALANCE_MAX: the most restarts (default 20).
max=${REEF3D_REBALANCE_MAX:-20}
n=0
while true; do
    "$@"
    rc=$?
    if [ ! -f REEF3D_Checkpoint/rebalance ]; then
        exit $rc
    fi
    n=$((n+1))
    if [ $n -gt $max ]; then
        echo "reef3d_rebalance.sh: more than $max rebalance restarts -- stopped"
        exit 3
    fi
    echo "reef3d_rebalance.sh: rebalance restart $n"
done
