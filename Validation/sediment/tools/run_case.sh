#!/bin/bash
# run_case.sh <REEF3D> <DiveMESH> <case dir> <run dir> [mpirun prefix]: copy the inputs, run DiveMESH and REEF3D
set -e
mkdir -p "$4"; cp "$3"/control.txt "$3"/ctrl.txt "$4"/; [ -f "$3"/geo.dat ] && cp "$3"/geo.dat "$4"/
cd "$4"; "$2" > divemesh.log 2>&1; ${5:-} "$1" > reef3d.log 2>&1
