#!/usr/bin/env bash
# Architect: Hans Bihs
# run_case.sh <REEF3D binary> <DiveMESH binary> <case dir> <run dir> [mpirun prefix]
# Copies control.txt/ctrl.txt into <run dir>, runs DiveMESH and REEF3D there.
# The regression dump (src/regression_dump.cpp) is written to <run dir>/reg; the CFD analysis reads it.
set -euo pipefail
BIN="$1"; DM="$2"; CASE="$(cd "$3" && pwd)"; RUN="$4"; MPI="${5:-}"
rm -rf "$RUN"; mkdir -p "$RUN/reg"; RUN="$(cd "$RUN" && pwd)"
cp "$CASE"/control.txt "$CASE"/ctrl.txt "$RUN"/
cd "$RUN"
"$DM" > divemesh.log 2>&1
REEF3D_REGRESSION_DIR="$RUN/reg" $MPI "$BIN" > reef3d.log 2>&1 || echo "REEF3D exit $?" >> reef3d.log
tail -n 2 reef3d.log
