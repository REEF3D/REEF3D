#!/bin/bash
# Run one level of the benchmark suite and make the plots.
#
#   tools/run_suite.sh nightly            # all nightly cases
#   tools/run_suite.sh release 'cfd_*'    # release level, CFD cases only
#
# Environment: REEF3D (binary), DIVEMESH (binary), REEF3D_MPIRUN (launcher), OUT (output root),
#              MAX_NP (cap on ranks, optional)
set -e
LEVEL=${1:-nightly}
shift || true
HERE=$(cd "$(dirname "$0")/.." && pwd)
REEF3D=${REEF3D:-$HOME/Codelite/REEF3D/bin/REEF3D}
DIVEMESH=${DIVEMESH:-$HOME/Codelite/DIVEMesh/bin/DiveMESH}
OUT=${OUT:-$HOME/reef3d_bench/$LEVEL-$(date +%y%m%d-%H%M)}
EXTRA=()
[ -n "$MAX_NP" ] && EXTRA+=(--max-np "$MAX_NP")
[ $# -gt 0 ] && EXTRA+=(--cases "$@")
cd "$HERE"
./benchmark.py run --level "$LEVEL" --reef3d "$REEF3D" --divemesh "$DIVEMESH" --out "$OUT" "${EXTRA[@]}" --check || rc=$?
./benchmark.py plot "$OUT" > /dev/null 2>&1 || true
echo "report: $OUT/benchmark.md"
exit ${rc:-0}
