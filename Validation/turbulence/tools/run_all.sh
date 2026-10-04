#!/usr/bin/env bash
# Architect: Hans Bihs
# run_all.sh <REEF3D binary> <DiveMESH binary> <cases dir> <output dir>   - runs every case in <cases dir> (single rank)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$1"; DM="$2"; CASES="$3"; OUT="$4"; mkdir -p "$OUT"
for c in "$CASES"/*/; do
  n=$(basename "$c"); echo "== $n"
  "$HERE/run_case.sh" "$BIN" "$DM" "$c" "$OUT/$n"
done
