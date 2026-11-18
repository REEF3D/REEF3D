#!/usr/bin/env bash
# Architect: Hans Bihs
# run_slide_cases.sh <REEF3D binary> <DiveMESH binary> <run dir> [out dir]
# Runs the two-layer landslide cases 16-19 (SFLOW slide layer W 150) on both grids, with the single-layer
# cases 01 and 02 for the comparison, and evaluates them with slide_eval.py (figures and JSON in <out dir>).
set -euo pipefail
BIN="$1"; DM="$2"; RUN="$3"; OUT="${4:-$RUN/results}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; V="$(cd "$HERE/.." && pwd)"
MPI="${REEF3D_MPIRUN:-mpirun} -np 2"
mkdir -p "$RUN" "$OUT"

# a copy of a case with the cell size of the grid (B 1) and the width of a 1D grid (B 10 ye) halved
fine() {
    local src="$1" dst="$2"
    rm -rf "$dst"; mkdir -p "$dst"; cp "$src"/control.txt "$src"/ctrl.txt "$dst"/
    python3 - "$dst" <<'EOF'
import sys, re
d = sys.argv[1]
c = open(d + "/control.txt").read().split("\n")
dx = None
for k, l in enumerate(c):
    s = l.split()
    if len(s) == 3 and s[0] == "B" and s[1] == "1":
        dx = float(s[2]); c[k] = "B 1 %g" % (dx/2)
for k, l in enumerate(c):
    s = l.split()
    if len(s) == 8 and s[0] == "B" and s[1] == "10" and abs(float(s[5]) - float(s[4]) - dx) < 1e-12:
        s[5] = "%g" % (float(s[4]) + dx/2); c[k] = " ".join(s)
open(d + "/control.txt", "w").write("\n".join(c))
# P 51 gauges of a 1D grid: y in the middle of the new cell
t = open(d + "/ctrl.txt").read().split("\n")
for k, l in enumerate(t):
    s = l.split()
    if len(s) == 4 and s[0] == "P" and s[1] == "51" and dx is not None and float(s[3]) < dx:
        s[3] = "%g" % (dx/4); t[k] = " ".join(s)
open(d + "/ctrl.txt", "w").write("\n".join(t))
EOF
}

run() {   # run <case dir> <run name> [mpirun prefix]
    "$HERE"/run_case.sh "$BIN" "$DM" "$1" "$RUN/$2" "${3:-}" > /dev/null
    echo "$2: $(tail -n 1 "$RUN/$2/reef3d.log" | tr -d '\n') $(grep -m1 'total time' "$RUN/$2/reef3d.log" >/dev/null && grep 'total time' "$RUN/$2/reef3d.log" | tail -n 1 | awk '{print "("$3" s)"}')"
}

C=$V
# 18: the slide alone, and the single-layer runs it must reproduce
run $C/01_coulomb_dambreak_mangeney/cases/mangeney single_mangeney
run $C/02_viscoplastic_arrest/cases/hb_core single_arrest
run $C/16_slide_alone/cases/mangeney s_mangeney
run $C/16_slide_alone/cases/arrest s_arrest
fine $C/16_slide_alone/cases/mangeney $RUN/_c/mangeney_fine; run $RUN/_c/mangeney_fine s_mangeney_fine
fine $C/16_slide_alone/cases/arrest $RUN/_c/arrest_fine;     run $RUN/_c/arrest_fine s_arrest_fine

# 17: moving bed, 2D, two ranks
run $C/17_slide_moving_bed/cases/rest rest "$MPI"
run $C/17_slide_moving_bed/cases/subaerial subaerial "$MPI"
fine $C/17_slide_moving_bed/cases/rest $RUN/_c/rest_fine;           run $RUN/_c/rest_fine rest_fine "$MPI"
fine $C/17_slide_moving_bed/cases/subaerial $RUN/_c/subaerial_fine; run $RUN/_c/subaerial_fine subaerial_fine "$MPI"

# 18: submarine
run $C/18_slide_submarine/cases/dambreak dambreak
run $C/18_slide_submarine/cases/slope35 slope35
fine $C/18_slide_submarine/cases/dambreak $RUN/_c/dambreak_fine; run $RUN/_c/dambreak_fine dambreak_fine
fine $C/18_slide_submarine/cases/slope35 $RUN/_c/slope35_fine;   run $RUN/_c/slope35_fine slope35_fine

# 19: the same with the semi-implicit water (A 210 11)
S=$C/19_slide_semi_implicit/cases
run $S/rest si_rest "$MPI"
run $S/subaerial si_subaerial "$MPI"
run $S/dambreak si_dambreak
run $S/slope35 si_slope35
fine $S/rest $RUN/_c/si_rest_fine;           run $RUN/_c/si_rest_fine si_rest_fine "$MPI"
fine $S/subaerial $RUN/_c/si_subaerial_fine; run $RUN/_c/si_subaerial_fine si_subaerial_fine "$MPI"
fine $S/dambreak $RUN/_c/si_dambreak_fine;   run $RUN/_c/si_dambreak_fine si_dambreak_fine
fine $S/slope35 $RUN/_c/si_slope35_fine;     run $RUN/_c/si_slope35_fine si_slope35_fine

E="python3 $HERE/slide_eval.py"
for g in "" _fine; do
    mkdir -p "$OUT/grid$g"
    $E mangeney  $RUN/s_mangeney$g $RUN/single_mangeney $OUT/grid$g > $OUT/grid$g/16_mangeney.json
    $E arrest    $RUN/s_arrest$g   $RUN/single_arrest   $OUT/grid$g > $OUT/grid$g/16_arrest.json
    $E rest      $RUN/rest$g       $OUT/grid$g > $OUT/grid$g/17_rest.json
    $E subaerial $RUN/subaerial$g  $OUT/grid$g > $OUT/grid$g/17_subaerial.json
    $E dambreak  $RUN/dambreak$g   $OUT/grid$g > $OUT/grid$g/18_dambreak.json
    $E slope35   $RUN/slope35$g    $OUT/grid$g > $OUT/grid$g/18_slope35.json
    mkdir -p "$OUT/grid$g/19_si"
    $E rest      $RUN/si_rest$g       $OUT/grid$g/19_si > $OUT/grid$g/19_si/17_rest.json
    $E subaerial $RUN/si_subaerial$g  $OUT/grid$g/19_si > $OUT/grid$g/19_si/17_subaerial.json
    $E dambreak  $RUN/si_dambreak$g   $OUT/grid$g/19_si > $OUT/grid$g/19_si/18_dambreak.json
    $E slope35   $RUN/si_slope35$g    $OUT/grid$g/19_si > $OUT/grid$g/19_si/18_slope35.json
done
echo "results in $OUT"
