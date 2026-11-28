#!/usr/bin/env bash
# Architect: Hans Bihs
# run_slide_cases.sh <REEF3D binary> <DiveMESH binary> <run dir> [out dir]
# Runs the two-layer landslide cases 16-19 (SFLOW slide layer W 150) on both grids, with the single-layer
# cases 01 and 02 for the comparison, and evaluates them with slide_eval.py (figures and JSON in <out dir>);
# then 18 with the non-hydrostatic water and case 28 (NTHMP benchmark 4) against the measurements (bm4_eval.py);
# then case 29 (the slide layer under NHFLOW), case 32 (NTHMP benchmark 5, subaerial) and case 35 (dilatancy, W 164).
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

# 18 with the non-hydrostatic water (A 220 2, W 162), one grid
run $C/18_slide_submarine/cases/dambreak_nh nh_dambreak
run $C/18_slide_submarine/cases/slope35_nh nh_slope35

# 28: NTHMP benchmark 4, and the non-hydrostatic case on the fine grid
B4=$C/28_slide_nthmp_bm4/cases
for c in hydro si nh nh_w128; do run $B4/$c bm4_$c; done
fine $B4/nh $RUN/_c/bm4_nh_fine; run $RUN/_c/bm4_nh_fine bm4_nh_fine

# 29: under NHFLOW
N4=$C/29_slide_nhflow/cases
for c in dambreak_h dambreak_nh bm4_h bm4_nh rest_slope dambreak_h_si dambreak_nh_si bm4_nh_si rest_slope_si bm4_nh_w163; do run $N4/$c nhf_$c; done

# 32: NTHMP benchmark 5
B5=$C/32_slide_nthmp_bm5/cases
for c in case1_sflow_h case1_sflow_nh case1_nhflow_nh case2_sflow_h case2_sflow_nh case2_nhflow_nh case1_nhflow_nh_si case2_nhflow_nh_si case1_nhflow_nh_w163 case2_nhflow_nh_w163; do run $B5/$c bm5_$c; done

# 35: dilatancy of the slide layer (W 164), a uniform submerged layer on 25 deg against the ODE of the closure
mkdir -p "$RUN/dil"
for c in loose dense coulomb dense_nhflow; do run $C/35_slide_dilatancy/cases/$c dil/$c; done

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
mkdir -p "$OUT/nh" "$OUT/28" "$OUT/28_grid"
$E dambreak  $RUN/nh_dambreak  $OUT/nh > $OUT/nh/18_dambreak.json
$E slope35   $RUN/nh_slope35   $OUT/nh > $OUT/nh/18_slope35.json
M=$C/28_slide_nthmp_bm4/test17_48-68s.txt
python3 $HERE/bm4_eval.py $M $OUT/28 hydrostatic=$RUN/bm4_hydro semi-implicit=$RUN/bm4_si "NH, W 162 0.02"=$RUN/bm4_nh \
        "NH, W 162 0.02, W 128 1"=$RUN/bm4_nh_w128 > $OUT/28/28_bm4.json
python3 $HERE/bm4_eval.py $M $OUT/28_grid "NH dx 5 mm"=$RUN/bm4_nh "NH dx 2.5 mm"=$RUN/bm4_nh_fine > $OUT/28_grid/28_bm4.json
mkdir -p "$OUT/29"
$E dambreak    $RUN/nhf_dambreak_h  $OUT/29 > $OUT/29/dambreak_h.json;  mv $OUT/29/18_dambreak.png $OUT/29/29_dambreak_h.png
$E dambreak    $RUN/nhf_dambreak_nh $OUT/29 > $OUT/29/dambreak_nh.json; mv $OUT/29/18_dambreak.png $OUT/29/29_dambreak_nh.png
$E rest_nhflow $RUN/nhf_rest_slope  $OUT/29 > $OUT/29/rest_slope.json
$E dambreak    $RUN/nhf_dambreak_h_si  $OUT/29 > $OUT/29/dambreak_h_si.json;  mv $OUT/29/18_dambreak.png $OUT/29/29_dambreak_h_si.png
$E dambreak    $RUN/nhf_dambreak_nh_si $OUT/29 > $OUT/29/dambreak_nh_si.json; mv $OUT/29/18_dambreak.png $OUT/29/29_dambreak_nh_si.png
$E rest_nhflow $RUN/nhf_rest_slope_si  $OUT/29 > $OUT/29/rest_slope_si.json
python3 $HERE/bm4_eval.py $M $OUT/29 "NHFLOW hydrostatic"=$RUN/nhf_bm4_h "NHFLOW NH, 10 layers, W 162 0.02"=$RUN/nhf_bm4_nh "NHFLOW NH, semi-implicit"=$RUN/nhf_bm4_nh_si "NHFLOW NH, W 163 1"=$RUN/nhf_bm4_nh_w163 \
        "SFLOW NH, W 162 0.02 (case 28)"=$RUN/bm4_nh > $OUT/29/29_bm4.json
mv $OUT/29/28_bm4_gauges.png $OUT/29/29_bm4_gauges.png
mkdir -p "$OUT/32"
for n in 1 2; do
    h=0.148; f=bm5_case1_d1.5mm.txt; if [ $n = 2 ]; then h=0.15; f=bm5_case2_d10mm.txt; fi
    python3 $HERE/bm4_eval.py --bench 5 --h $h --align crest --name 32_bm5_case$n $C/32_slide_nthmp_bm5/$f $OUT/32 \
            "SFLOW hydrostatic"=$RUN/bm5_case${n}_sflow_h "SFLOW NH, W 162 0.1"=$RUN/bm5_case${n}_sflow_nh \
            "NHFLOW NH, 10 layers"=$RUN/bm5_case${n}_nhflow_nh "NHFLOW NH, semi-implicit"=$RUN/bm5_case${n}_nhflow_nh_si "NHFLOW NH, W 163 1"=$RUN/bm5_case${n}_nhflow_nh_w163 \
            > $OUT/32/32_bm5_case$n.json
done
mkdir -p "$OUT/35"
python3 $HERE/dilatancy_eval.py $RUN/dil $OUT/35 > /dev/null
echo "results in $OUT"
