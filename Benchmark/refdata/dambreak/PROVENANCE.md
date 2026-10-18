# Dam-break reference data -- provenance

Fetched 2026-10-03. All numbers copied verbatim from the fetched sources (each numeric block fetched at least twice
via different URLs and compared); no rescaling. Conversions appear only in file headers.

| File | Rows | Content | Source | Confidence these are original experimental values |
|---|---|---|---|---|
| `martin_moyce_1952_surge_n2_2_a1.125in_pysph.dat` | 10 | M&M surge front, n^2=2, a=1.125 in; T=t*sqrt(2g/a), Z=z/a | PySPH `pysph/examples/db_exp_data.py` (`mm_data_1`) | Medium: experimental, but **digitised from M&M Fig. 3** by PySPH ("extracted from plots") |
| `martin_moyce_1952_surge_n2_2_a2.25in_pysph.dat` | 15 | Same, a=2.25 in | PySPH `mm_data_2` | Medium (digitised, Fig. 3) |
| `martin_moyce_1952_surge_n2_2_lethe.dat` | 14 | M&M surge front; tau=t*sqrt(2g/L), delta=x/L | Lethe `examples/multiphysics/dam-break/dam-break-2d.py` (`x_exp`,`y_exp`) | Medium: delta values sit on regular steps (1.11, 1.23, 1.44 ...), suggesting fixed-position readings, but Lethe does not say which a or whether tabulated/digitised |
| `koshizuka_oka_1996_surge_pysph.dat` | 9 | K&O experiment, T=t*sqrt(2g/L), Z=z/L | PySPH `ko_data` | Medium (digitised) |
| `koshizuka_oka_1996_surge_dualsphysics.dat` | 15 | K&O tip position, columns `time Xmax` | DualSPHysics `examples/main/01_DamBreak/EXP_X-DamTipPosition_Koshizula&Oka1996.txt` | Medium-low: units not stated; appears rescaled to a 1 m x 2 m column (dimensional s, m) |
| `kleefsman_2005_P1_pysph.dat` | 47 | MARIN P1 pressure; T=t*sqrt(g/H), p/(rho g H), H=0.55 m | PySPH `kleefsman_exp_data_p1` | Medium-low: digitised, coarse; PySPH comment says sqrt(2g/H) but its own plotting uses sqrt(g/H) (adopted here) |
| `kleefsman_2005_P3_pysph.dat` | 49 | MARIN P3 pressure, same scaling | PySPH `kleefsman_exp_data_p3` | Medium-low (same caveats) |
| `kleefsman_2005_geometry.txt` | -- | Tank, box, sensor coordinates from 5 sources, with discrepancies noted | Veldman ComFLOW page, LS-DYNA example, SPHERIC PDF, paper text, Whiterose eprint | High for tank/box/H2/H4/P1/P3/P5/P7 (ComFLOW page, co-author). P1/P3 z conflict (0.025/0.099 vs 0.021/0.101) unresolved |

## URLs
- PySPH: https://raw.githubusercontent.com/pypr/pysph/master/pysph/examples/db_exp_data.py ;
  https://github.com/pypr/pysph/blob/master/pysph/examples/db_exp_data.py?plain=1 ;
  https://github.com/pypr/pysph/blob/main/pysph/examples/dam_break_3d.py?plain=1
- Lethe: https://raw.githubusercontent.com/chaos-polymtl/lethe/master/examples/multiphysics/dam-break/dam-break-2d.py
- DualSPHysics: https://github.com/DualSPHysics/DualSPHysics/blob/master/examples/main/01_DamBreak/EXP_X-DamTipPosition_Koshizula%26Oka1996.txt ;
  https://github.com/DualSPHysics/DualSPHysics/blob/master/examples/main/01_DamBreak/CaseDambreakVal2D_Def.xml
- ComFLOW (Veldman): https://www.math.rug.nl/~veldman/comflow/dambreak.html
- SPHERIC Test 2: https://www.spheric-sph.org/tests/test-02 (zip with `test_case_2_exp_data.xls`, unreachable here)
- LS-DYNA: https://www.dynaexamples.com/icfd/intermediate-examples/dam4

## Nondimensionalisation notes
- Martin & Moyce: with a = initial column width and n^2 = height/width (n^2 = 2 -> height 2a), T = t*sqrt(2g/a), Z = z/a
  (as stated by PySPH and used by Lethe). Not confirmed against the original paper (Royal Society/JSTOR blocked).

## Not found
- Martin & Moyce residual column height vs time: no machine-readable source found.
- Original tabulated M&M values (paper inaccessible).
- Kleefsman H1-H4 water-height time series and P2, P4-P8 pressure series (only in SPHERIC zip/xls; ComFLOW data
  links are dead).
- Unambiguous H1 and H3 coordinates.
