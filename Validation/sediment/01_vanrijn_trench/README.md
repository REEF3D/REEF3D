# Trench migration in a current: van Rijn (1986) flume tests

Validation of REEF3D's suspended-load models against measurements:
- NHFLOW suspended load (S 12 1);
- SFLOW depth-averaged suspended load (S 120 1);
- the new sediment supply at the inflow (S 129, S 130, S 124).

## Experiment

Delft Hydraulics flume tests, 1980 (van Rijn & Tan 1985, Rijkswaterstaat Communication 41; van Rijn 1986, J. Hydraul. Eng. 112(6), Figs. 16–17; van Rijn 1987, 1993).

| | value |
|---|---|
| flume | 30 m × 0.5 m |
| approach flow | h0 = 0.39 m, U0 = 0.51 m/s (q = 0.199 m²/s) |
| sand | d50 = 0.16 mm, d90 = 0.20 mm; suspended sand 0.13 mm, ws = 0.013 m/s |
| supply at the inlet (equilibrium) | 0.04 kg/s/m: 0.03 suspended + 0.01 bed load |
| bed roughness | ks = 0.025 m (small bedforms) |
| duration | bed measured at t = 0 and 15 h |

Three trenches were tested, with side slopes 1:3 (test 1, depth 0.175 m), 1:7 (test 2, 0.170 m) and 1:10 (test 3, 0.150 m). Over the 15 hours, the trench fills and migrates downstream, and the bed downstream of it erodes because the trench traps sediment.

`refdata/van_rijn_1986_fig16_initial.csv` and `refdata/van_rijn_1986_fig16_measured.csv` hold the initial trench and the measured bed after 15 h. Coordinates are x from station 1 (1 m upstream of the lip in test 1) and z below the flat bed.

- **Source:** digitised by the ANUGA contributors from a 600 dpi scan of van Rijn (1986) Fig. 16 (anuga-community/anuga_core, `validation_tests/sediment/van_rijn_trench`, PR #393).
- **Accuracy:** about ±0.1 m in x and ±0.003 m in z.
- **Coverage:** only the downstream slope and the reach beyond the trench were measured.

## Set-up (`../tools/mk_trench.py`)

**Domain.** x = −20 … 20 m: a 20 m approach, the trench at 0 … 6 m, and 14 m beyond. Δx = 0.1 m, flat bed at z = 0.30 m, trench from `geo.dat` (DIVEMesh G 10).

**Flow.**
- Discharge inflow: B 60 1, W 10 = q·0.5 m.
- Water level 0.69 m at the outflow (F 62); the outflow bed is fixed (S 50 2).
- Morphological factor S 35 10: 300 s spin-up + 5400 s = 15 h.
- Supply at the inlet (S 129 1):
  - bed load S 130 = 0.01/2650 m²/s;
  - concentration S 124 = the suspended supply divided by the discharge (SFLOW), or by the flux β·q with β = 0.73, the ratio ∫uc dz/(q C̄) of the approach flow (NHFLOW).

| | NHFLOW | SFLOW |
|---|---|---|
| hydrodynamics | hydrostatic, semi-implicit (A 534 1, θ = 1), 20 σ-layers, k-ε (A 560 1, A 567 1), bed friction A 519 1 with the hydraulic roughness B 50 = 0.025 m (S 28 2, new) | hydrostatic (A 220 0), semi-implicit (A 226 1, Δt ≤ 0.1 s), Manning from B 50 = 0.025 m |
| bed load | van Rijn (S 11 1) | van Rijn (S 11 1) |
| suspended load | S 12 1, van Rijn reference concentration | S 120 1, van Rijn (2007) c_a + Rouse |
| settling velocity | S 23 0.013 (S 25 0) | S 23 0.013 |
| grain roughness S 21 | **30** (k_s' = 4.8 mm) | **8** (k_s' = 1.3 mm) |

**Calibration.** As in van Rijn's and ANUGA's studies, the grain roughness S 21 was set so that the approach flow carries the measured supply. It was tuned on the flat flume without the trench (S 35 0):

| | bed load [kg/s/m] | suspended [kg/s/m] |
|---|---|---|
| measured supply | 0.010 | 0.030 |
| NHFLOW, S 21 30 | 0.008 | 0.034 |
| NHFLOW, S 21 3.75 (k_s' = 3 d90) | 0.0005 | 0.003 |
| SFLOW, S 21 8 | 0.010 – 0.012 | 0.025 – 0.030 |

In NHFLOW the sediment's bed shear stress comes from the bottom-layer velocity, with the log law of the grain roughness. The bottom velocity, however, belongs to the profile of the bedform roughness, 0.025 m. Without calibration, the grain-roughness log law at the bottom cell gives about 0.2 Pa instead of the 0.5 Pa of van Rijn's depth-averaged skin friction, so k_s' has to absorb the difference. A skin-friction efficiency factor μ = (C/C')² would remove this calibration.

## Results (15 h, `results/fig8_trench_validation.png`)

RMS of the bed level at the measured points, and the deepest point after 15 h:

| test | measured deepest | NHFLOW RMS | NHFLOW deepest | SFLOW RMS | SFLOW deepest |
|---|---|---|---|---|---|
| 3 (1:10) | 0.079 m at 6.57 m | **2.5 cm** | 0.088 m at 7.4 m | 5.0 cm | 0.075 m at 10.9 m |
| 2 (1:7) | 0.082 m at 5.96 m | **1.9 cm** | 0.090 m at 6.5 m | 5.2 cm | 0.076 m at 10.7 m |
| 1 (1:3) | 0.100 m at 6.11 m | **3.3 cm** | 0.111 m at 7.2 m | 6.1 cm | 0.092 m at 10.4 m |

For comparison:
- van Rijn's own 2DV model (SUTRENCH): about 1 cm;
- ANUGA (depth-averaged, two-layer closure, tuned): 1.7–2.5 cm, and 2.9–3.9 cm with instantaneous exchange;
- SCHISM (2017): about 2 cm.

**NHFLOW.**
- The trench fills from upstream and is about half its depth at 15 h. The deepest point is 0.4–1 m too far downstream and 1 cm too deep.
- Downstream of the trench the bed erodes as measured.
- The deposit on the downstream slope arrives as a steep front; in the flume it was spread over the slope.
- The approach reach aggrades by 2–4 cm because the calibrated capacity is slightly below the supply.

**SFLOW.**
- The trench depth after 15 h is right, but it has migrated about 4 m too far: the upstream slope and the floor fill completely, and the remaining trench sits at 10–11 m.
- This matches the earlier finding that the equilibrium-profile factor F adapts about twice too fast: the suspended sediment deposits as soon as the flow slows, instead of staying in the water column over the trench. The depth-averaged model needs a separate adaptation time (Galappatti & Vreugdenhil) to keep the sediment high in the column; the ANUGA study shows the same, and its two-layer closure was the fix.
- S 121 is not the knob for this: it sets F in the exchange and in C_eq = c_a/F together. S 121 1 and 2 give 16 and 12 cm.

**Bug found by this case.** The sediment model had zero face velocities at the domain edges, in SFLOW and in NHFLOW. Bed load could neither enter nor leave through an in- or outflow edge, and in SFLOW the first cell carried almost no bed load. With the closed inflow face, the supply S 129 had no effect. Before the fix, SFLOW gave 7–8 cm with a deposit 4 cm above the original bed behind the front.

## Running

```bash
../tools/mk_trench.py <dir> nhflow 3 "S 21=30" "S 124=7.8e-5"     # writes control.txt, ctrl.txt, geo.dat (= cases/nhflow_test3)
../tools/mk_trench.py <dir> sflow 3 "S 21=8"                      # = cases/sflow_test3
../tools/run_case.sh <REEF3D> <DiveMESH> cases/nhflow_test3 <run dir>
python3 ../tools/sed_eval.py trench <run dir> 3                    # RMS, deepest point
python3 ../tools/flux.py <run dir> -15 0 15                        # transport along the flume (NHFLOW needs P 26 1)
```

Run time on one core: NHFLOW about 40 min, SFLOW about 15 min per test.

Builds: hans_dev 3f7c738 + `261031_REEF3D_sediment_supply_validation-on-3f7c738.patch`.
