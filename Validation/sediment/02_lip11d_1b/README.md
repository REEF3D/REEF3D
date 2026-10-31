# Wave-driven beach profile: LIP 11D test 1B (Delta Flume, 1993)

This case validates the waves in the sediment transport (S 110, S 11 6–9, S 120) with REEF3D::SEASTATE in surfbeat mode, coupled with SFLOW and NHFLOW. It tests whether these models reproduce the measured offshore migration of a breaker bar under erosive waves.

## Experiment

Delta Flume, 230 m × 5 m × 7 m (Roelvink & Reniers 1995; Arcilla et al. 1994; van Rijn 1999, 2007).

| | test 1B |
|---|---|
| waves | JONSWAP, Hm0 = 1.4 m, Tp = 5 s |
| still water | 4.1 m above the flume floor |
| sand | d50 = 0.22 mm (0.20 mm in van Rijn 1999), d90 = 0.33 mm |
| duration | 18 h |

**Reference profiles** (`refdata/`), from the XBeach skillbed input files (`svn.oss.deltares.nl/repos/xbeach/skillbed/input/DeltaflumeLIP11D/`, `h_1B.dep`, `h_1C.dep`):
- `LIP11D_1B_initial_profile.csv`: the profile at the start of 1B.
- `LIP11D_1C_initial_profile_approx_1B_end.csv`: the profile at the start of 1C, which is the end of 1B because each LIP test started from the previous end profile.
- Both files have x from the wave board and z above the flume floor.
- In 1B, the values beyond x = 200 m are a linear extrapolation in the model file.

The measured Hrms, undertow and intermediate profiles are in the binary files `skillbed/data/DeltaflumeLIP11D/1B.mat` and `1C.mat`. These could not be downloaded here (`svn export https://svn.oss.deltares.nl/repos/xbeach/skillbed/data/DeltaflumeLIP11D/`).

## Set-up (`../tools/mk_lip.py`)

**Domain.** x = 0 … 210 m, Δx = 1 m, profile from `geo.dat`.

**Waves (SEASTATE surfbeat).**
- Surfbeat mode (A 770 1), boundary JONSWAP 1.4 m / 5 s;
- Roelvink breaking with γ = 0.40 (A 740 2, A 741);
- roller (A 748 1);
- radiation stress (A 751 1); bound long waves at x- (A 774 1).

**Sediment.**
- Soulsby–van Rijn (S 11 6) with the suspended part in the depth-averaged suspended load (S 120 2);
- Ruessink/Abreu wave shape (S 116 1), asymmetry S 117 0.1 and breaking stirring (S 118 1);
- Eulerian current (S 119 0), sandslide (S 90 1);
- morphological factor 10 from 600 s, i.e. 600 s + 6480 s.

**Flow models.**
- SFLOW: Manning bed friction with B 50 = 0.01 m.
- NHFLOW: explicit, 5 σ-layers, k-ε, A 519 1 with B 50 (S 28 2).

## Results (`results/fig9_lip1b.png`)

Brier skill score of the profile after 18 h for x = 100 … 190 m, with the initial profile as the baseline:

| run | BSS | RMS |
|---|---|---|
| SFLOW, γ 0.40 | 0.02 | 8.9 cm |
| NHFLOW, γ 0.40 | −0.03 | 9.1 cm |
| SFLOW, γ 0.45 (default case) | −0.09 | 9.3 cm |
| SFLOW, γ 0.35 / 0.55 | 0.01 / −0.59 | |
| SFLOW, no asymmetry S 117 0 | −0.25 | |
| SFLOW, no roller / no stirring S 118 0 / no suspended model S 120 0 | −0.07 / −0.05 / −0.13 | |

**Not validated yet.** The models have no skill on LIP 1B.

What the measurements show:
- The bar moves about 8 m offshore.
- Up to 0.27 m of sand is deposited at x = 128 – 138 m.
- The trough at x = 140 – 156 m is eroded by up to 0.15 m.

What the models do:
- Bed changes reach only about 0.1 m.
- The bar crest is eroded at x = 140 – 145 m.
- Sand is deposited shoreward of the bar, at x = 150 – 158 m, instead of offshore.

So the net transport over the bar is too weak and has the wrong distribution. In a time average the undertow is offshore, as expected: the host mass flux averages to about 0 in the closed flume, and M/h is about 0.3 m/s at the bar.

Next steps:
1. Compare Hrms and the undertow with the measured data in `1B.mat`, to separate the hydrodynamics from the transport.
2. Check how the instantaneous surfbeat velocities (long waves) combine with the stirring in the wave-averaged formulas.
3. Run the van Rijn 2007 bed and suspended load (S 11 7 + S 120 1) with a calibrated wave-related suspended transport. XBeach reaches positive skill with its sediment transport parameter facua and the roller.

## Running

```bash
../tools/mk_lip.py <dir> sflow "A 741=1.0 0.40"                          # = cases/sflow_1b
../tools/mk_lip.py <dir> nhflow "A 741=1.0 0.40" "A 534=del" "A 535=del" "A 536=del"   # = cases/nhflow_1b
../tools/run_case.sh <REEF3D> <DiveMESH> cases/sflow_1b <run dir>
python3 ../tools/sed_eval.py lip <run dir>                              # BSS, RMS for x = 100 .. 190 m
```

Run time on one core: SFLOW about 3 min, NHFLOW about 40 min.
