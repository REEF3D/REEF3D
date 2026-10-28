# Cavitation (W 60): validation

Homogeneous liquid-vapour mixture with Schnerr-Sauer mass transfer, `src/cavitation.{h,cpp}`
(formulation in the header), mixture face density `src/density_cav.{h,cpp}`, the phase-change
divergence implicit in the pressure correction (`pjm_corr` -> `cavitation::projection`).

Keywords

| | | default |
|---|---|---|
| `W 60` | 1: Schnerr-Sauer | 0 |
| `W 61` | vapour pressure p_v [Pa], in the pressure scale of the solver | 0 |
| `W 62` | vapour density [kg/m³] | 0.0173 |
| `W 63` | vapour kinematic viscosity [m²/s] | 5.6e-4 |
| `W 64` | nuclei per volume of liquid n0 [1/m³] | 1.6e13 |
| `W 65` | nuclei diameter [m] | 2e-6 |
| `W 66` | evaporation and condensation coefficients C_v C_c | 1 1 |
| `W 67` | initial vapour box xs xe ys ye zs ze | |
| `W 68` | initial vapour sphere x y z r | |
| `W 69` | cavitation number control sigma U x y z tau: from t = tau on p_v = p_ref - sigma rho U²/2, p_ref the pressure at (x, y, z) averaged over tau; before that p_v = W 61 | off |
| `W 70` | Reboud correction of the eddy viscosity (RANS): mu_t with f = rho_v + (1 - alpha)^n (rho_l - rho_v) instead of rho, n (10) | 0 (off) |

Setup notes
- The pressure level must be fixed by a Dirichlet boundary: the outlet side a wall `C 14 21`
  with an outlet patch over it, `B 441 id 4 xs xe ys ye zs ze` and `B 412 id 0.0` (pressure
  correction = patch pressure − cell pressure), with `B 60 1`. p_v is relative to it:
  W 61 = -sigma rho U²/2 for the cavitation number sigma = (p_inf - p_v)/(rho U²/2), p_inf = 0.
  Alternative without patches (needed with mesh refinement, which does not take B 440-442): the
  plain outflow `C 14 2`, `B 60 1`, `B 77 10`; with B 77 10 the outflow correction is
  outflow pressure − cell pressure (Dirichlet, `poisson_pcorr`), case 01 then collapses at
  0.1835 s (+0.6 %). Not suitable: `B 60 0` (no outlet pressure at all) and the gauge `B 30 1`
  (shifting the level after each projection against the implicit phase-change term: the vapour
  vanishes within a few ms).
- With `W 69` the run holds the wanted sigma whatever the pressure level: the drag of the body
  lowers the upstream pressure against the outlet (head form D/40: -15.5 Pa with the outlet patch,
  -3.5 to -5 Pa with the outflow B 77 10 of the refined runs), at fixed W 61 sigma was 0.39 and 0.42
  for the same nominal 0.4. The control must be slow and start after the start-up: with tau 1 s
  from t = 0 (first version) the start-up pressure (+950 Pa upstream at D/20 + AMR) went into p_v,
  the growing cavity pushed the water out and kept the upstream pressure up, 14 m³ of vapour and a
  damped oscillation of 7 s. Now from t = tau on, tau 5 s, W 61 a good estimate: D/10 sigma 0.25,
  0.282 at 5 s, 0.255 at 10 s. p_v is written as the last column of the cavitation log.
- Gravity off (`W 22 0`) unless the hydrostatic pressure is meant to be part of p - p_v.
- Single phase (`F 30 0`, `F 80 0`); a fluid at rest needs a fixed time step (`N 48 0`, `N 49 dt`).
- Output: `alpha_v` in the VTU files (`P 24 1` adds the mixture density), and
  `REEF3D_CFD_Cavitation/REEF3D-CFD-Cavitation.dat`: time, vapour volume, equivalent radius
  (3V/4π)^(1/3), vapour production rate, minimum pressure.
- Restart: the checkpoints (`P 48 n`, restart `I 45 1`) hold alpha_v with every other field; the
  restart of case 01 from step 100 is bitwise identical to the continuous run, the .dat file is
  appended. The state files (`P 40`, `I 40`) do not hold alpha_v.
- Mesh refinement (`A 10 6`, `G 1 > 0`, `G 15` boxes): supported with one time step for all
  levels (`G 7 0`), the pressure level by `B 77 10` (above). One cavitation object per grid,
  alpha_v filled around the patches and restricted, the phase-change rows of every grid in the
  composite projection; alpha_v is not flux-matched at the patch edges. The patches are not in
  the .dat file (level 0 only, the restricted alpha_v).

## 01 Collapse of a vapour layer

`cases/01_slab_collapse`: channel 1 m × 0.02 m × 0.02 m (dx = 5 mm), closed at x = 0, outlet
patch at x = 1 m (p = 0); a vapour layer 0.2 m thick at the closed end with p_v = -10 kPa. The liquid column
between the layer and the outlet is a piston, ρ (L − h) ḧ = −(p_out − p_v), the 1D Rayleigh problem.
4 ranks, 800 steps, ~1 min.

    DiveMESH && mpirun -np 4 REEF3D
    python ../../scripts/slab_compare.py --plot

Result (`cases/01_slab_collapse/REEF3D-CFD-Cavitation.dat`): collapse at 0.1830 s against
0.1825 s of the piston equation (+0.3 %), the layer thickness within 1.2 mm over the whole
collapse (`slab_collapse.png`). The collapse is held by the implicit coupling: the pressure at the
interface stays at p_v and the column accelerates under the full p_out − p_v.

Short regression version: `Regression/cases/cfd_1d_cavitation_slab` (20 steps).

## 02 Hemispherical head form (Rouse & McNown 1948)

`cases/02_headform`: cylinder of diameter D = 1 m with a hemispherical nose at zero yaw,
Re = 1.36e5 (U = 1 m/s), sheet cavitation behind the shoulder. `make_case.py` writes the body
(`floating.stl`, closed, 29808 triangles), the inputs and a Betzy job per cavitation number into
`d<D/dx>/<name>/` (git-ignored, 7 MB per case):

    python3 make_case.py 40 128 15        # D/40, 3.0 M cells, 1 Betzy node, 15 s: nocav, sigma 0.5, 0.4, 0.3
    python3 make_case.py                  # D/20, 0.74 M cells, 8 ranks, 20 s (local)

Copy `REEF3D` and `DiveMESH` into each case directory, set the account in `job.slurm`, `sbatch`.
Domain [-4, 12] × [-4, 4]², inflow U at x = -4, outlet pressure patch p = 0 on the x+ wall;
p_v = -sigma rho U²/2; k-epsilon with wall functions at the direct-forcing body (`T 10 1`);
density ratio rho_l/rho_v = 1000 (`W 62 1`, `W 63 1e-5`) and CFL 0.2 (`N 47`): with the real
ratio 58000 the D/40 sigma 0.5 run diverged at t = 3.3 s (lateral velocity 567 m/s in the vapour,
k up to 4 m²/s²); at D/20 sigma 0.4 the two ratios give the same cavity (s/D 0.725-2.225), the
same vapour volume (0.43 m³) and Cp within 0.006, as reported for sheet cavitation in the
literature (Kunz et al. 2000, Senocak & Shyy 2002 use 1000).
VTU frames over the last 5 s (`P 35`), probe 1 upstream for the reference pressure.

    python ../../scripts/headform_analysis.py nocav sigma0.5 sigma0.4 sigma0.3 --exp ../../rouse_mcnown.csv --plot headform_cp.png

(from `d40/`; the csv is `cases/02_headform/rouse_mcnown.csv`). Cp is the wall pressure, extrapolated
from two bands 0.5-1.5 dx and 1.5-2.5 dx off the surface (at the nose the normal pressure gradient
ρu²/R makes Cp 0.05 D off the wall ~0.3 higher), averaged over time and azimuth; sigma from
probe 1; the cavity from alpha_v > 0.1 next to the wall.

Measured data: `rouse_mcnown.csv`, read from the scanned figure of the hemispherical head in the
report (K = sigma ≥ 0.80, 0.50, 0.40, 0.30 at the 15 piezometers, about ±0.02, see the header):
no cavitation Cp_min = -0.74 at s/d = 0.66 (on the nose, before the shoulder at 0.785); with
cavitation Cp = -K along the cavity, closure at s/d ≈ 1.0 (K 0.5), 1.3-1.5 (K 0.4), 1.8-2.1 (K 0.3)
with a small overshoot.

Local D/20 runs (12 s, window 7-12 s), no cavitation: the nose up to s/D = 0.5 agrees with the
measurement; the suction peak is too weak, laminar Cp_min = -0.45 with reverse flow behind the
shoulder to s/D ≈ 2, k-epsilon -0.49 at s/D = 0.63 with the separation to s/D ≈ 1.1. Ten cells
over the nose radius are not enough for the direct-forcing body; D/40 on Betzy is the comparison.

D/40 on Betzy (3.0 M cells, 128 ranks, 15 s, 2974 steps, rho_l/rho_v = 1000), sigma 0.5 case,
window 10-15 s: p_ref = -15.5 Pa, so the actual sigma = 0.469. Cavity from s/D = 0.71 (just before
the shoulder, as measured) with Cp = -0.469 = -sigma along it, alpha_v up to 0.88 at the wall,
closure (Cp = -0.1) at s/D ≈ 1.40 with an overshoot +0.125 at 1.5; vapour volume 0.12 m³
oscillating at about 0.2 Hz. Measured at K = 0.50: plateau -0.50 over s/d 0.66-0.92, closure
≈ 1.15, no overshoot. The cavity is ~0.25 D longer; part of it is the lower sigma (the measured
closure moves to ≈ 1.4 at K = 0.40). Nose: Cp(0.40) = +0.11 against -0.01 measured, the suction
develops later. To hit sigma exactly, shift p_v by the upstream pressure of the run:
W 61 = p_ref - sigma rho U²/2 (here -265.5 for 0.5): `python3 make_case.py 40 128 15 --pref -15.5`.

D/40, all four cases to 15 s with `--pref -15.5` (user, 2026-10-19): actual sigma 0.506, 0.392,
0.309; the upstream pressure is -12.7, -19.7, -10.9 Pa (no cavitation +3.6). Only sigma 0.4 is
steady after 10 s (vapour volume 0.366 m³ ± 0.0007); sigma 0.5 (0.16 -> 0.10 m³ from 5 to 15 s) and
sigma 0.3 (0.99 -> 0.53 m³) still shrink, and both have a mean side force (0.5: Fy -26, Fz -20 N;
0.3: Fy -50 N; no cavitation and 0.4: a few N), an asymmetric cavity. Drag (6DOF, frontal area)
Cd 0.33 without cavitation, 0.44-0.48 with. nocav and sigma 0.3 ran with 5e7ceda (before the
hans_dev 6DOF and flux fixes), sigma 0.5 and 0.4 with 0b88701. 
D/40 to 30 s, all four with 6dcd406 (user, 2026-10-20), means over 20-30 s:

| | sigma | vapour volume [m³] | Cd | mean Fy, Fz [N] |
|---|---|---|---|---|
| no cavitation | | | 0.287 | 0.1, 0.0 |
| 0.5 | 0.501 | 0.1031 (steady from 15 s) | 0.340 | 0.2, -2.2 |
| 0.4 | 0.391 | 0.364 (steady from 10 s) | 0.476 | 0.2, -2.0 |
| 0.3 | 0.307 | 0.762 (steady from 10 s) | 0.449 | -2.5, -3.0 |

All steady and axisymmetric in the mean: the shrinking cavities and the mean side forces of the
15 s runs were the start-up and the code before the hans_dev flux and 6DOF fixes. With cavitation
the forces alternate from step to step (lag-1 correlation -0.9 to -0.95, side force 40-60 N rms,
13 N or less after averaging two steps): a numerical odd-even oscillation in time of the
implicit phase-change coupling, not shedding; harmless for the means, to be removed before
pressure signals of cavitating flow are used for acoustics.
Surface pressure, 20-30 s (40 VTU frames, `headform_analysis.py --t0 20 --t1 30`), against
Rouse & McNown; closure: s/D where Cp is halfway between -sigma and the recompression peak:

| | Cp_min / plateau CFD | exp. | cavity start CFD | closure CFD | closure exp. |
|---|---|---|---|---|---|
| no cavitation | -0.61 at 0.66 | -0.74 at 0.66 | | | |
| 0.5 | -0.50 | -0.50 | 0.71 | 1.35 | 1.04 |
| 0.4 | -0.39 | -0.40 | 0.69 | 1.83 | 1.43 |
| 0.3 | -0.31 | -0.30 | 0.69 | 2.47 | 2.09 |

The cavity holds p_v (plateau at -sigma) and its length grows with falling sigma as measured, but
it is 0.3-0.4 D too long, starts downstream of the measured one (the plateau of the experiment
covers s/D 0.54) and closes sharply with a recompression peak of Cp 0.2-0.5 (measured 0.05-0.09,
the time average of an unsteady closure). Without cavitation the suction peak is 0.13 too weak
and the recovery behind the shoulder too slow (-0.47 against -0.17 at s/D 0.92): D/40 does not
resolve the flow at the shoulder, which sets the separation and the cavity detachment.
Next, the grid at the body: `make_case.py ... --amr` refines the box x in [-0.75, 2.45],
|y|, |z| < 0.7 (nose, shoulder, the cavities down to sigma 0.3) by one level (G 1 1), written to
`d<D/dx>_amr/`. The x+ side becomes the outflow `C 14 2` with `B 77 10` (patch boundaries are not
allowed with the refinement), N 40 3, G 40 1 (the patches placed for the load), G 7 0. Output in
the LAGOON store (P 18 1) with alpha_v on the patches; `headform_analysis.py` reads the finest
patches from it (with `lagoon_read.py`). Two runs: D/20 + AMR (D/40 at the body) against the
uniform D/40 above checks the refinement (`python3 make_case.py 20 64 30 --pref -15.5 --amr`), D/40 +
AMR (D/80 at the body, ~4.5 M patch cells, dt about half of D/40) is the resolution study
(`python3 make_case.py 40 256 30 --pref -15.5 --amr`). Local check (D/20 + AMR, sigma 0.4, 8
ranks): 8 patches, 0.56 M patch cells, load max/mean 1.01, ~7 s per step, dt 0.0064 s.

D/20 + AMR (D/40 at the body), sigma 0.4, with the condensation limit and W 69 (9fe9fa4, user,
2026-10-26): sigma 0.400 +- 0.001 from 5 s on, vapour volume 0.128 m³ steady from 15 s (+- 0.0002),
side force rms 0.0 N (no 2 dt oscillation), Cd 0.348. Cp plateau at -sigma from s/D 0.75 to 1.45,
closure (halfway between -sigma and the recompression peak) s/D 1.57 against 1.43 measured; peak
Cp +0.04 at 1.72 (measured +0.05 at 1.79). Before the condensation limit: closure 1.83 (uniform D/40,
sigma 0.39), the spurious pinning of traces of vapour behind the cavity to p_v made it 0.3-0.4 D too
long.

D/20 + AMR, all three sigma (9fe9fa4, user, 2026-10-27), means over 20-30 s; closure = s/D where
Cp is halfway between -sigma and the recompression peak (sigma 0.5 has no peak: halfway to the
value at s/D 1.3):

| sigma | closure CFD | measured | uniform D/40 before the condensation limit | vapour [m³] | Cd |
|---|---|---|---|---|---|
| 0.500 | 1.10 | 1.04 | 1.35 | 0.008 (alpha_v at the wall <= 0.42) | 0.314 |
| 0.400 | 1.56 | 1.43 | 1.83 | 0.128 | 0.348 |
| 0.301 | 2.43 | 2.09 | 2.47 | 0.63, still growing 1 % per s | 0.345 |

The trend with sigma is right and sigma is held to 0.001; the cavities are 0.06, 0.13 and 0.34 D too
long, the longest the most (sigma 0.5 corrected from 1.14: that value took the measured
threshold). The CFD cavities are steady (side force rms 0.1 N) with a sharp closure,
the measured long cavity closes over s/D 1.8-2.4 (a time average of an unsteady closure: the
re-entrant jet sheds the rear of the cavity). k-epsilon without a correction for the mixture
(Reboud et al. 1998: mu_t f(rho), f = rho_v + (1 - alpha)^n (rho_l - rho_v), n = 10) is known to damp
the re-entrant jet and give too long, steady cavities.

With the Reboud correction (W 70 10, D/20 + AMR, 66867d6, user, 2026-10-28), means over 15-30 s:

| sigma | closure CFD | without W 70 | measured | vapour [m³] (without) | alpha_v at the wall | Cd |
|---|---|---|---|---|---|---|
| 0.501 | 1.10 | 1.10 | 1.04 | 0.004 (0.008) | <= 0.24 | 0.313 |
| 0.401 | 1.56 | 1.56 | 1.43 | 0.06 +- 0.007 (0.128) | ~0.5 | 0.349 |
| 0.302 | 2.44 | 2.43 | 2.09 | 0.30 +- 0.03 (0.63) | ~0.48, falling from s/D 1.6 | 0.345 |

(closure here for all runs: halfway between -sigma and the recompression peak, for sigma 0.5 the
value at s/D 1.3.) The correction halves the vapour volume and makes the cavities unsteady (smooth
low-frequency changes of the volume and side forces of 2-13 N rms at sigma 0.4 and 0.3, no 2 dt
oscillation; the time-averaged alpha_v at the wall ~0.5) and the rise of the pressure at the closure
more gradual, as measured, but the mean closure point does not move. The remaining 0.06-0.35 D are
not the eddy viscosity of k-epsilon at this resolution; next: D/40 + AMR (D/80 at the body).

A restart from the state files (P 40/I 40) does not help, they do not hold the vapour fraction;
the checkpoints (P 48/I 45) do, a run with `P 48` can be continued.

## 03 Vapour layer with mesh refinement

The layer of 01 in 2D (y one cell, dx = 5 mm) with the outflow `B 77 10` instead of the patch,
fixed dt 2.5e-4 s, 0.25 s. Collapse time against the piston equation (0.1825 s):

| grid | cells | collapse | wall time (1 rank) |
|---|---|---|---|
| uniform 5 mm | 1600 | 0.1825 s (+0.0 %) | 7 s |
| uniform 2.5 mm | 6400 | 0.1827 s (+0.2 %) | 27 s |
| 5 mm + `G 1 1`, `G 15 0 0.3 0 0.005 0 0.04` (2.5 mm patch over the layer) | 1600 + 1920 | 0.1817 s (−0.4 %) | 22 s |

The same on 2 ranks (the patch on one of them): 0.1817 s. The layer thickness follows the piston
curve within 0.6 mm in all three. Regression: `Regression/cases/cfd_2d_cavitation_slab_amr`
(2 ranks, 40 steps).
