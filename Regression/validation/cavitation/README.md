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

Setup notes
- The pressure level must be fixed by a Dirichlet boundary: the outlet side a wall `C 14 21`
  with an outlet patch over it, `B 441 id 4 xs xe ys ye zs ze` and `B 412 id 0.0` (pressure
  correction = patch pressure − cell pressure), with `B 60 1`. p_v is relative to it:
  W 61 = -sigma rho U²/2 for the cavitation number sigma = (p_inf - p_v)/(rho U²/2), p_inf = 0.
  Not suitable: the plain outflow with `B 77 10` (since d9ac54c40 the correction there is
  outflow pressure − ghost pressure = 0, the level stays free: in case 01 the column then sits at
  p_v and the layer shrinks at a constant 0.33 m/s), `B 60 0` (no outlet pressure at all) and
  the gauge `B 30 1` (shifting the level after each projection against the implicit phase-change
  term diverges at the first step).
- Gravity off (`W 22 0`) unless the hydrostatic pressure is meant to be part of p - p_v.
- Single phase (`F 30 0`, `F 80 0`); a fluid at rest needs a fixed time step (`N 48 0`, `N 49 dt`).
- Output: `alpha_v` in the VTU files (`P 24 1` adds the mixture density), and
  `REEF3D_CFD_Cavitation/REEF3D-CFD-Cavitation.dat`: time, vapour volume, equivalent radius
  (3V/4π)^(1/3), vapour production rate, minimum pressure.

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
W 61 = p_ref - sigma rho U²/2 (here -265.5 for 0.5).
