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
