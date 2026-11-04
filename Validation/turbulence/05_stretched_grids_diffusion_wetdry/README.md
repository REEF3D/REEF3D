# 05 Stretched grids, NHFLOW momentum diffusion, wet/dry turbulence (patches 0015–0019)

## What it checks

| patch | change | case |
|---|---|---|
| 0015 | Implicit first-order upwind (CFD `ifou`, NHFLOW conservative `ifou`, SFLOW `ifou`) per face and conservative with the cell's own spacing (was one upwind direction per cell and DX[i−1] for positive flow) | `cfd_2d_channel_ke_uniform` / `_stretched` (T 12 1) |
| 0016 | NHFLOW momentum diffusion: σ metrics at the faces, DZP/DYP in the central differences, factor 2 in w only on σz², `nhflow_idiff_v_2D` filled, `nhflow_ediff` rewritten (A 512 1) with a viscous time-step limit | `nhflow_2d_bump_visc_a512_1` / `_2` |
| 0017 | CFD wall details: VRANS wall-law cells, no wall cells at open boundaries, T 44 sampling height, ε*/ω* clamp, T 31, D 20 0 warning | (regression A/B) |
| 0018 | NHFLOW turbulence: zero gradient towards dry columns (was k = ε = 0), implicit boundary elimination, relaxation of ε/ω and EV0, inflow ghosts, σ update | `nhflow_3d_bank_ke` |
| 0019 | `nhflow_ediff`: the wall rule A 513 of the implicit form; dry cells skipped | `nhflow_2d_bump_visc_a512_1` / `_2` |

- **before:** hans_dev 4c91acff2
- **after:** + 0015–0019
- Built in a Linux sandbox, g++ 13 `-O2`, single rank.

## Cases

- `cfd_2d_channel_ke_uniform` / `_stretched`: CFD 2D open channel, 12 m × 0.24 m water, U = 0.5 m/s, ks = 1 mm,
  k-ε with `T 12 1` (implicit first-order upwind for k and ε), 40 s. The stretched case uses `B 101 1`, `B 111 2.0`.
- `nhflow_2d_bump_visc_a512_1` / `_2`: NHFLOW 2D, 300 m × 2 m, bump 0.5 m high at x = 100–200 m, laminar
  ν = 0.01 m²/s (`W 2`), bed shear A 519 1, 600 s; explicit (A 512 1) against implicit (A 512 2) diffusion. Laminar
  because A 512 is forced to 2 with any turbulence model (`control_parse.cpp`).
- `nhflow_3d_bank_ke`: NHFLOW 3D, 200 m × 40 m, flat bed for y < 20 m, bank (T 62) rising to 2 m at y = 40 m,
  h = 1.5 m (shoreline at y = 35 m), 41 m³/s, k-ε, A 519 1, 600 s.

## Results

**NHFLOW explicit vs implicit diffusion** (`../tools/nhf_compare.py`, x = 100–250 m, relative to the maximum of the
implicit run):

| | u max / mean | w max / mean | free surface (of the depth) |
|---|---|---|---|
| before | 50 % / 19 % | 71 % / 14 % | 5.9 % |
| + 0015–0018 | 50 % / 19 % | 71 % / 14 % | 5.9 % |
| + 0019 | **1.6 % / 0.5 %** | 2.5 % / 0.5 % | 0.03 % |

At x = 60 m: explicit before u_bed = 0.065 m/s, u_surface = 1.30 m/s, h = 2.195 m; after 0.366 / 1.066 / 2.061 m;
implicit 0.360 / 1.072 / 2.061 m. The explicit form used the bed ghost value of A 518 2 (zero), so the bed was a
no-slip wall on top of the bed shear. The stencil corrections of 0016 alone change the explicit run by < 0.01 % in u and
0.4 % in w, and the implicit run by < 0.002 % here (2D, nearly uniform σ). The remaining 1.6 % is the time discretisation
(explicit RK stages against implicit Euler per stage).

**CFD ifou on a stretched grid** (`../tools/cfd_channel_x.py`, depth-averaged ν_t/ν_eq and k/k_eq):

| x (m) | 2 | 6 | 10 | 11.8 |
|---|---|---|---|---|
| uniform, before = after | 0.950 / 0.964 | 0.916 / 0.962 | 0.927 / 0.982 | 0.923 / 1.044 |
| stretched, before = after | 0.950 / 0.963 | 0.920 / 0.964 | 0.940 / 0.989 | 0.938 / 1.053 |

Largest field change before → after: k 4.2e-4, ν_t 6.8e-4 of the maximum on the stretched grid (8e-5 / 1.6e-4 on
the uniform grid). In a channel near equilibrium the streamwise convection of k and ε is small, so the
non-conservative form had little effect here. The fix matters where k or ε change strongly along a stretched
direction (wakes, inflow transitions). The 1.5 % difference between the grids at x = 10 m is the coarser
resolution, not `ifou`.

**NHFLOW wet/dry front, 3D bank** (`../tools/nhf_bank.py`, depth-averaged over x = 100–180 m):

| y (m) | h (m) | u (m/s) | ν_t before | ν_t after | k before | k after |
|---|---|---|---|---|---|---|
| 20 | 1.45 | 0.998 | 5.70e-3 | 5.70e-3 | 5.55e-3 | 5.55e-3 |
| 28 | 0.72 | 0.682 | 2.70e-3 | 2.70e-3 | 3.73e-3 | 3.74e-3 |
| 30 | 0.52 | 0.449 | 1.23e-3 | 1.17e-3 | 1.76e-3 | 1.71e-3 |
| 32 | 0.32 | 0.237 | 2.07e-4 | 1.44e-4 | 3.99e-4 | 3.40e-4 |
| 34 | 0.12 | 0.126 | 9.0e-6 | 7.6e-6 | 5.05e-4 | 5.14e-4 |

The change is confined to the last two to three wet rows. Before, the dry neighbour columns acted as k = ε = 0
Dirichlet values. ε = 0 at the shoreline lowered ε in the shallow cells and raised ν_t = Cμ k²/ε. Zero gradient
removes that, and ν_t drops by 30 % at y = 32 m. The flow is unchanged within 0.2 %. Neither run is close to the
local equilibrium ν_t,eq = κu*h/6 in the shallowest rows (0.6 before, 0.4 after at y = 32 m). The shoreline is
only two cells wide, and the vertex output averages over the dry cells. This case documents the effect; it is not
an accuracy test.

## How to run

From this folder:

```
../tools/run_all.sh <REEF3D> <DiveMESH> cases run
python3 ../tools/cfd_channel_x.py run/cfd_2d_channel_ke_uniform run/cfd_2d_channel_ke_stretched
python3 ../tools/nhf_compare.py run/nhflow_2d_bump_visc_a512_2 run/nhflow_2d_bump_visc_a512_1 100 250
python3 ../tools/nhf_bank.py run/nhflow_3d_bank_ke
```

About 15–20 min per CFD case, 1.5 min per NHFLOW 2D case and 4 min for the 3D bank on one core.
