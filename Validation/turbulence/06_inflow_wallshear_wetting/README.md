# 06 Inflow turbulence, wall shear per RK stage, NHFLOW wetting fronts (patches 0020–0025)

## What it checks

| patch | change | case |
|---|---|---|
| 0020 | CFD discharge inflow (B 60): k, ε/ω profile per inflow column with the depth to the interpolated interface (was the top face of the highest wet cell over the whole inflow), bed ks as `roughness::ks_val` (B 55, S 28) | `cfd_2d_channel_ke` |
| 0021 | CFD explicit wall shear with the velocity of the RK stage (was uⁿ in every stage; PLIC: a partly updated u in stage 3) | `cfd_2d_channel_ke` (steady: no change expected) |
| 0022 | IO flags at patch BC faces (the cell was taken from the previous loop) | regression `cfd_2d_patch_inflow_komega` |
| 0023 | Input checks T 10, T 12, T 21, A 560 | – |
| 0024 | NHFLOW turbulence output (vertex averaging was shifted by half a cell in x), print_2D, CDS2 cut-offs, unused functions | – |
| 0025 | NHFLOW k-ε/k-ω: newly wetted columns start from the wet neighbours; length-scale limit l ≤ κh (A 564 ≥ 1) | `nhflow_3d_bank_ke`, `nhflow_2d_bore_beach_ke` |

- **before:** hans_dev ed3c2c66b
- **after:** + 0020–0025
- Built in a Linux sandbox, g++ 13 `-O2`, single rank.

## Cases

- `cfd_2d_channel_ke`: CFD 2D open channel, 12 m × 0.24 m water, U = 0.5 m/s, ks = 1 mm, k-ε, RK3, 40 s.
- `nhflow_3d_bank_ke`: NHFLOW 3D, 200 m × 40 m, flat bed for y < 20 m, bank (T 62) rising to 2 m at y = 40 m,
  still water 1.5 m, 41 m³/s, k-ε, A 519 1, 600 s. The water level rises during the spin-up, so the columns
  between y ≈ 33 and 36 m are wetted during the run.
- `nhflow_2d_bore_beach_ke`: NHFLOW 2D, 40 m, 1:20 beach from x = 20 m, still water 0.5 m, a 1.0 m box at
  x < 5 m (F 72) released as a bore that runs up and over the beach, k-ε, 40 s. Analysed from the regression-dump
  states (cell values) written every 100 steps (`REEF3D_REGRESSION_EVERY=100`).

## Results

**NHFLOW bank, columns wetted during the run** (`../tools/nhf_bank.py`, depth-averaged ν_t over x = 100–180 m
divided by the local equilibrium κu*h/6 with u* = u κ/ln(11h/ks)):

| y (m) | h (m) | before | after | upstream x = 10–40 m (after) |
|---|---|---|---|---|
| 20–26 | 1.45–0.92 | 1.09–1.25 | 1.09–1.25 | 0.88–0.96 |
| 28 | 0.72 | 1.38 | 1.38 | 1.00 |
| 30 | 0.52 | 1.20 | **1.40** | 1.14 |
| 32 | 0.32 | 0.42 | **1.01** | 1.17 |
| 34 | 0.12 | 0.09 | 0.12 | 0.66 |

Before, a column wetted during the run started with k = ε = 0. That is a fixed point of k-ε: ν_t = 0 gives no
production. The column only got turbulence by diffusion from its neighbours. After 200 steps (regression case),
columns near the inflow with 0.1–0.5 m of water still had k = ν_t = 0 in the top layer. Further downstream, after
600 s, ν_t at y = 32 m was 0.42 of the local equilibrium; it is now 1.01, close to the upstream value. The flow
(u, h) changes by less than 0.4 %. In the 200-step regression run the wetting start alone gives this change; the
length-scale limit alone changes ν_t by up to 8 % of its maximum.

**NHFLOW bore on a dry beach** (`../tools/nhf_front.py`, columns with 1 cm < h < 10 cm):

| | before | after |
|---|---|---|
| largest l/h = ν_t/(cμ^¼ √k h) | 0.71 (t = 17.8 s) | 0.40 (the limit κ) |
| median over time of the largest l/h | 0.19 | 0.19 |
| largest ν_t, wet / shallow | 6.5e-3 / 2.3e-3 m²/s | 6.5e-3 / 2.3e-3 m²/s |
| largest run-up x, shoreline at t = 25 / 30 / 40 s | 39.95 m, 32.35 / 28.85 / 32.05 m | the same |

The bore carries its turbulence with it, so here the limit clips one peak in the thin swash layer; run-up and ν_t
are unchanged.

**CFD channel** (`../tools/cfd_channel_x.py`, depth-averaged ν_t/ν_eq and k/k_eq, U = 0.5 m/s, h = 0.24 m):

| x (m) | 0.5 | 2 | 6 | 11.8 |
|---|---|---|---|---|
| before | 0.999 / 1.003 | 0.951 / 0.964 | 0.918 / 0.962 | 0.924 / 1.044 |
| after | 1.016 / 1.011 | 0.963 / 0.971 | 0.922 / 0.964 | 0.924 / 1.043 |

The water depth at the inflow is 0.2431 m. Before, H was the top face of the highest wet cell (0.24 m), which
gave a profile slightly too thin in the upper layers (k at z = 0.22 m 4.7 % lower). The change fades out
downstream (< 0.3 % from x = 6 m). The wall shear with the stage velocity does not change the steady state
(u near the bed +0.04 %).

## How to run

From this folder:

```
REEF3D_REGRESSION_EVERY=100 ../tools/run_all.sh <REEF3D> <DiveMESH> cases run
python3 ../tools/cfd_channel_x.py run/cfd_2d_channel_ke
python3 ../tools/nhf_bank.py run/nhflow_3d_bank_ke
python3 ../tools/nhf_front.py run/nhflow_2d_bore_beach_ke
```

About 15 min for the CFD case and 4 min for each NHFLOW case on one core.
