# The 90 deg open-channel junction of Weber et al. (2001) on three blocks (P6a, P6b, P6c)

The validation case of the multi-block grids (P5 gate): the IIHR flume of Shumate (1998) and Weber, Schumate & Mawer
(2001, J. Hydraul. Eng. 127(5)) as published by Huang, Weber & Lai (2002, J. Hydraul. Eng. 128(3), 268-280).

- main channel 21.946 m (72 ft), branch 3.658 m (12 ft) joining at 90 deg 5.486 m (18 ft) downstream of the main
  inlet, both W = 0.914 m (3 ft) wide, horizontal bed
- total discharge Qt = 0.170 m3/s; RUN2 q* = Qm/Qt = 0.25 (Qm 0.043, Qb 0.127), RUN5 q* = 0.75 (Qm 0.127, Qb 0.043)
- depth at the downstream end of the main channel (7 W below the junction, Huang's H0): 0.3054 m (RUN2), 0.3100 m
  (RUN5); Fr 0.35
- coordinates of the measurements: x/W positive upstream, the branch at -1 < x/W < 0; y/W from the bank of the
  branch (0) to the opposite bank (1)

## Files

| file | |
|---|---|
| `weber_case.py` | the case on three blocks: `blocks.json` (LAGOON's block topology, checked by `lagoonviz.core.multiblock.check`), `blocks.dat`, `control.txt`, `ctrl.txt`; `--layout a`: block 0 the main channel upstream (native inflow), block 1 the junction and the main channel downstream (water level `G 54 1 2`); `--layout b`: block 0 the main channel downstream (native outflow), block 1 the main channel upstream and the junction (inflow `G 54 1 1`); the branch block 2 (inflow `G 54 2 1`) on a part of block 1's +j face, a vertex of three blocks at a corner of the branch; `--h0` the water level at the outlet, `--xctrl` more lines (`"A 260 3"`) |
| `weber_cart.py` | the same flume on ONE Cartesian grid with the same grid lines (native SFLOW): solids (`S 10`) beside the branch, the branch inflow a patch (`B 440` face 2, `B 411 -Qb`: the discharge of a patch counts along the axis) |
| `weber_compare.py` | runs against the measurements: the water surface along the main channel (Fig. 5 of Huang et al., RUN5), the rise across the junction, the minimum in the separation zone, the reversed depth-averaged flow; `--png` the profiles |
| `data/run5_surface_fig5.json` | the measured water surface of RUN5 digitized from Fig. 5 of Huang et al. (2002): 7 lines y/W 0.020 ... 0.979, 252 points, about +-0.003; x/W from the 3 in spacing of the points (-4 + (px - 211)/234.4), checked against the contours of Fig. 3 |
| `weber_lib.py` | the loaders of both kinds of run (the AMR dump of the blocks, the VTP of the Cartesian grid) in Huang's coordinates |
| `weber_fig8.py` | NHFLOW runs against the streamwise velocity in the cross section x/W = -2 (Fig. 8 of Huang et al., RUN5 and RUN2): rms (the model clipped to the colour bar), the width of the slow zone (u/U0 = -0.5) near the bed and the surface; `--png` the measured and the computed section and three profiles |
| `data/fig8_u_xm2.json` | u/U0 at x/W = -2 digitized from the colours of Fig. 8 (Experiment), RUN5 495 points, RUN2 416 points, 16 bins of 0.1 (about +-0.1), saturated at -1.4; U0 = 0.628 m/s |
| `weber_fig910.py` | NHFLOW runs against the vertical velocity profiles of Figs. 9 (streamwise) and 10 (transverse) of Huang et al. (RUN5, stations x/W 1, 0, -1, -2, -6 at y/W 0.25, 0.5, 0.75 and three in the branch): rms and mean of model - measured per station and figure; `--png` all profiles |
| `digitize_fig9_10.py` | digitizes the dots of Figs. 9 and 10 from the 300 dpi images of the PDF (`pdfimages -f 8 -l 9 -png`): the curve of the model removed row by row, isolated dots and dots on the curve; writes `data/fig9_10_profiles.json` |
| `data/fig9_10_profiles.json` | the measured profiles of Figs. 9 and 10 (Shumate 1998, RUN5): 33 stations, 9-13 points each, about +-0.05 U0 (Fig. 9) and +-0.015 U0 (Fig. 10); Fig. 10 read with the scale 0.6 (its bar reads 6.0, which gives 3.4 U0 at the junction, ten times the branch inflow), x/W = 1 of Fig. 10 not digitized (the dots lie on the zero line) |
| `make_regression_cases.py` | the regression cases of P6a, P6b and P6c (below); keeps the references of the cases it rewrites |

## Results (SFLOW, W/dx = 16, 7168 cells, 400 s)

Water level at the outlet set so that the depth 7 W below the junction is Huang's H0 (`--h0 0.313` RUN5,
`0.3084` RUN2); z/H0 over the mean depth across x/W = -4; rise: mean z/H0 across x/W = +1; the surface minimum below
the corner cell of the branch.

| run | eddy viscosity | rise | minimum | reversed flow (depth-averaged) | Fig. 5: mean / rms / max |
|---|---|---|---|---|---|
| RUN5 measured | | +5.0 % (text), 1.057 (Fig. 5) | 0.95 | | |
| RUN5 blocks | none (100 s, tailwater 0.296) | +9.5 % | 0.91 (corner) | 13 W long, 0.25 W wide | +0.017 / 0.026 / 0.071 |
| RUN5 blocks | parabolic `A 260 4` | +7.4 % | 0.986 | 2.8 W, 0.06 W | +0.012 / 0.016 / 0.044 |
| RUN5 blocks | mixing length `A 260 3` | +6.2 % | 0.974 | none | +0.006 / 0.012 / 0.062 |
| RUN5 blocks, layout a | mixing length | +6.2 % | 0.974 | none | +0.006 / 0.012 / 0.062 |
| RUN5 one Cartesian grid | mixing length (native) | +6.7 % | 0.984 | none | +0.009 / 0.013 / 0.038 |
| RUN2 measured | | +8.8 % | 0.92 | | |
| RUN2 blocks | mixing length | +14.5 % (+11.8 % over the depth at x/W = -7) | 0.955 | 4.0 W, 0.19 W | |
| RUN2 blocks | parabolic | unstable after 110 s: the separation reached the outlet (17 W), water flowed back in | | | |

- the blocks and the single Cartesian grid with the same grid lines agree (rise 6.2 / 6.7 %; the Cartesian
  mixing length takes its gradients from the face velocities, the profiles of the Cartesian grid are the node
  values of the VTP output); the water level on the face of a block (layout a) and the native outflow (layout b)
  give the same solution
- RUN5 (main channel dominant): the shape and the level of the surface across the junction as measured with the
  mixing length; the depression in the separation zone (0.95 - 0.97 at x/W -2 ... -1) is too shallow (0.98 - 0.99)
- RUN2 (branch dominant): the rise across the junction 1.3 - 1.6 times the measured one, the depression too
  shallow; the depth-averaged model has none of the helical flow below the junction (Huang et al., Fig. 11), the
  3D NHFLOW run with turbulence is the next comparison

## Results P6b

SFLOW k-epsilon (`A 260 1`) on the Cartesian grid (`weber_cart.py --turb 1`, W/dx 16, 400 s), NHFLOW k-omega on the
three blocks (`weber_case.py --solver nhflow --layout b --n 12 --nz 8 --xctrl "A 560 2"`, hydrostatic, 300 s):

| run | model | rise | minimum | reversed flow | Fig. 5: mean / rms / max | Fig. 8 (x/W = -2): rms |
|---|---|---|---|---|---|---|
| RUN5 measured | | +5.0 % | 0.95 | none at x/W = -2 | | |
| RUN5 one Cartesian grid | SFLOW k-epsilon | +5.4 % | 0.967 | none | -0.001 / 0.008 / 0.025 | |
| RUN5 blocks | NHFLOW k-omega, 8 layers | +5.7 % | 0.967 | none | +0.005 / 0.013 / 0.073 | 0.25 |
| RUN2 measured | | +8.8 % | 0.92 | at the surface, y/W < 0.2 | | |
| RUN2 one Cartesian grid | SFLOW k-epsilon | +7.5 % | 0.876 (corner cell) | 0.9 W, 0.09 W | | |
| RUN2 blocks | NHFLOW k-omega, 8 layers | +11.4 % | 0.933 | top layer 1.0 W x 0.08 W, bottom layer 1.5 W x 0.17 W | | 0.18 |

- the depth-averaged k-epsilon is the closest of the SFLOW models (RUN5 Fig. 5 rms 0.008 against 0.012 with the
  mixing length; RUN2 rise +7.5 % against +14.5 %); before P6b it went unstable (below)
- NHFLOW RUN5: the surface as SFLOW k-epsilon; in the section x/W = -2 the slow zone next to the branch-side bank
  has the measured width at the bed (u/U0 = -0.5 at y/W 0.13) but not at the surface (0.04 against 0.18), the core
  is slower than measured (-1.18 against -1.3 ... -1.4): the shear layer diffuses over 2 - 3 cells (W/dx 12)
- NHFLOW RUN2: the reversed flow is larger at the bed than at the surface, measured is the opposite (Fig. 8: reversed
  at the surface for y/W < 0.2, none at the bed); the hydrostatic model with 8 layers misses the vertical structure
  of the separation (the helical flow below the junction); next: non-hydrostatic (`A 520 2`) and W/dx 16

## Results P6c

The runs of P6a and P6b had no bed friction: REEF3D's defaults `A 218 0` (SFLOW) and `A 519 0` (NHFLOW) leave the
momentum without it, `B 50` acts only on the turbulence models (REEF3D now warns, `!!!`).  P6c: `weber_case.py` and
`weber_cart.py` set `A 218 1` (SFLOW, Manning with the `ks` of the same Chezy at h0 as the log law of the Nikuradse
`--ks` 1 mm, C 64) and for NHFLOW `A 519 1` (the log law of the bottom layer), `B 61 2` (log-law profiles at both
inflows) and `A 567 3` (the free-surface value of omega, Celik & Rodi y' 0.07 h: without it k-omega kept its eddy
viscosity up to the surface, 0.004 m2/s behind the junction, and mixed the shear layer over the depth).  NHFLOW
k-omega, W/dx 12, 8 layers, 300 s; SFLOW W/dx 16, 400 s:

| run | model | rise | minimum | reversed flow (top / bottom layer) | Fig. 5 rms | Fig. 8 rms | slow zone u/U0 = -0.5 at y/W, bed / top (measured 0.13 / 0.18, RUN2 0.16 / 0.31) | Fig. 9 / 10 rms (U0) |
|---|---|---|---|---|---|---|---|---|
| RUN5 measured | | +5.0 % | 0.95 | none at x/W = -2 | | | | |
| RUN5 | SFLOW k-epsilon, Cartesian | +5.5 % | 0.970 | none | 0.008 | | | |
| RUN5 | SFLOW mixing length, blocks | +6.2 % | 0.977 | none | 0.012 | | | |
| RUN5 | NHFLOW non-hydrostatic | +6.1 % | 0.975 | 0.46 W x 0.08 W / none | 0.013 | 0.23 | 0.11 / 0.19 | 0.094 / 0.051 |
| RUN5 | NHFLOW hydrostatic | +5.8 % | 0.974 | none | 0.013 | 0.22 | 0.13 / 0.17 | 0.090 / 0.046 |
| RUN5 | P6b: NHFLOW hydrostatic, no friction, no surface omega | +5.7 % | 0.967 | none | 0.013 | 0.25 | 0.13 / 0.04 | 0.095 / 0.049 |
| RUN2 measured | | +8.8 % | 0.92 | larger near the surface (Figs. 7, 8) | | | | |
| RUN2 | SFLOW k-epsilon, Cartesian | +7.5 % | 0.885 | 0.9 W x 0.09 W | | | | |
| RUN2 | SFLOW mixing length, blocks | +13.6 % | 0.956 | 3.4 W x 0.19 W | | | | |
| RUN2 | NHFLOW non-hydrostatic | +11.3 % | 0.942 | 2.6 W x 0.17 W / 1.5 W x 0.08 W | | 0.145 | 0.20 / 0.29 | |
| RUN2 | NHFLOW hydrostatic | +11.0 % | 0.936 | 2.3 W x 0.17 W / 1.3 W x 0.17 W | | 0.146 | 0.20 / 0.29 | |
| RUN2 | P6b: NHFLOW hydrostatic, no friction, no surface omega | +11.4 % | 0.933 | 1.0 W x 0.08 W / 1.5 W x 0.17 W | | 0.18 | 0.23 / 0.27 | |

- with the bed friction, the log inflow and the surface value of omega the separation has the measured vertical
  structure: larger near the surface than near the bed (RUN2: 2.6 W against 1.5 W; P6b had it the other way round),
  and the slow zone at x/W = -2 has the measured width at the bed and at the surface for both runs
- non-hydrostatic and hydrostatic agree within the scatter of the comparison (the non-hydrostatic RUN5 has a short
  reversed zone at the surface); the profiles upstream and in the branch (Fig. 9 a, b, f) within 0.02 - 0.05 U0
- left: the core of the main channel at x/W = -2 is 0.15 - 0.2 U0 slower than measured (Fig. 9 d, y/W 0.25 and 0.5),
  the transverse velocity at x/W = -1 near the bed 0.1 U0 too small (Fig. 10 c), the rise of RUN2 +11 % against
  +8.8 %; the next steps are W/dx 16 - 24 and more layers (the shear layer over 2 - 3 cells), a check of the depth
  H0 (all runs are about 1.5 % deeper than measured at x/W = -7) and the side-wall friction (`A 519 2`, refused on
  blocks)
- SFLOW: the bed friction changes the results little (RUN5 k-epsilon +5.4 -> +5.5 %, mixing length +6.2 % both)
- the P6b wall lists made the NHFLOW surface dip by 6 mm in the first column of block 1 behind the interface across
  the channel (x/W = -1); all NHFLOW runs above are with the P6c wall lists (`../nhflow_multiblock/README.md`)

## Notes

- the channel sloshes after the start (still water, the inflows switched on): a period of about 70 s, little
  damping; steady after 250 - 300 s with the eddy viscosity
- without eddy viscosity the depth-averaged flow separates at the downstream corner of the branch and stays separated
  for about 12 W (measured: 2 - 3 W); the flow squeezed past it, standing waves downstream, the rise across the
  junction twice the measured one.  Where the separation reaches the outlet, water flows back in through it and the
  run fails, with the water level on the face of a block (layout a, RUN5 without eddy viscosity) as with the native
  outflow (layout b, RUN2 with the parabolic eddy viscosity): neither outflow limits the backflow
- the native k-epsilon of SFLOW (`A 260 1`) on the Cartesian reference went unstable (6.8 m/s after 50 s): the
  face velocity of the branch inlet (a patch, `B 440`) stayed 0 (the wall treatment of the faces), the shear
  across the inlet row produced k (0.4 m2/s2 and nu_t 0.14 m2/s in the branch at W/dx 8, where the bed friction
  gives 0.003); P6b sets the face velocities of the patch in- and outlets from the cell and its ghost cell
  (`sflow_momentum_func::face_velocities`): k 0.003 in the branch, stable
- backflow through the outflow (RUN2 parabolic, layout b; RUN5 without eddy viscosity, layout a): P6b clips the
  inflowing normal velocity in the ghost cells of the SFLOW outflow (`outflow_clip`, with `B 60` >= 1, the native
  outflow and the water level on the face of a block `G 54 b 2`); RUN2 parabolic now runs through 400 s (umax
  1.6 m/s), but its separation zone still reaches the outlet (the parabolic model is not suited to RUN2); RUN5
  without eddy viscosity on layout a through 300 s
- the wall lists of the blocks (P6b, `../nhflow_multiblock/README.md`): the corner of block 1's +j face where the
  wall ends and the branch begins; NHFLOW drove a jet of 4 m/s along it before P6b, SFLOW changes there slightly
  (RUN5 mixing length after 400 s: the minimum 0.9741 -> 0.9743, the rise and Fig. 5 unchanged)

## Regression cases (`make_regression_cases.py`)

| case | what |
|---|---|
| `sflow_2d_curv_channel_rot30_mixlen` | P6a: the curvilinear stress and the mixing length on a rotated river reach |
| `sflow_2d_mb_amr_parab`, `sflow_2d_mb_amr_mixlen_mpi2` | P6a: the eddy viscosity on four turned blocks with two refinement levels |
| `sflow_2d_mb_weber_junction`, `_mpi2` | P6a: the junction on three blocks (W/dx 8, parabolic), 200 steps |
| `nhflow_3d_mb_quadrants_kw` | P6b: k-omega on four turned blocks (the fills of the blocks, the composite implicit diffusion) |
| `nhflow_3d_mb_amr_ke_mpi2` | P6b: k-epsilon on the blocks and a refinement box over their vertex, 2 ranks |
| `nhflow_3d_amr_adaptive_kw` | P6b: k-omega on regridded patches (fresh patches from old ones or from the parent) |
| `nhflow_3d_mb_weber_junction_kw` | P6b: the junction with NHFLOW and k-omega (W/dx 6, 4 layers), the wall lists at the corner of the branch |
