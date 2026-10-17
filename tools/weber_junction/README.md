# The 90 deg open-channel junction of Weber et al. (2001) on three blocks (P6a)

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
| `make_regression_cases.py` | the regression cases of P6a (below) |

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

## Notes

- the channel sloshes after the start (still water, the inflows switched on): a period of about 70 s, little
  damping; steady after 250 - 300 s with the eddy viscosity
- without eddy viscosity the depth-averaged flow separates at the downstream corner of the branch and stays separated
  for about 12 W (measured: 2 - 3 W); the flow squeezed past it, standing waves downstream, the rise across the
  junction twice the measured one.  Where the separation reaches the outlet, water flows back in through it and the
  run fails, with the water level on the face of a block (layout a, RUN5 without eddy viscosity) as with the native
  outflow (layout b, RUN2 with the parabolic eddy viscosity): neither outflow limits the backflow
- the native k-epsilon of SFLOW (`A 260 1`) on the Cartesian reference went unstable (6.8 m/s after 50 s)
