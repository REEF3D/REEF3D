# Open faces of the blocks and the Nidelva - harbour case (P5d)

Phase P5d of the multi-block curvilinear concept: in- and outflow on the blocks 1.. of a multi-block grid
(`G 50 1`, DIVEMesh `R 20 1`) for REEF3D::SFLOW and REEF3D::NHFLOW, and the application to the lower
Nidelva and the harbour of Trondheim.

## Open faces: `G 54 b t v`

| key | meaning |
|---|---|
| `G 54 b 1 Q` | block `b` (1..): discharge `Q` (m3/s) through its -i face |
| `G 54 b 2 WL` | block `b` (1..): water level `WL` (m) at its +i face (start ramp `F 64` from `F 60` as the native outflow) |

One line per face; the face must have no interface (a wall otherwise).  The in- and outflow of the native grid
(block 0) stay as before (`B 60`, `W 10`, `F 62`, DIVEMesh `C 11`, `C 14`).  An open face is treated as the
native in- and outflow with `B 60 1`:

- 2D lists: the boundary cells of the face get bc 1 / 2 (`reefamr::build_bc2D`, as `gcsl_setbc` on a native
  grid whose first cell is the first cell of the block), so `gcslin`, `gcslout` and `IOSL` hold them.
- `ioflow_amr_block` (one per block with an open face) does the in- and outflow of the block as `ioflow_f`:
  SFLOW the velocity `Q / A` along the face normal in the ghost cells and the water level in the last cells and
  their ghost cells; NHFLOW the plain (`B 61 1`) or log profile (`B 61 2, 4, 5`, scaled to `Q`), the inflow
  level from the mean surface of the wet inflow cells (`F 50 2, 4`), the outflow level, and `Uo` of the river
  outflow of `start4V`.  A block lies on one rank: the sums over the face are local.
- SFLOW: the HLL fluxes of the in- and outflow faces from the ghost state (`sflow_HLL::open_io`).
- NHFLOW: the patch lexer of the block has `B 60 1` with `river_xm`, `river_xp` (new lexer flags, 1 by default)
  telling which of its i-faces are in- and outflow; `start*V` and `nhflow_curv` (pair_bc, iwall, bc_fluxes) read
  them; the 3D lists `gcin`, `gcout` and the `IO` marks of the ghost cells (pressure rows) of the open faces
  (`nhflow_amr::open_lists`).
- No refinement within `ioband` (4) cells of an open face (`reefamr::blocks_refine_setup`).
- Sediment (P8a): the cells beyond an open face carry the sediment flag, as the ghost cells of the native in- and
  outflow (`sediment_f::dfbed_open`, kept by `active_sflow`): the bed load leaves through the outflow face with the
  upwind flux of the last cell (before P8a the face was closed for the sediment and the end of an outflow block
  filled with sand).  As on the native grid, no sediment enters through an inflow face (its face velocity of the
  sediment is 0); hold the bed of the inflow reach with `S 71`, or (P8b) with `S 76 xs xe ys ye`: the bed of the
  cells in the box is held (their Exner change is 0, `sediment_f::ini_guard`), their bed load is computed as
  everywhere and carried into the cells downstream: the supply of the reach at its capacity.  Any number of boxes,
  physical coordinates of the cell centres (curvilinear grids and the patches of the blocks).  As for the inflow
  guard `S 78`, the slope flux of `S 86` is closed at the faces of the held cells; otherwise a box at the inflow
  is the x window `S 71` (regression `sflow_2d_mb_sed_supply`, the three-block bend of
  `sflow_2d_mb_sed_open_faces_noneq` with `S 76 -1 2 -1 1` for `S 71 2.0`: the bed differs by 2e-6 of its l2
  norm, the slope flux at x = 2 m).  The sand slide (`S 90`) and the bed filter (`S 100`) do not look at the boxes.

Not with: wave generation on a block (`B 98`), hydrographs (`B 60 2-4`), the potential flow start (`I 11`
treats the open faces of the blocks as walls).

## Tests

Binaries from the environment: `REEF3D_BIN`, `DIVEMESH_BIN`.  Runs go to `./runs`.

| script | what | result |
|---|---|---|
| `t1_open_faces.py [steps]` | river reach (skewed grid turned by 30 deg, a mound, friction): one block with the native in- and outflow against `in` (block 1 the first 16 cells, `G 54 1 1 Q`), `out` (`G 54 1 2 h` on the last 16 cells), `both` (inflow block 1, middle block 0, outflow block 2); `SOLVER`, `NP`, `LAYOUTS`, `PAIRS`, `XCTRL` | SFLOW: 200 steps eta 3e-16 (`in`), 2e-13 (`out`, `both`); NHFLOW: `in` = `out` = `both` to 2e-11, all of them = one block to round-off at 20 steps (2e-16); 2 ranks = 1 rank.  The flow over the mound amplifies round-off differences (1e-10 after 100 steps, 1e-5 after 400; the same between one block on 1 and 2 ranks) |
| `nidelva_harbour_grid.py geo.dat out` | the two-block (`--split 2`: three-block) grid of the river and the harbour from a geo.dat raster | |
| `t2_harbour.py geo.dat rundir` | the Nidelva - harbour runs (`SOLVER`, `NP`, `TEND`, `OUT`, `Q`, `WL`, `AMR`, `SED`, `MF`; P8b: `SEDSTART`, `SPIRAL`, `SPIRALLAG`, `SUSP`, `SUSPCA`, `NONEQ`, `SUPPLY`) | see P8b below |
| `plot_harbour.py run [run ...] --out prefix` | speed and water level maps, the profile along the river, water volume and inflow | |
| `plot_bed.py run out.png [vmax]` | bed change map of a sediment run | |
| `plot_sed.py run [run ...] --labels .. --out prefix` | P8b: along the river the speed, the mean bed change of the cross-sections and the bed change of the right minus the left bank; the bed volume budget of the held reach, the free river and the harbour blocks (`prefix_budget.txt`) | |
| `nidelva_blocks_json.py griddir casedir` | the grid of `nidelva_harbour_grid.py` as the block topology of LAGOON (`blocks.json`: 10 corners, 12 edges, 3 blocks with their open faces) and a case folder (control.txt, ctrl.txt, the schematic bed as geo.dat); LAGOON's Set-up -> Blocks -> Make grid writes `blocks.dat`, `R 20 1`, `G 50 1`, `G 54` | LAGOON: no rule broken, 66,450 cells, smallest angle 68 deg; DIVEMesh and 60 s SFLOW on 2 ranks |
| `make_regression_cases.py` | `sflow_2d_mb_open_faces`, `sflow_2d_mb_open_faces_amr`, `nhflow_3d_mb_open_faces_mpi2` (P8b: `sflow_2d_mb_sed_supply` made by hand from `sflow_2d_mb_sed_open_faces_noneq`) | |

## The Nidelva - harbour grid

From the 5 m raster of the city (LAGOON map, EPSG:25832, `georef.json`; local frame): water `z <= 0`, land
`z > 0`; land specks in the water (single raster points) and isolated ponds are removed.

- Block 0 (native grid): a curvilinear corridor along the Nidelva from the upstream end of the raster to the line
  y = -300 m, where the river enters the harbour block: centreline by a least-cost path along the middle of the
  water, smoothed; cross-sections along the normals to the banks plus one dry cell (limited to 0.7 of the bend
  radius on the inner side), the last 60 m straight to the north onto the nodes of the harbour block; Laplace
  smoothing of the interior.  4.35 km, 435 x 10 cells.  Inflow: `B 60 1`, `W 10 Q` (DIVEMesh `C 11 1`).
- Blocks 1, 2: the harbour and the fjord, x -400 .. 2300 m, y -300 .. 2000 m, 10 m cells, i to the north; the
  river corridor ends on the south face of block 1 (point matching), the north faces (in the fjord) carry the
  water level (`G 54 1 2 WL`, `G 54 2 2 WL`).  Two blocks side by side (east, west) for the load of two ranks.
- The bed is SCHEMATIC: the raster has no water depths.  River 4 m in the middle (from the distance to the bank),
  harbour and fjord 8 m at the quays and deeper offshore (up to 25 m), smoothed over 15 m; the datum is moved up
  by the still water level (10 m).  Replace by surveyed depths when available (`bed` in `nidelva_harbour_grid.py`).
- Bed friction (P6c): `t2_harbour.py` and `nidelva_blocks_json.py` set `A 218 1` (SFLOW, Manning n = ks^(1/6)/20 with
  `B 50 0.05`: n 0.030) and `A 519 1` (NHFLOW, the log law of the bottom layer).  The P5d runs of this case had no
  bed friction (REEF3D's defaults `A 218 0`, `A 519 0`; `B 50` alone acts only on turbulence and sediment): their
  river velocities and the SFLOW/NHFLOW agreement are those of a frictionless river; SFLOW rerun in P8b (below).

## P8b: the Nidelva with bed friction and sediment (SFLOW, 2 ranks, `Q=100 WL=10 TEND=3600`)

Flow (`A 218 1`, n 0.030; `p8b_flow_*.png`): the long wave of the discharge reaches the harbour after 1300 s
(until then all 100 m3/s go into storage); after that the water volume still grows by 20 m3/s (0.15 % of the
volume over the last 2300 s): not quite steady at 3600 s.  The inflow discharge is 100.0 m3/s from 300 s on.  The
water level falls by 0.57 m along the 4.3 km of the river (0.72 m at the inflow, 0.15 m at the mouth); speeds
0.2 - 0.6 m/s in most of the river, 1.0 - 1.9 m/s at the narrow sections of the schematic bed (190 m, 350 m and
2300 m from the inflow) and 1.3 m/s at the mouth.  Without friction (`XCTRL="A 218 0"`, as the P5d runs) the
river lies 0.10 - 0.15 m lower (0.58 m at the inflow, 0.05 m at the mouth: about the same fall, 0.53 m, the losses
of the narrow sections and the flow around the schematic bed), speeds about 10 % higher (median 0.44 against
0.41 m/s, largest 2.1 against 1.7 m/s), and the water volume is steady at 3600 s (-1.7 against +19 m3/s); the
harbour stays at the level of the fjord faces in both (eta 1e-4 m).

Sediment (`SED=1`, d50 0.5 mm, MPM, MF 50, the bed moving from 1200 s: 33 h of morphology; `S 86 1`; supply
`S 76 850 1000 -2200 -2075`, the first 76 m of the corridor in full; `p8b_sed_profile.png`, `p8b_sed_budget.txt`,
`p8b_bed_sed*.png`):

| run | river (m3: net, erosion, deposition) | harbour blocks | wall time (s/step, 3600 s) |
|---|---|---|---|
| flow only | | | 0.148 |
| `S 95 1` | -704, -8736, +8033 | +838, -1192, +2030 | 0.188 |
| `S 95 2` (`S 97 1`) | -704, -8839, +8135 | +841, -1170, +2011 | 0.218 |
| `S 95 2` + `S 64 0.5` (`S 65 2.3`) | -637, -8381, +7744 | +773, -1103, +1876 | 0.250 |
| `S 95 1`, the held reach to 200 m (`SUPPLY="760 1000 -2200 -1900"`) | +291, -6330, +6622 | +829, -1178, +2007 | |

- The bed of the held reach does not move (0 to round-off) and the reach feeds the river: net supply 134 m3 (bed
  volume, `S 95 1`).  A held reach that ends at a constriction supplies the capacity of the constriction: the
  200 m box ends at the 1.9 m/s of the narrow section at 190 m and puts 1120 m3 into the river (0.65 m of deposit
  after its end, with `S 64` more: the lag carries the suspended part of the high capacity past the end).  The
  default of `t2_harbour.py` ends the box in uniform flow (0.55 m/s).
- The bed changes where the flow changes: scour in the narrow sections (0.6 - 0.7 m) and deposits (up to 1.8 m)
  after them, on the inner side of the bends that follow (`p8b_bed_sed2.png`: the two constrictions at x 640 and
  800 m); the S-bend at 2300 m (0.36 m deposit, 0.4 m between the banks); the mouth (-0.3 m) and a deposit in the
  river channel of the harbour block (y -295 .. 485 m).  The bends of the lower river carry little load at
  0.2 - 0.5 m/s (Shields number near the threshold): the bed of their two banks differs by 0.25 m at most.
- Nothing reaches the open faces in the fjord: no bed change north of y 1800 m (below 0.01 m north of 1000 m),
  speeds at the faces 0.04 m/s.  The sediment outflow through `G 54` (P8a) is not exercised by this case at the
  mean discharge.
- `S 95 2` against `S 95 1`: the same budget; the bed differs by 0.05 m rms (bed change 0.14 m rms) over the
  corridor, up to 1.5 m in the deposit after the upstream constriction where the spiral flow builds up over
  `L_a` (up to 79 m in the river, 212 m in the harbour) instead of at once; beyond 400 m from the inflow 0.009 m rms
  (change 0.06 m rms), at most 0.16 m.  Adding `S 64 0.5`: 0.03 m rms, 0.008 m rms beyond 400 m; `L_s` up to 23 m
  (w_s 0.072 m/s).  The lags across the corridor and the two harbour blocks: 2 - 15 passes per step, mostly 3.
- NHFLOW (`A 519 1`) was not rerun: 0.95 s per step with 5 layers on 2 cores, 6 - 7 h for the 3600 s.

## Refinement at the shoreline (lake at rest, `Q=0 AMR=1 TEND=300`)

Before the shoreline fix of SFLOW AMR (`../sflow_curvilinear/t8_shoreline_amr.py`) the refined box started
waves over the whole harbour (99 % of the speeds in the box 3.6e-2 m/s, the largest 0.64 m/s) and kept a
0.7 m/s eddy at the quay where the west edge of the box crosses the shore (x 1100, y 780; also in the runs with
the river, `AMR=1`).  Now 5.2e-3 m/s (99 %) and 0.20 m/s, quieter than one level (1.8e-2 m/s, 0.57 m/s); the eddy
is gone (0.002 m/s with the river).  The largest speeds left lie in single cells at steep quays of the
schematic bed, with and without refinement: the wetting and drying of one level.
