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

Not with: wave generation on a block (`B 98`), hydrographs (`B 60 2-4`), the potential flow start (`I 11`
treats the open faces of the blocks as walls).

## Tests

Binaries from the environment: `REEF3D_BIN`, `DIVEMESH_BIN`.  Runs go to `./runs`.

| script | what | result |
|---|---|---|
| `t1_open_faces.py [steps]` | river reach (skewed grid turned by 30 deg, a mound, friction): one block with the native in- and outflow against `in` (block 1 the first 16 cells, `G 54 1 1 Q`), `out` (`G 54 1 2 h` on the last 16 cells), `both` (inflow block 1, middle block 0, outflow block 2); `SOLVER`, `NP`, `LAYOUTS`, `PAIRS`, `XCTRL` | SFLOW: 200 steps eta 3e-16 (`in`), 2e-13 (`out`, `both`); NHFLOW: `in` = `out` = `both` to 2e-11, all of them = one block to round-off at 20 steps (2e-16); 2 ranks = 1 rank.  The flow over the mound amplifies round-off differences (1e-10 after 100 steps, 1e-5 after 400; the same between one block on 1 and 2 ranks) |
| `nidelva_harbour_grid.py geo.dat out` | the two-block (`--split 2`: three-block) grid of the river and the harbour from a geo.dat raster | |
| `t2_harbour.py geo.dat rundir` | the Nidelva - harbour runs (`SOLVER`, `NP`, `TEND`, `OUT`, `Q`, `WL`, `AMR`, `SED`, `MF`) | |
| `plot_harbour.py run [run ...] --out prefix` | speed and water level maps, the profile along the river, water volume and inflow | |
| `nidelva_blocks_json.py griddir casedir` | the grid of `nidelva_harbour_grid.py` as the block topology of LAGOON (`blocks.json`: 10 corners, 12 edges, 3 blocks with their open faces) and a case folder (control.txt, ctrl.txt, the schematic bed as geo.dat); LAGOON's Set-up -> Blocks -> Make grid writes `blocks.dat`, `R 20 1`, `G 50 1`, `G 54` | LAGOON: no rule broken, 66,450 cells, smallest angle 68 deg; DIVEMesh and 60 s SFLOW on 2 ranks |
| `make_regression_cases.py` | `sflow_2d_mb_open_faces`, `sflow_2d_mb_open_faces_amr`, `nhflow_3d_mb_open_faces_mpi2` | |

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

## Refinement at the shoreline (lake at rest, `Q=0 AMR=1 TEND=300`)

Before the shoreline fix of SFLOW AMR (`../sflow_curvilinear/t8_shoreline_amr.py`) the refined box started
waves over the whole harbour (99 % of the speeds in the box 3.6e-2 m/s, the largest 0.64 m/s) and kept a
0.7 m/s eddy at the quay where the west edge of the box crosses the shore (x 1100, y 780; also in the runs with
the river, `AMR=1`).  Now 5.2e-3 m/s (99 %) and 0.20 m/s, quieter than one level (1.8e-2 m/s, 0.57 m/s); the eddy
is gone (0.002 m/s with the river).  The largest speeds left lie in single cells at steep quays of the
schematic bed, with and without refinement: the wetting and drying of one level.
