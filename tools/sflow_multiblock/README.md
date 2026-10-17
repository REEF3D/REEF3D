# SFLOW on multi-block curvilinear grids (G 50 1, DIVEMesh R 20 1): tests

Phase P2 of the multi-block curvilinear concept: several structured curvilinear blocks that meet
along their faces form the level-0 grid of REEF3D::SFLOW, with the patches of the mesh refinement
(REEFAMR) inside the blocks and across their interfaces.

## Running a case

`blocks.dat` in the case folder: the number of blocks, then per block a line `ni nj` and the
`(ni+1)(nj+1)` nodes `x y zb`, i fastest (as `curvgrid.dat`).  Every block must be right-handed
(i to j counter-clockwise).  Blocks touch where their face nodes coincide; DIVEMesh finds the
interfaces as runs of coincident face nodes (any of the eight index orientations) and writes them
to `DIVEMesh_Log/multiblock.log`, the blocks to `DIVEMesh_Paraview/multiblock-<b>.vtk`.

DIVEMesh `control.txt`: `R 20 1` (multi-block grid), `R 21 <ring>` (default 8), the usual `B 1`,
`B 10` and `M 10` (they are overwritten by the box of block 0).  Block 0 becomes the grid of the
run (the native grid of REEF3D): its own cells plus a ring of `R 21` cells on every face with an
interface, whose nodes are the nodes of the neighbour blocks.  The interfaces of block 0 have to
cover whole faces.  The blocks, their interface segments with the index transforms and the nodes
of blocks 1.. go into the section `BLKS` of `DIVEMesh_Grid/grid-geometry.dat`.

REEF3D `ctrl.txt`: `G 50 1`; nothing else marks the multi-block grid.  The ring must hold at least
`EXT + 3` cells (5 for the hydrostatic SFLOW), every neighbour block must be as deep behind its
face; REEF3D stops with a message otherwise.

Available: hydrostatic SFLOW (`A 220 0`) with what P1 offers on curvilinear grids (WENO or limited
reconstruction, RK2/RK3, Manning friction, wetting and drying, walls, `F 72` boxes in physical
coordinates), several ranks, mesh refinement (`G 1`, static boxes `G 10` and no-refinement boxes
`G 11` in physical coordinates, regridding `G 2` with `G 20` surface jumps and `G 22` shoreline,
placement `G 40 1`).  Sediment transport (`S 10 1`, bed load and the sand slide, phase P5a): `../sediment_curvilinear/README.md`.  Refused with a message: `A 220 >= 1` (non-hydrostatic pressure and Boussinesq
across blocks: phase P2c), subcycling `G 7 1`, moving bodies `X 10`.  Faces of the blocks 1.. without
an interface are walls; in- and outflow boundaries, solids and topography of the DIVEMesh grid
belong to block 0.  Discharge inflow and water level on the -i / +i faces of the blocks 1.. (`G 54`, phase
P5d): `../harbour_multiblock/README.md`.

Output: the native grid (block 0 with its ring) as before; the blocks 1.. and the patches in the AMR
output (`REEF3D_SFLOW_AMR`, VTS: level 0 for the blocks).  The AMR gauges (`P 51`,
`REEF3D_SFLOW_AMR_gauges.dat`) find the cell by its quadrilateral on curvilinear grids (also on one
block); the native gauges and probes still use the logical box.

## How it works

- Blocks 1.. are level-0 patches of REEFAMR (`lev[0]`, one per block), placed on the ranks for the
  load (largest block onto the rank with the fewest cells); block 0 is the native grid, distributed
  as before.  Level-0 index space: the native box, then the blocks side by side with gaps of 64
  cells (nothing is allocated there), each block starting on a multiple of 8 (tile alignment).
- A cell beyond a block face is found through the interface transform of its segment, on level l
  acting on the 2^l x 2^l children of the level-0 cells; corners through two faces (a regular vertex
  of four blocks).  Velocities are Cartesian components, so the filled values need no rotation.
- Before every stage the cells around the blocks and the ring of the native grid are filled from
  the stage input of the blocks that hold them (`fill_run0`); after the step, after a regrid and
  before the tags the ring and the cells around the blocks are filled again (the outer ring cells
  drift from the wall of the native box during a step, the CFL number of level 0 runs over the ring;
  a regrid needs the restricted values next to the interfaces).
- Patches: the tile maps cover the block index space; marked tiles are merged into rectangles and
  cut at the block faces (and at the rank boxes of block 0 without `G 40 1`).  A marked tile around a
  block (in the ring of block 0, in the gap of the others) marks the tiles of the cells it stands
  for, so the refined region continues across the interfaces with the buffer and the nesting; a
  tile around a block stays only where all of those are marked.
- The cells around a patch are filled from the cell they stand for (same level, or prolonged from
  its parent in the neighbour block).  Flux matching across an interface: the coarse face in the
  neighbour block, with the sign of the fluxes where its normal points against the fine faces
  (the face depth keeps its sign).  The fine bed uses the beds and node geometry of the block that
  holds the cell; native beds beyond the rank arrays come from the nodes.
- Singular vertices (three or five blocks around a node): the corner of the ring has no cells
  there; degenerate ring cells become solid, the other corner cells are taken through one of the
  faces, and no patch comes nearer than `max(4, (EXT+3+1)/2 + nest + 1)` level-0 cells.

Limitation: at a singular vertex level 0 conserves the volume only to about 1e-10 (relative, after
60 steps of a dam break over the vertex): the reconstruction next to the vertex reaches the corner
cell beyond both faces, which does not exist there and is taken through one face on each side.

## Tests (`REEF3D_BIN`, `DIVEMESH_BIN` select the binaries; runs go to `./runs`)

| script | checks |
|---|---|
| `t1_mb_identity.py [steps] [np]` | dam break across the interface of two blocks, the second with its indices turned by 90 deg: both blocks equal one block per cell to round-off (1e-15), 1 and 2 ranks |
| `t2_mb_quadrants.py [steps] [np]` | four blocks around a regular vertex, turned by 0, 90, 180, 270 deg, wavy skewed grid, dam break over the vertex: 2e-15 after 40 steps on 1-4 ranks, volume to round-off (`SKEW=0` Cartesian nodes, `RW` ring width) |
| `t3_mb_wetdry.py` | as t1 with a slope running dry (600 dry cells): 3e-15 on 1 and 3 ranks |
| `t4_mb_amr.py [steps] [np]` | t2 with a static refinement box over the vertex (`LEV` levels, `PLACE` for `G 40`): level 0 and every fine cell equal one block with the same box to 3e-14 after 40 steps, the same on 1-4 ranks, placed or not (bit for bit between rank counts) |
| `t5_mb_regrid.py [steps] [np]` | t4 with patches following the surface jumps (`G 20`), regrid every 4 steps, 60 steps: the same layouts, 5e-14 |
| `t6_mb_singular.py [steps] [np]` | a hexagon split into three blocks (singular vertex), dam break over it, with and without two levels: no patch within 3.6 m of the vertex, no spikes, volume to 1e-10 (see the limitation) |
| `t7_mb_shore.py [steps] [np]` | wetting and drying with shoreline (`G 22`) and surface-jump refinement across a turned interface, 2 levels, regrid every 4 steps: 7e-14 against one block, 1 and 3 ranks |
| `make_regression_cases.py` | writes the `Regression/cases/sflow_2d_mb_*` cases (with their `blocks.dat`) |

The identity tests need tile-aligned blocks (sizes and origins in multiples of the tile width
over 2^l in the one-block grid) to give the same refined region in both runs; the multi-block
runs are correct without that, the patches then only cover slightly different cells.
