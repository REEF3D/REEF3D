# NHFLOW on multi-block curvilinear grids (G 50 1, DIVEMesh R 20 1): tests

Phase P4 of the multi-block curvilinear concept: several structured curvilinear blocks that meet along
their faces form the level-0 grid of REEF3D::NHFLOW, hydrostatic and non-hydrostatic, with the patches
of the mesh refinement (REEFAMR) inside the blocks and across their interfaces.

## Running a case

As for SFLOW (`../sflow_multiblock/README.md`): `blocks.dat` in the case folder (number of blocks, per
block `ni nj` and the nodes `x y zb`, i fastest; every block right-handed), DIVEMesh `control.txt` with
`R 20 1` (`R 21` the ring width, default 8), `B 1` and `B 10` (the vertical: `B 10` z range over `B 1`
gives the sigma layers), REEF3D `ctrl.txt` with `G 50 1`.  Block 0 becomes the native grid of the run
(its cells plus a ring of `R 21` cells on every face with an interface); the faces of block 0 with an
interface may touch several blocks, together they have to cover the whole face.  The faces of block 0
without an interface keep the DIVEMesh boundary types `C 11` (-i), `C 14` (+i), `C 13` (-j), `C 12`
(+j): in- and outflow, open edges and wave generation belong to block 0.  The faces of the blocks 1..
without an interface are walls.

Available: everything P3 offers on curvilinear grids (`../nhflow_curvilinear/README.md`) in every block:
hydrostatic and non-hydrostatic (`A 520 0`, `1`, `2`), HLL, all reconstructions, RK2/RK3, wetting and
drying, explicit and implicit diffusion (`A 512 1`, `2`: one composite solve over the blocks), bed
friction, `F 72` boxes in physical coordinates, several ranks; in- and outflow, tides and open edges,
wave generation and absorption on the faces of block 0; the potential flow initialisation (`I 11`) as
one solve over the blocks; mesh refinement (`G 1`) inside the blocks and across their interfaces (static
boxes `G 10`, `G 11`, regridding `G 2` with the surface flags `G 20`, `G 21`, vertical refinement
`G 6 1`, wetting and drying in the patches `G 30 1`, placement `G 40 1`).  Refused with a message:
subcycling (`G 7 1`), floating bodies (`X 10`), turbulence (`A 560`, not in the mesh refinement
either), breaking (`A 550`).  Output: the native grid (block 0 with its ring) as before, the blocks
1.. and the patches in the AMR output (`REEF3D_NHFLOW_AMR`, `.vts`, level 0 for the blocks); the AMR
gauges (`P 51`) find the cell by its quadrilateral in every block.

Keys: `G 51 n` passes of the block-local V-cycles in the preconditioner (default 2), `G 53 w` the
coarse space across the blocks, aggregates of `w x w` columns (default 4, 0 off), `G 52 1` the
isotropic walls (set on every multi-block grid, see below; give it to a one-block run to compare).

## How it works

- Blocks 1.. are level-0 patches of REEFAMR (`lev[0]`, `reefamr_block.cpp` of P2), each with its own
  lexer (the block nodes), `fdm_nhf` and NHFLOW objects; block 0 is the native grid.  Before every
  stage the cells around the blocks and the ring take the stage input of the cells they stand for
  (eta, water level, wet and deep, U, V, W, UH, VH, WH per layer, P per node; `fill_level0`), the bed
  of these cells once at the start (`fill_bed0`).  Velocities are Cartesian components: nothing turns.
  Each block then runs its stage as one grid; the ring is refilled at the end of the step.
- Composite pressure (`nhflow_amr_press.cpp`): the rows of every block, the ring no unknown (filled
  from the vectors, `pr_fill0`); BiCGStab with a preconditioner of three parts: the coarse correction
  across the blocks (`G 53`: aggregates of 4 x 4 columns of every block, all layers, the coarse
  operator summed from the rows of the blocks and factorized on every rank), then a REEFMG V-cycle per
  block on the residual (block 0 over the rank grid without the ring, the blocks 1.. patch-local), `G 51`
  passes with the neighbour corrections filled in between; the patches follow as before (FAC).  The
  same machinery solves the potential flow (`I 11`, `nhflow_amr::potential`).
- Walls along i and j: the Cartesian NHFLOW treats x-walls (ghost velocities 0, the pressure correction
  Dirichlet) and y-walls (slip, Neumann) differently, so a wall would act differently in a block whose
  indices are turned.  `G 52 1` (set on multi-block grids) makes the i-walls act as the j-walls: the
  velocity pair mirrored at the face normal (`nhflow_curv::pair_bc`), the other fields symmetric, the
  pressure correction Neumann; the bed smoothing of `A 509` from the old values (independent of the
  sweep direction); the porosity of the partition halo below the bed and above the surface set as for
  the own columns, so that a block split over ranks along i or j acts as one rank.  On Cartesian runs
  without `G 52` none of this changes (the regression cases stay bit for bit).
- Ghost nodes beyond a wall next to an interface are taken through the interface first and reflected at
  the wall of the neighbour block (`reefamr::block_node_any`), as the one-block grid reflects them.
- Mesh refinement: as for SFLOW (P2); the coarse faces next to fine patches across a turned interface
  take the fine fluxes with the sign of their normal (`flux_hook`, `reefamr_match::sgn`).

## Tests (`REEF3D_BIN`, `DIVEMESH_BIN` select the binaries; runs go to `./runs`)

| script | checks |
|---|---|
| `t1_mb_identity.py [steps] [np]` | dam break across the interface of two blocks, the second turned by 90 deg: both blocks equal one block (`G 52 1`) per cell, all fields, `A520="0 1 2"`: 1e-15 in U, 5e-12 in P (|P| 200), 1 and 2 ranks |
| `t2_mb_quadrants.py [steps] [np]` | four blocks around a regular vertex turned by 0/90/180/270 deg, wavy skewed grid, bumpy bed, dam break over the vertex: one block to round-off for `A 520 0, 1, 2` (`A 509 0`), volume to round-off |
| `t3_mb_count.py [steps] [np]` | 1, 4, 9, 16 blocks (NB x NB, each turned): the same solution as one block (to the solver tolerance), pressure iterations per solve against the block count |
| `t4_mb_potential.py [np]` | potential flow (`I 11`): channel (block 0, inflow and outflow) blocked by a dry obstacle, the bypass above it split into 1, 3, 8, 15 turned blocks: initial velocities against one block (to the solver tolerance), iterations against the block count |
| `t5_mb_amr.py [steps] [np]` | the four blocks of t2 with a static refinement box over the vertex (`LEV` levels) or patches following the surface (`ADAPT=1`): level 0 and every fine cell against one block with the same refinement |
| `make_regression_cases.py` | writes the `Regression/cases/nhflow_3d_mb_*` cases (with their `blocks.dat`) |

`p4lib.py`: `run4` (DIVEMesh and REEF3D in a fresh folder, `REEF3D_AMR_DUMP` set: `nhflow_amr` writes the
state of every grid of a rank, own cells, to `amr/amr_r<rank>.bin` at every output step), `amr_grids`,
`level0`.
