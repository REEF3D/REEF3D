# 6DOF bodies on curvilinear grids, mesh refinement and multi-block grids (P5c): tests

Phase P5c of the multi-block curvilinear concept: the floating and moving bodies of REEF3D::NHFLOW and
REEF3D::SFLOW on curvilinear grids (`G 50 2`, one block from `curvgrid.dat`), on the patches of the mesh
refinement (`G 1`, REEFAMR) and on multi-block grids (`G 50 1`, DIVEMesh `R 20 1`).

## Running a case

The body keys are those of the Cartesian runs; the hull (`X 110`, `X 133`, STL) and the positions
(`X 23`, `X 350` ship data, probes) are in physical coordinates.

| | one curvilinear block | mesh refinement (`G 1`) | multi-block (`G 50 1`) |
|---|---|---|---|
| NHFLOW floating `X 10 1` | yes | yes (`G 12` zone, `G 13` wake, `G 40 1` placement) | yes (not with `G 40`) |
| NHFLOW prescribed `X 10 2` | yes | yes | yes |
| NHFLOW ship pressure `X 10 3` | yes | as before (level 0) | refused |
| NHFLOW ship loads `X 350` | yes | yes | yes |
| SFLOW prescribed `X 10 2` | yes | yes | yes |
| SFLOW ship pressure `X 10 3` | yes | yes | yes |

Refused on curvilinear grids with a message: the level-set force integration `X 60 2` (use `X 60 1`, the
STL integration), the viscous forces `X 38`, the pierced-column water level `X 13`, the nets `X 320`.

## How it works

- Point location (`grid_locate.cpp`): a bucket table over the rank holds the quadrilaterals of the cell
  centres (`geo_quad`: the bilinear coordinates of a point, the interpolations `ccipol*`, `ccslipol4`,
  `posc_sig`), the cells as quadrilaterals of their nodes (`geo_cell`: the cell that holds a point,
  half-open), the cell centres (`geo_cells`: the cells near a triangle).  `geo_width` is the physical width
  of a cell (area over the mean face lengths), `geo_grad` the chain-rule gradient from the four neighbours,
  `cellarea` the area (bit for bit `DXN DYN` on Cartesian grids).
- The hull on the sigma grid (`geo_raycast_sigma.cpp`): vertical rays and band distances on the physical
  cell centres; the rank box from the nodes.  The normal of the forcing (`X 15 1`, `2`) is the physical
  gradient of the level set along the sigma level with the sigma correction.  Smoothing widths (`A 526`,
  `X 41`) and the band `DSM` from `geo_width`.
- Loads: the triangle belongs to the rank of the cell that holds its centroid (`geo_cell`), the pressure is
  interpolated bilinearly in the cell centres.  The surface integrals of the body (volume, centre, inertia)
  are taken about the centre of the hull box: at UTM coordinates about the origin they lose 9 digits.
  Remeshing (`X 185`): the target length from the mean physical cell width (rounded to 6 digits, so that
  the hull does not depend on the ranks), the hull remeshed about the centre of its box.
- Ship loads (`X 350`): the actuator cells by `geo_cell`, positions `pos_x`, `pos_y`.
- SFLOW (`6DOF_obj_ray2D_curv.cpp`): the 2D level set is the exact signed distance to the waterline (the
  section of the hull with the still water level): inside by the parity of the waterline segments to the
  right of the centre, distance to the nearest segment within 10 cell widths, the draft from vertical rays.
  The Cartesian rays along the grid lines and the re-initialisation do not apply on a curvilinear grid.
- Mesh refinement: the zone around the hull (`G 12`, `G 13`) flags the cells by their physical centres.
  A point of the hull belongs to the finest grid in whose logical range it lies: `lexer::curv_logical`
  gives the global logical coordinates of a physical point (bilinear inverse over the global nodes, the
  same on every rank); the patch cells subdivide the level-0 cells bilinearly, so the logical range of a
  patch is exactly its physical area.  The placed patches (`G 40 1`) take the triangles by these ranges
  (`reefamr::finest_rank`).  SFLOW patches build the waterline level set and the draft from the hull on
  their own cells (as level 0).
- Multi-block grids: the hull is cast on every block (the blocks 1.. as level-0 patches, the native grid
  with its ring); a triangle belongs to the rank of the native cell (not the ring) or of the block piece
  that holds its centroid, and every quadrature point of it is sampled on the grid of the rank that holds
  the point (a hull triangle reaches across the block faces).  SFLOW: the initial draft depression
  (`X 10 3`, `X 207 0`) and the first-step pressure of the blocks as on level 0.

## Tests

Binaries from the environment: `REEF3D_BIN`, `DIVEMESH_BIN` (as `../sediment_curvilinear`).  Runs go to
`./runs`.  Fixed time step (`N 48 0`, `N 49 0.01`): NHFLOW takes the adaptive step from the speed on
curvilinear grids.

| script | what | result |
|---|---|---|
| `t1_6dof_rotation.py [steps]` | one block: Cartesian, curvilinear on the same nodes (curv_id), turned by 90 deg about (1000, 2000) (rot90), turned by 30 deg at UTM coordinates or about (1000, 2000) (`ROT30=1`, `ROT30=small`); `CASE` floating, towed, thrust, propeller, shipwave, shipwave_stl (NHFLOW), towed, shipwave, shipwave_stl (SFLOW); `NP`, `YAW`, `LIN`, `ANG`, `FIXCG`, `XCTRL` | NHFLOW towed, thrust, propeller, shipwave_stl: curv_id = Cartesian to 2e-15, rot90 to 1e-11; floating: 1e-4 (the remeshing of the hull breaks ties differently; `X 185 0`: round-off); rot30 with `XCTRL "A 514 0;A 524 1;X 15 2"`: 1e-12 about (1000, 2000), 1e-9 at UTM (`X 185 0`); SFLOW: rot90 1e-12, rot30 (`A 211 0`) 1e-8 at UTM; curv_id vs Cartesian 1e-3 (the waterline level set is built differently); 2 ranks = 1 rank |
| `t2_6dof_amr.py [steps]` | the body of T1 with one refinement level around the hull (`G 12`), Cartesian, curv_id, rot90; per global cell of all grids; `CASE`, `SOLVER`, `NP`, `LEV`, `WAKE` (`G 13`), `PLACE` (`G 40 1`) | NHFLOW: the same patches, curv_id = Cartesian to 2e-15, rot90 to 4e-14 (`X 185 0`); placed patches on 2 ranks = 1 rank; SFLOW towed and shipwave_stl: rot90 to 4e-11 |
| `t3_6dof_mb.py [steps]` | the skewed wavy grid of P5a T5 as one block and as four quadrants turned by 0, 90, 180, 270 deg, the hull over the vertex; `CASE`, `SOLVER`, `NP`, `LEV`, `NOTURN`, `ROTONE`, `SHIFTX/Y`, `TOP` | NHFLOW, blocks not turned (`NOTURN`, `X 185 0`): 3e-9 after 20 steps, the loads of the first step to 6 digits; SFLOW `X 10 3` and `X 10 2` with `LEV 1`: four = one to 1e-15; SFLOW `X 10 3`, `LEV 0`: 4e-16 |
| `t4_canal_bend.py` | P5 gate: a floating box (0.5 x 0.5 m, draft 0.15 m; surge, sway, heave, yaw free) drifting with the current (0.4 m/s) through a canal of 3 m with a 90 deg bend (centreline radius 3 W): three blocks (inlet reach with `G 54 1 1`, bend, outlet reach native) against one block of the same nodes (native inflow); `RUN`, `TEND`, `NP`, `PNG` | the box crosses both block interfaces and the bend in 25 s and follows the axis (0.1 m to the inside in the bend, 0.4 m to the outside 9 m below it); three vs one: 1.1 mm after 79 s and 23 m, the heading 0.2 deg (P6b: 1 cm and 1.3 deg, from the wall lists of the blocks, see `../nhflow_multiblock/README.md`) (it keeps its heading, 0 -> 9 deg, while the axis turns by 90 deg); hydrostatic (`A 520 0`) the heave grew in the bend |
| `cmp_runs.py a npa b npb` | two runs cell by cell (level, block, I, J) | |

## Notes

- Rotation by angles other than 90 deg: the WENO weights (`A 514 4`, SFLOW `A 211 4`) act on the Cartesian
  components and the default direct forcing `X 15 1` weights every Cartesian component by `|n_x|`, `|n_y|`
  of the hull normal: neither is invariant under such a turn (rot30 2e-3 with the defaults).  `X 15 0`
  and `X 15 2` are vector forms.
- NHFLOW bodies and the orientation of the index frame: with the hull in a block whose indices are turned
  (multi-block, or one block with its index frame turned, `ROTONE`), the solution differs from the unturned
  frame by 1e-4 near the hull (eta 2e-4 after one step).  The cause is in the flux scheme at the faces next
  to and inside the body (`DF < 0` is treated as dry): in `nhflow_signal_speed.cpp` the x faces with a dry
  right neighbour take `Sn = us + cs` while the y faces take `Sw = ve + 2 ce`, and faces between two body
  cells get zero signal speeds, where `hll_sweep_x/_y` take the left flux.  Both give a preferred index
  direction (also on Cartesian grids: x against y, + against -).  Making them symmetric changes the
  Cartesian results of every case with a floating body or a dry bed, so it is left to a separate change
  (measured: 2e-4 -> 2.5e-5 with both changes; more such conventions remain).
- Multi-block grids: the band `DSM` and the remeshing length are the mean width of the native grid (block 0
  with its ring), so the remeshed hull of a multi-block run differs from the one-block run; compare with
  `X 185 0`.
- Large hull triangles (`X 185 0`) with placed patches (`G 40 1`) or several ranks on a multi-block grid:
  the quadrature points of a triangle far from its centroid lie outside the grids of the owning rank and
  are clamped to its outermost cells.  With the remeshed hull (default) the triangles are about one cell.
- SFLOW `X 10 2` on the blocks 1.. without refinement: the body of the previous stage (level 0 moves the
  body in its stage, after the blocks): 6e-4 against one block.  The refined patches do the same, so with
  the zone around the hull (`G 12`) four blocks = one block.
- SFLOW mesh refinement on a curvilinear grid: the patches build the waterline level set from the hull
  (exact distance); the Cartesian patches rebuild it from their draft (interface half way between the
  centres).  curv_id against Cartesian with refinement: 2e-2 near the hull.
