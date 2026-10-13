# NHFLOW on curvilinear grids (G 50): tests

Phase P3 of the multi-block curvilinear concept: REEF3D::NHFLOW on a boundary-fitted curvilinear
grid in the horizontal plane, the sigma grid in the vertical, with mesh refinement (REEFAMR).
The grid is set up as for SFLOW (`../sflow_curvilinear/README.md`): `G 50 1` with a DIVEMesh
river corridor, or `G 50 2` with `curvgrid.dat` (nodes `x y zb`, i fastest) and a DIVEMesh grid of
the logical box (`B 2 ni nj nz`, `B 10 0 ni*h 0 nj*h 0 h`; the vertical as usual). The cell bed is
the mean of its four node beds.

## What runs on curvilinear grids

- hydrostatic and non-hydrostatic NHFLOW (`A 520 0`, `1`, `2`), HLL (`A 511 1`), all reconstructions
  (`A 514`), RK2 and RK3, wetting and drying, walls, `F 72` boxes in physical coordinates (cell
  centres), several ranks
- boundaries: river inflow (`B 60`, discharge along the inflow face normal, `B 61` plain or log
  profile) and stage outflow (`F 62`); tides and currents at open edges (`B 510`-`B 525`, Riemann,
  Flather, clamped level or discharge, waves on the background); wave generation by relaxation
  (`B 98 2`, `B 99 1/2`, zones by `B 107` / `B 108` lines in physical coordinates, the wave frame
  by `B 105`), Dirichlet (`B 98 3`) and active generation (`B 98 4`), active absorption (`B 99 3-5`)
- diffusion (`A 512 1` explicit, `A 512 2` implicit) and turbulence (`A 560 1/21` k-epsilon,
  `2/22` k-omega, `31` LES), bed friction (`A 519 1`)
- mesh refinement (`G 1`): static boxes (`G 10`, `G 11`, physical coordinates), regridding with the
  surface flags (`G 2`, `G 20`, `G 21`), vertical refinement (`G 6 1`), subcycling (`G 7 1`),
  wetting and drying in the patches (`G 30 1`), the composite pressure and implicit diffusion,
  placement (`G 40`), AMR output as `.vts` on the physical nodes, AMR gauges (`P 51`) by the cell
  whose quadrilateral holds the point

- the potential flow initialisation (`I 11`): the compact Laplacian, the inflow flux along the
  normal of the inflow faces, the velocities by the chain rule (phase P4)
- multi-block grids (`G 50 1` with DIVEMesh `R 20 1`): `../nhflow_multiblock/README.md` (phase P4)
- `G 52 1`: the i-walls act as the j-walls (the velocity pair mirrored at the face normal, the
  pressure correction Neumann, the bed smoothing `A 509` from the old values), so that a wall acts
  the same in any index direction; set on multi-block grids, by default the walls along i keep the
  Cartesian treatment (ghost velocities 0)

Refused with a message: HLLC (`A 511 2`), `A 520 3`, 2D slices (`j_dir 0`), side-wall friction
(`A 519 2`), 6DOF, membranes, sediment, FEM structures, solids (`A 580`-`599`), wind, SEASTATE,
the B 96 relaxation zones (an x range of the box: give `B 107` / `B 108`), Dirichlet or
active generation on a current (`W 10`), `B 97 1`, the LAGOON output of the mesh refinement
(`P 18`). The native probes and line probes still refer to the logical box.

## How it works

As SFLOW (`sflow_curv.h`): the state stays cell-centred with Cartesian velocity components; the
faces carry the normals of the geometry layer. The code is in `src/nhflow_curv.h/.cpp`, every
curvilinear path is a branch on `p->geo_curv` (the Cartesian kernels are not changed).

- HLL: after the reconstruction (both components on both face directions) and the wet-dry fluxes,
  Us/Un and Ve/Vw become the normal velocities; momentum flux of UH: UH (u.n) + P n_x (VH: n_y,
  WH: no pressure), continuity UH.n; divergence: sum of face flux times length over the area;
  bed-slope source g eta sum(d_face S_face)/A (well balanced on any grid)
- boundary faces (`bc_fluxes`, after `ghostcell::start1V/2V`): walls the ghost pressure along the
  face normal; in- and outflow faces and open edges the flux of the ghost state through the face
- side walls (j-faces): the ghost cells mirror the normal velocity of the wall face
  (`pair_bc`, called by `ghostcell::start4V` once U (UH) and V (VH) have their ghost cells); the
  river outflow ghost cells act on the normal velocity (positive outwards)
- sigma metrics: gradients of eta and the depth limited along the grid lines (`A 524 2`, the
  gradient that has the limited line derivatives as directional derivatives; identical to the
  Cartesian one on Cartesian nodes), central (`A 524 1`); the Laplacians compact
- Poisson rows: compact horizontal part (face length over area and the centre distance along the
  normal), sigma cross terms by the chain rule, the non-orthogonal part deferred with the pressure
  of the right-hand side; right-hand side divergence and pressure gradients by the chain rule
- diffusion: the stress form 2 nu u_xx + nu u_yy + nu v_xy written as nu Lap_h(u) +
  nu d/dx(div_h): the compact Laplacian (implicit with `A 512 2`) plus the rest explicit; scalar
  diffusion and the k, epsilon/omega convection on the face normals and lengths
- iowave: positions from `position::pos_x/pos_y` (the physical cell and face centres), the
  in- and outflow velocities along the face normals, the open edges decomposed into the normal and
  tangential parts of the edge face
- mesh refinement: the patch nodes subdivide the level-0 nodes bilinearly, so a coarse face is two
  collinear halves and the flux matching stays the mean of the fine faces; restriction and the
  conservative prolongation per area of the children (`cw4`); the fine bed is shifted so that its
  area mean is the parent's bed (a lake at rest stays at rest); refluxing (`G 7 1`) with face length
  over area

## Tests (`REEF3D_BIN`, `DIVEMESH_BIN` select the binaries; runs go to `./runs`)

| script | checks |
|---|---|
| `t1_nhf_rotation.py [steps]` | dam break in a closed basin (`CTRL` the ctrl lines): Cartesian, curvilinear on the same nodes, turned by 180 and by 90 deg: per cell index for `A 520 0`, `1`, `2`: on the nodes eta 4e-16, U 1.5e-15, P 4e-12; turned and shifted to (1000, 2000) m 2e-13 (the rounding of the coordinates; `OX=1 OY=1`: 2e-15) |
| `t2_nhf_lake.py [steps]` | lake at rest with an island (wetting and drying) on a wavy, skewed grid (`A520="0 2 1"`): 2e-5 velocities; the same lake submerged: exactly 0 (the residual comes from the NHFLOW wetting and drying, the Cartesian grid gives 4e-4) |
| `t3_nhf_channel.py` | gate test: linear waves (kh 2) along a channel, Cartesian, turned by 30 deg at UTM coordinates, orthogonal mapped (stretched, `map0`) and turned (`map30`); `WG=dirichlet` (`B 98 3`, `B 99 3`) or `WG=relax` (`B 98 2`, `B 99 1`, zones by `B 107/108`), `AMR=adapt` with patches: turned = Cartesian (map30 = map0) to 0.1 % of the wave height, Poisson iterations within 0.5 % (3 % with relaxation zones) |
| `t4_nhf_river.py [steps]` | river reach (`B 60 1`, `F 62`, bed friction) Cartesian, on the same nodes, turned by 30 deg; `VAR=ediff, idiff, ke, kw, les`; `XCTRL="A 514 0;A 524 1;"` (linear reconstruction): turned = on the nodes to 1e-10 (1e-8 with turbulence), the WENO weights give 1e-3 |
| `t5_nhf_tide.py [steps]` | oblique tide through a basin with Riemann and Flather edges, turned by 30 deg (tide direction and origin turned): 1e-13 with the linear reconstruction |
| `t6_nhf_mpi.py [steps]` | 2 and 4 ranks: curvilinear on the Cartesian nodes = Cartesian per rank (NH, 4e-12 in P); hydrostatic on the turned grid bit for bit against one rank |
| `t7_nhf_amr.py [steps]` | mesh refinement, dam break with patches following the surface: Cartesian, on the same nodes, turned by 90 deg; `VAR=sub` (G 7 1), `vref` (G 6 1), `two` (two levels), `static` (G 10), `place` (G 40 1, with `NP=2`), `cdiff` (composite implicit diffusion, A 512 2); `NP=2`: level 0 and every patch to 2e-15 (on the nodes) and 3e-13 (turned), the same layouts, volume to round-off, the AMR gauges turned = on the nodes; with `cdiff` on the nodes = Cartesian to truncation (the isotropic diffusion), turned = on the nodes to 3e-13 |
| `t8_nhf_amr_wavy.py [steps]` | mesh refinement on the wavy grid: lake at rest over a bump under six patches (`A 520 0`, `2`): velocities 2e-15 on both levels; dam break with regrid every step (G 7 0 / 1): volume to 6e-15 |
| `t9_nhf_amr_shore.py [steps]` | wetting and drying in the patches (`G 30 1`), water running up a plane beach under a static patch: turned by 90 deg = on the nodes to 1e-15 (both against Cartesian to the bed rounding of DIVEMesh, 1e-10) |
| `t10_nhf_potential.py` | potential flow initialisation (`I 11`) on a river reach (discharge inflow, stage outflow): curvilinear on the nodes = Cartesian bit for bit; with a bump, turned by 30 deg at UTM coordinates = on the nodes to the solver tolerance (4e-10) |
| `make_regression_cases.py` | writes the `Regression/cases/nhflow_3d_curv_*` cases (with their `curvgrid.dat`) |

Not identical to Cartesian on the same nodes, by design: the diffusion (the isotropic form above;
equal to truncation, the turned grid equal to the grid on the nodes), and the WENO reconstruction
on grids turned by angles other than multiples of 90 deg (the nonlinear weights act on the
Cartesian components).
