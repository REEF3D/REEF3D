# SFLOW on curvilinear grids (G 50): tests

REEF3D::SFLOW runs on a boundary-fitted curvilinear grid in the horizontal plane (phase P1 of the
multi-block curvilinear concept). The grid is a logical box of `ni x nj` cells whose nodes have
arbitrary (x, y) positions.

## Running a case

One DIVEMesh run (river corridor, `R 1 2`): the Cartesian setup (geodat, solids, `C 11`-`C 14`
with the inflow and outflow sides) plus `R 1 2`, `R 2 <design water level>`, `R 4 0 <nj>`.
DIVEMesh builds the corridor grid, then writes the logical box of it (`ni x nj` cells of the mean
spacing, the vertical of the Cartesian grid) as the grid of the run: inflow at i = 0, outflow at
i = ni, banks at j = 0 and nj, the CURV section in `DIVEMesh_Grid/grid-geometry.dat`.
`ctrl.txt`: `G 50 1`. Cell bed levels are the mean of the four node bed levels.

Grids from elsewhere (`G 50 2`): `curvgrid.dat` in the case folder, first line `ni nj`, then
`x y zb` of node `(i,j)` for `j = 0..nj`, `i = 0..ni` (i fastest), and a DIVEMesh grid of the
logical box: `B 2 ni nj 1`, `B 10 0 ni*h 0 nj*h 0 h` with `h` about the mean spacing, `C 11 1`
(inflow at i = 0), `C 14 2` (outflow at i = ni), `C 12 21`, `C 13 21` (banks). DIVEMesh `R 1 1`
writes `curvgrid.dat` and these keys (`DIVEMesh_Log/corridor.log`) without changing its grid.

Available on curvilinear grids so far: hydrostatic SFLOW (`A 220 0`), non-hydrostatic pressure
(`A 220 1`, `2`, `3`), Boussinesq (`A 220 4`), WENO or limited reconstruction (`A 211`), RK2/RK3,
Manning friction (`A 218`), wetting and drying, walls, discharge inflow (`B 60`) along the inflow
face normal, stage outflow (`F 50`, `F 62`), `F 72` boxes in physical coordinates, several ranks,
mesh refinement (`G 1`: static boxes `G 10` and no-refinement boxes `G 11` in physical
coordinates (cell centres), regridding `G 2`, subcycling `G 7 1`, placement `G 40 1`, with the
composite non-hydrostatic and Boussinesq solves), VTP output on the physical nodes, AMR output as
VTS (structured grids on the physical nodes), AMR gauges (`P 51` with `G 1`) at physical points,
explicit horizontal diffusion (`A 212 1`, `W 2`) with the algebraic eddy viscosity (`A 260 3` mixing
length, `A 260 4` parabolic; P6a), also on the patches and the blocks. Refused with a message:
`A 220 >= 5`, the implicit diffusion `A 212 2`, the turbulence transport models `A 260 1`, `2`,
wave generation, patch boundary conditions, 6DOF, SEASTATE, `I 11`. Sediment transport (`S 10 1`, phase P5a): `../sediment_curvilinear/README.md`. The native probes and other output that use x/y positions still refer to the logical box;
the LAGOON print of the AMR hierarchy is not written on curvilinear grids.

Numerics (P1b, P1c):
- non-hydrostatic pressure: velocity correction and divergence by central differences with the
  chain rule through the cell centres; Poisson operator compact (face length over area times the
  centre distance along the face normal), the non-orthogonal part of the face fluxes as a deferred
  correction from the previous pressure (exact for a Cartesian or orthogonal grid)
- Boussinesq: the inversion V = H (u_a + V1') runs on the grid-aligned components a = u_a.e1,
  b = u_a.e2 (e1 through the i-neighbour centres), per component with the other one lagged as on
  the Cartesian grid, but the own component with all its terms, also the non-orthogonal ones:
  BiCGStab with the matrix-free operator and the line solves of the normal-derivative rows as
  preconditioner (a few tens of iterations).  With the rows alone and the non-orthogonal own terms
  lagged, the scheme went unstable on skewed grids; both components as one system converge badly
  for h/dx > 2 (grad div).  Stable at 45 deg non-orthogonality; on the Cartesian grid the same as
  the line solver to the solver tolerance (`N 44`)
- mesh refinement: restriction and the conservative fix-ups area weighted, refluxing with the face
  lengths over the cell areas (bilinear subdivision: the child faces are exact halves of the parent
  face, the flux matching stays exact); the composite NH pressure with the curvilinear rows
- mesh refinement on curvilinear grids: the fine bed is prolonged with area-weighted children
  means (a lake at rest stays at rest on skewed grids); the composite Boussinesq solve exchanges
  the vector a e1 + b e2 between the grids (the bases of parent and children differ)
- mesh refinement at a shoreline (the coarse grid does not resolve it: a coarse cell is dry while some of
  its children lie below the water, or wet with a mean water level plus mean bed above the water):
  - the initial hierarchy fills the children of a dry parent whose bed lies above the surface of its wet
    neighbours up to that surface (as the initial state of level 0 is a water surface), a shallow parent
    with its surface; the patches of a regrid and the ghost cells keep the mass of the parent (dry stays dry)
  - a covered cell at a shoreline shows the uncovered neighbours (reconstruction, pre-balanced bed term,
    the slopes of the prolongation) the surface of its wet children, not its mean water level plus mean
    bed; with subcycling (`G 7 1`) through the steps of the coarser level (`sfc`, an offset carried with
    the water level)
  - the fine faces on a patch boundary that replace a coarse face (refluxing): no mass flux into a dry
    cell of the patch from a wet ghost cell (as `aij_E` does between two cells of a grid), and a face
    between two dry cells whose ghost cell has a wet parent is a wall of the parent (the pressure of the
    parent's surface; the flux between the two dry surfaces, bed + film, moved the water at rest)
  - before: 0.1 m/s at a lake at rest with an emerged island where a patch boundary crossed or ran near
    the shoreline (0.7 m/s at a steep quay of the Nidelva - harbour run); now the same as on one level
    (`t8_shoreline_amr.py`)
- horizontal diffusion (P6a, `sflow_ediff::curv_stress`): the depth-integrated stress nu h (grad u + grad u^T)
  through the four faces of a cell over its area (finite volume, conservative), nu h of a face the mean of
  its two cells; the normal derivative from the two centres (as the compact Laplacian of the NH pressure),
  the non-orthogonal part k.grad u and the transposed part from the face mean of the cell gradients; no
  stress through walls and towards dry cells (free slip, as the dry neighbours of the Cartesian operator).
  The Cartesian operator keeps its form h nu (2 u_xx + u_yy + v_xy); on the same nodes the two agree to
  second order (the effect of the diffusion on a small dam break within 6 %, `t9b`), not to round-off; with
  large depth gradients the conservative form differs more (the gradient of h nu)
- eddy viscosity on curvilinear grids and on the grids of REEFAMR (`sflow_ediff::algebraic_eddyv`): every
  grid computes `A 260 3` / `4` from the state of the stage in its cells and its first ghost layer (both
  cells of a face of two grids see the same viscosity); the native turbulence module stays void there
- `A 209 0` is recommended (the default in-place bed smoothing depends on the partition)

Limitation: Boussinesq with mesh refinement on a skewed grid went unstable when the finest cells
were much smaller than the depth (dx < h/4 at 20 deg non-orthogonality; dx = h/4 was stable): the
coarse-fine coupling of the composite operator is not symmetric, and the lagged other component
then amplifies.  On orthogonal grids (and single-level skewed grids) no such limit was found.

Several blocks (multi-block grid, DIVEMesh `R 20 1`): see `../sflow_multiblock`.

## Tests (`REEF3D_BIN`, `DIVEMESH_BIN` select the binaries; runs go to `./runs`)

| script | checks |
|---|---|
| `t9_diffusion_rotation.py` | horizontal diffusion (`A 212 1`, `W 2`), mixing length (`A 260 3`), parabolic (`A 260 4`): river reach with an emerged island on the logical grid and rotated by 30 deg (`A 211 0`): eta per cell index to 1.0e-10 (as without diffusion); a skewed, wavy grid runs |
| `t9b_diffusion_cartesian.py` | dam break of T1 with `W 2 0.05`: the effect of the diffusion on eta on the curvilinear and on the Cartesian grid (same nodes) differs by 6 % for a small dam (`HDAM=1.05`), by 54 % for the 2 m dam (the conservative stress with the gradient of h) |
| `t1_dambreak_rotation.py` | dam break in a closed basin: Cartesian kernels, curvilinear kernels on the same nodes, and on the grid rotated by 90 deg and shifted agree per cell index to round-off (1.8e-15) |
| `t1b_channel_rotation.py` | river reach with discharge inflow and stage outflow: curvilinear = Cartesian on the same nodes (3.6e-15); rotated by 30 deg at UTM coordinates with an emerging island: equal to truncation level with WENO (the nonlinear weights act on the Cartesian components), to 8e-11 with a linear reconstruction (`A 211 0`) |
| `t2_lake_at_rest.py`, `t2b_variants.py` | lake at rest on a wavy, skewed grid: velocities 2.5e-15 fully wet; with an emerged island 1e-5 at the shoreline, the same as on the Cartesian grid (wet/dry treatment of SFLOW) |
| `t4_nh_bous_rotation.py` | `A 220 1`, `2`, `3`, `4`: dam break, curvilinear = Cartesian on the same nodes and rotated by 90 deg to 1e-15 (`A 220 4`: curvilinear = Cartesian to the solver tolerance, rotated = curvilinear to 1e-15); river reach with an island rotated by 30 deg (`A 211 0`) to 1e-10 |
| `t5_amr_rotation.py` | mesh refinement, dam break with patches following the surface (`G 20`, regrid every 4 steps): hydrostatic, NH, Boussinesq (small step), subcycling, two ranks with placement: curvilinear = Cartesian = rotated by 90 deg per level-0 cell to round-off, identical patches, volume to round-off |
| `t6_lake_amr.py` | lake at rest on the wavy, skewed grid with static and regridded patches (shoreline flagging), hydrostatic, NH, Boussinesq, subcycling: velocities at round-off (2e-15), volume constant |
| `t8_shoreline_amr.py` | lake at rest with an EMERGED island (1.6 m, water level 1.0 m) and mesh refinement: patch boundaries through the island and next to the shoreline, the whole island refined, two levels, subcycling, NH, Boussinesq, regridding along the shoreline, two ranks: see the table below |
| `t7_skew_stability.py` | Boussinesq on strongly skewed grids (up to 45 deg non-orthogonality), one level, and with patches: a small disturbance stays bounded |
| `t3_bend.py` | 90 deg bend, Manning friction: discharge through the outlet, superelevation 0.0099 m against U^2 B/(g R) = 0.0102 m; 2 ranks equal 1 rank to round-off |
| `make_regression_cases.py` | writes the `Regression/cases/sflow_2d_curv_*` cases |
| `read_curv.py` | reads the CURV section of a geometry file |

`t8_shoreline_amr.py` (400 steps, 9.6 s; the largest velocity of level 0 in m/s and the change of the water volume;
one level: 5.9e-7, and 1.3e-4 for one level at the dx of level 2, the film of the wetting and drying):

| case | before | now | volume change before / now |
|---|---|---|---|
| `cut`: patch boundary through the island | 4.1e-2 | 7.7e-7 | 2e-14 / 3e-14 |
| `cut_y`: the same along y | 9.2e-2 | 8.8e-7 | 2e-14 / 2e-14 |
| `edge`: patch boundary 0.6 m off the shore | 4.4e-2 | 8.9e-7 | 2e-14 / 2e-14 |
| `all`: the whole island refined | 4.7e-2 | 8.9e-7 | 3e-14 / 2e-14 |
| `sub`: `cut` with `G 7 1` | 5.9e-2 | 1.7e-6 | 3e-14 / 3e-14 |
| `lev2`: `cut` with two levels | 3.6e-2 | 2.5e-4 | 8.2e-7 / 2e-14 |
| `nh`, `bous`: `A 220 1`, `A 220 4` | 2.9e-2, 4.5e-2 | 7.8e-7, 7.4e-7 | 3e-14 / 2e-14 |
| `shore`: regridded along the shoreline (`G 22 1`, `G 2 5`) | 2.6e-2 | 8.4e-7 | 2e-14 / 2e-14 |
| `shore2_sub`: two levels, `G 7 1` | 1.9e-2 | 5.7e-4 | 1.9e-6 / 3e-14 |
| `mpi2`, `shore_mpi2`: two ranks | 4.1e-2, 2.6e-2 | 7.7e-7, 8.4e-7 | 2e-14 / 3e-14 |
