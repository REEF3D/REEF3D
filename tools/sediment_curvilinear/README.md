# Sediment transport on curvilinear grids, mesh refinement and multi-block grids (P5a, P5b, P7, P8a): tests

Phase P5a of the multi-block curvilinear concept: the sediment module of REEF3D (`S 10 1`) with
REEF3D::NHFLOW and REEF3D::SFLOW on curvilinear grids (`G 50 2`, one block from `curvgrid.dat`), on the
patches of the mesh refinement (`G 1`, REEFAMR) and on multi-block grids (`G 50 1`, DIVEMesh `R 20 1`).

## Running a case

The sediment keys are those of the Cartesian runs.  On a curvilinear grid the `F 72` boxes, the `S 71`/`S 72`
and `S 77` windows and the `S 74` zones are in physical coordinates (cell centres).

Available on one curvilinear block: the bed-load formulas (`S 11`), the conservative Exner forms `S 31 1`
and `S 31 2` with `S 32 1` (FOU) and `S 32 2` (CDS), the bed filter (`S 100`, `S 101`), the bed-slope and
critical-shear reductions (`S 80`-`S 84`), the bed-load direction (`S 82`), the non-equilibrium bed load
(`S 33`), all sand slides (`S 90 1`-`5`), the time-step modes (`S 15`), the multi-fraction bed (`S 51`),
the NHFLOW suspended load (`S 12`) and the coupled RK2 step (`S 10 12`), several ranks.  On the mesh refinement and on multi-block grids: the
same without the suspended load and the multi-fraction bed.  Refused with a message: the advective Exner
form `S 31 3`, the discretisations `S 32 3`-`5`, the VRANS bed (`S 10 2`); on the mesh refinement and on
multi-block grids also `S 12` (NHFLOW) and `S 51`.

NHFLOW on curvilinear grids takes the time step from the speed `sqrt(U^2+V^2)` (independent of the
direction of the grid; before P5a from the larger component, which depended on the turn of the grid).

## How it works

- Geometry (`sediment_curv.h`): the face normals and lengths, cell areas and centres of the geometry
  layer.  The bed velocity of a cell `ub`, `vb` (Cartesian components) gives the face velocities `P`, `Q`
  along the face normals (`face_normals`, boundary conditions of `P`, `Q` as on the Cartesian grid); the
  cell velocity for the bed-load direction is contravariant from the four faces (`cellvel`); the upwind
  sign of an Exner face is the normal velocity over the speed with the tangential part from the four
  crossing faces (`face_sign`); the bed gradient is the chain rule of the index derivatives with the
  one-sided stencils of the Cartesian code next to closed cells (`bedgrad`).
- The Cartesian functions are kept as they are: each starts with `if(p->geo_curv==1) { X_curv(...);
  return; }` and the `*_curv` copy follows it with the geometry of `sediment_curv.h` (bed slope, bed-load
  direction, bed shear of SFLOW, `bedload_EH`, Exner `topovel`, non-equilibrium, time step, filter, the
  sand slides, the multi-fraction bed).  Cartesian runs do not change.
- Exner (`sediment_exner_curv.cpp`): the divergence of the face fluxes `q*L` over the cell area; the time
  step from the cell widths (area over the mean length of the opposite faces).
- The filter on a curvilinear grid is the symmetric product `(Sy Sx + Sx Sy)/2` of the one-dimensional
  passes of the Cartesian filter (`sediment_xfilter`): conservative, the same transfer on a uniform
  grid, and independent of which index family is `i`.
- The sand slides use the distance of the cell centres and the area ratio of donor and receiver
  (`SLIDE_D`, `SLIDE_ARC`): the slid volume is conserved on any grid.
- Mesh refinement and multi-block grids (`nhflow_amr_sed.cpp`, `sflow_amr_sed.cpp`, `sediment_amr.cpp`,
  `sediment_exner_amr.cpp`): every grid (level 0, every block and patch) has its own `sediment_f`; a
  sediment step runs all grids: the physics of every grid, qb of the covered cells from their children
  (mean per area), qb around the patches (same level, else linear from the coarser level), the Exner face
  fluxes of every grid, the coarse faces next to a finer patch from the patch (flux matching: the bed
  volume over the leaf cells is conserved), one `dtsed` for all grids, the bed, the restriction of the bed
  to the covered cells, the bed to the flow and to the cells around the patches.
- The filter and the sand slide over all grids (`reefamr_sed.h`): the grids run each step of them in
  turn and exchange their ghost cells between the steps (`reefamr::xslice`; the material slid into the
  cells around a grid is added to the cell it stands for, `reefamr::xslice_add`, remote cells through
  `xrun_rev`).  The interfaces of the grids of a level (blocks, the ring of the native grid, patches that
  meet) are open as the inside of a grid; covered cells and cells around a patch with a coarser source are
  closed.  The intermediate fields of the symmetric filter belong to the index frame of their grid and
  are swapped in the cells around a grid whose source block is turned against it.
- The flow keeps the water level WL of a column through a bed update (eta is that of the old bed until
  the first stage, as on the uniform grid).  The NHFLOW patches therefore take the water level of the
  cells around them from the source with sediment (`fill_stage`), not from eta plus their depth.

## Bend morphology (P5b)

Two corrections of the bed-load direction for river bends, both off by default, on Cartesian and curvilinear
grids, the mesh refinement and multi-block grids (`sediment_exner_bend.cpp`):

| key | | default |
|---|---|---|
| `S 86 1` | transverse bed slope effect: `q = qb (s - grad(z)/f(theta))`, `f = A_sh theta^B_sh (d/h)^C_sh` (Talmon, van Mierlo and Struiksma 1995) | off |
| `S 87 S 88 S 89` | `A_sh`, `B_sh`, `C_sh` | 9, 0.5, 0.3 |
| `S 95 1` | spiral-flow correction (SFLOW): the transport direction turned towards the centre of curvature of the streamlines by `tan(delta) = alpha_I h/R_s`, `alpha_I = E_spir 2/kappa^2 (1 - sqrt(g)/(kappa C))`, `C = sqrt(g) U/u*` from the bed shear | off |
| `S 96` | `E_spir` | 1.0 |

- The slope part `-qb/f dz/dn` is a face flux of its own (central: the mean of `qb/f` of the two cells, `dz/dn`
  between the cell centres): conservative and diffusive.  It is closed at walls, at the domain edges, to cells
  without an erodible bed and to cells whose bed is held (inflow guard `S 78`).  It is an explicit diffusion of the
  bed: `dtsed` is kept below `0.4 (1-n)/(S35 max(qb/f (1/dx^2+1/dy^2)))` (all grids of the mesh refinement).
- The spiral angle of a face is the mean of its two cells, 0 at closed faces (the turn would open them) and where the
  two cells cancel to round-off (a symmetry line: the upwind cell of `S 31 2`/`S 32 1` must not follow the sign of
  round-off).  `S 31 2` turns the cell vectors as well.  The curvature `1/R_s = u x (u.grad)u/|u|^3` is signed
  (positive turning left); `tan(delta)` is clamped to 1.
- Mesh refinement and multi-block grids: the cells around a grid take `qb/f`, `tan(delta)`, the bed and the open
  flag of the slope flux from their source like qb (`sediment_exner::amr_field`), so that the two grids of a face
  compute the same flux (the flow bed of the cells around a fresh patch can lag behind the sediment bed after a
  regrid; their own inflow guard is not that of their source).
- NHFLOW resolves the secondary flow: `S 95` is refused there; `S 86` and `S 95` are refused with `S 31 3`.
- Equilibrium of a long bend: the transverse flux vanishes, `dz/dn = f(theta) alpha_I h/r` (t6).  With
  `E_spir 1` and the LFM flume the inner bank falls dry (`tan(delta)` about 0.4, `f` about 0.9); Delft studies use
  `E_spir` of about 0.4-0.5.

## Lags of the bend morphology (P7)

Two steady relaxations along the depth-averaged flow, both off by default (SFLOW; Cartesian and curvilinear grids,
the mesh refinement and multi-block grids; `sediment_exner_lag.cpp`):

    f + L s.grad(f) = f_e ,    s = (u, v)/|U|

| key | | default |
|---|---|---|
| `S 95 2` | spiral flow with adaptation (Kalkwijk & Booij 1986; the form of Delft3D and HEC-RAS 2D): the intensity `I` relaxes towards `I_e = h |U| K_c` (`K_c` the signed curvature of the streamlines) over `L_a = S97 (m-1) h/(2 kappa^2)`, `m = kappa C/sqrt(g)`, `C = sqrt(g) U/u*` from the bed shear (at most 100); the bed load turns by `tan(delta) = alpha_I I/|U|` (`alpha_I` of `S 95 1`).  `S 95 2` with `S 97 0` is `S 95 1` | off |
| `S 97` | factor of `L_a` | 1.0 |
| `S 64 f_s` | suspended part of the transport (0..1): the formula (`S 11`) gives the capacity `q_e`; `(1-f_s) q_e` is bed load (equilibrium, or `S 33`), `f_s q_e` suspended load that relaxes over `L_s = h |U|/(S65 w_s)` (Galappatti & Vreugdenhil 1985, depth-averaged; `w_s` the fall velocity, `S 23`/`S 25`) | 0 |
| `S 65` | `c_a/c_mean`, the near-bed over the depth-averaged concentration of the suspension: shortens `L_s`.  For a Rouse profile with `Z = w_s/(kappa u*)` and `a = 0.05 h`: 2.3 (Z 0.3), 4.2 (Z 0.57, t7), 8.8 (Z 1) | 1.0 |

- First-order upwind along `s` (curvilinear: the contravariant components of `s`), solved in every sediment step to
  convergence by Gauss-Seidel with alternating sweep directions (as the non-equilibrium bed load `S 33`: at most
  `S 49` sweeps, the change below 1e-6 of the largest `f_e`), started from the solution of the step before.  The upwind
  neighbour is a sediment cell (of the grid, the next rank or the cells around the grid): at the inflow, at walls and
  open edges the cell takes its own `f_e`.  Where `|U| < 1e-4 m/s` `I_e` and `L_a` are 0.
- The spiral flow (`S 95`) and the transverse slope (`S 86`) turn the total load, bed load plus the lagged suspended
  load (`qb/f(theta)` of the slope flux from the total; it is computed again after the passes below).  Turning the
  bed load only (tried in P7) moved the overshoot of t7 far
  downstream (section 30 against 16, `S 65 1`) or left the bed upstream of the bend and after the overshoot too
  steep (`S 65 4`).  The Exner fluxes carry bed load plus suspended load: the bed volume stays conserved.  Output:
  `ST_spiral` (`S 95 2`) and `ST_qbs` (the suspended load) with `P 76 1`.
- Mesh refinement and multi-block grids: every grid relaxes its own cells and holds the cells around it (filled from
  their sources: `Isp`, `qsn` are fields 5 and 6 of `amr_field`) and its covered cells (the mean of their children).
  After the fills the relaxations run again from the filled cells (`lag_refresh`), passes until the cells around the
  grids no longer move the solution (the change of the first sweep below the tolerance, at most `S 49` passes; the
  log prints the passes).  Along a river the passes follow the grids downstream (2-4 per step in t9); in recirculating
  flow with `L` longer than the domain (the dam break of t5 with `L_s` 2 m) up to 150.
- P8a: the non-equilibrium bed load (`S 33`) is the same relaxation and runs the same way on the mesh refinement and
  the multi-block grids: the cells around every grid and its covered cells held (`lagfix`, `qbn` is field 7 of
  `amr_field`), the passes of `lag_refresh` from the bed load on (non-equilibrium, suspended part, `qb/f`).  Before,
  every grid relaxed alone (the cells around it overwritten afterwards): four blocks differed from one by 1e-3 of the
  bed change.  Without the mesh refinement and blocks `S 33` is unchanged (bit for bit).
- P8b: `S 76 xs xe ys ye` holds the bed of the cells in a box (any number of boxes, physical cell centres) while
  their load is computed and carried downstream: the sediment supply of an inflow reach on grids where the x
  window `S 71` does not fit (the Nidelva, `../harbour_multiblock/README.md`).
- Refused with a message: `S 64` with NHFLOW (it resolves the suspended load, `S 12`), with the multi-fraction bed
  (`S 51`) and with `S 31 3`.

## Tests

Binaries: `REEF3D_BIN`, `DIVEMESH_BIN` (DIVEMesh with `R 20`).  Runs go to `./runs`.

| script | what | result (P5a) |
|---|---|---|
| `t1_sed_rotation.py [steps]` | bed-load dam break over a flat bed: Cartesian, curvilinear on the same nodes, turned by 90 deg (`SOLVER`, `A520`, `XCTRL`, `NP`) | NHFLOW and SFLOW: bed per cell index to 4e-15 of the change, curvilinear = Cartesian; the same with two fractions (`S 51`) and with `S 10 12` |
| `t2_sed_river.py [steps]` | river reach over a mound: discharge inflow, stage outflow, friction; grid on the Cartesian nodes, turned by 90 deg, turned by 30 deg at UTM, wavy; bed volume against the fluxes through the open faces (`MOUND`, `SLOPE`, `XCTRL "S 90 1"`) | rot90 to 2e-14 of the change; volume change 1e-18 (conserved) on all grids, also with the sand slide; rot30 at UTM differs by the round-off of the flow at these coordinates (4e-4 of the change, adaptive dt) |
| `t3_sed_mpi.py [steps]` | (a) curvilinear against Cartesian per rank on 2 and 4 ranks; (b) the river turned by 30 deg with the sand slide on 1, 2, 4 ranks | (a) 1e-18; (b) bit for bit with `A 509 0` (the Gauss-Seidel bed smoothing `A 509` of NHFLOW depends on the ranks) |
| `t4_sed_amr.py [steps]` | the dam break of t1 with mesh refinement: static box (`LEV` levels) or patches following the surface (`ADAPT=1`), Cartesian, curvilinear, turned by 90 deg; bed volume over the leaf cells (`SOLVER`, `NP`, `XCTRL`) | leaf volume to 2e-17 of the change for NHFLOW and SFLOW, 1-2 levels, adaptive, 2 ranks, with the sand slide; curvilinear and turned = Cartesian to 1e-17 (SFLOW static `LEV=1`: 1e-4 of the change in the coarse cells next to the coarse-fine face, see the notes) |
| `t5_sed_mb.py [steps]` | the four turned quadrants of P4 (multi-block grid) against one block: dam break with bed load over a bump (`BUMP`), the filter and the sand slides (`XCTRL "S 90 n"`), a static box over the vertex (`LEV`), several ranks (`NP`), the blocks not turned (`NOTURN`), one block turned (`ROTONE`) | NHFLOW and SFLOW, `LEV` 0 and 1, `S 90 1`-`5`, 1-3 ranks: level 0 and every fine cell to 2e-16 (round-off), the bed volume conserved; SFLOW `LEV=0` against the plain one-block run (no mesh refinement).  NHFLOW with refinement on 2 ranks: as the flow alone (P4, see the notes) |
| `t6_sed_bend.py [steps]` | P5b: the LFM bend of Delft Hydraulics (180 deg, R 4.25 m, B 1.7 m, h 0.2 m, Q 0.17 m3/s, d50 0.78 mm, Chezy 26.4), SFLOW, curvilinear grid, `S 86 1`, `S 95 1`, `E_spir 0.4`; the transverse slope of the cross-sections against the equilibrium of the model (`DS`, `MF`, `ESPIR`, `NOSPIRAL`, `NOSLOPE`) | `DS=0.2 MF=10` (2760 s flow, morphological factor 10; P6c with the bed friction `A 218 1`, Manning `ks` of the Chezy at `h0`, the P5b runs had none): the transverse slope of the second half of the bend (90-180 deg) is 1.02-1.28 of the equilibrium (mean 1.16; P5b without friction 0.90-1.15, mean 1.05), the outer bank 0.12 m deeper, the point bar 0.09-0.12 m higher; without `S 95` no transverse slope (0.02 of it), without `S 86` the outer bank scours without bound (0.65 m, 7.5 times the slope; P5b)|
| `t7_talmon_bend.py` | P6c: run 3 of Talmon & de Graaff (1989) in the curved flume of the LFM, TU Delft (`data/talmon1989_run3.json`: 180 deg, R 4.10 m, W 0.5 m, Q 5 l/s, a0 0.051 m, C 19.2, D50 90 um, 75 % suspended): SFLOW, Engelund-Hansen total load as bed load, `S 86`, `S 95` with `E_spir 1`, the measured a/a0 of 48 sections x 9 positions against the depth of the model averaged over the last 300 s (the bed forms and the migrating bars are averaged out in the measurements); `KS`, `DS`, `MF`, `ESPIR`, `TAVG`, `PNG` | `DS=0.05` (10 cells across, 900 s flow, MF 100): the second half of the bend (sections 25-45) within 0.09 a0 rms, outer - inner 0.59 a0 (measured 0.56); the overshoot after the entrance too weak and too early (0.78 a0 at 23 deg against 1.45 a0 at 50 deg: the spiral flow follows the curvature at once, no adaptation length).  `DS=0.025` (20 cells): a second transverse mode (a trough 0.1 m from the inner bank), outer - inner only 0.27 a0: the result depends on the grid; missing physics are the lag of the suspended load (adaptation length h u/ws about 1.3 m, 25 a0) and of the spiral flow; P7 adds them (below) |
| `t7_compare.py [--grids] [--runs] [--labels] [--png]` | P7: the t7 runs (`runs/t7_talmon<grid><run>/t7.npz`) against the measurements and each other: rms, outer - inner, the deepest point, the difference of the grids; a figure of the outer bank, outer - inner and the inner bank along the sections | see below |
| `t8_bend_lags.py equiv\|profile [steps]` | P7: the LFM bend of t6 (`DS 0.1`): `equiv` `S 95 2` with `S 97 0` against `S 95 1` after 3000 steps with the bed moving; `profile` the bed held, the spiral intensity `ST_spiral` on the centre line with `S 97 1` against the relaxation of `I_e` (the run with `S 97 0`) by the scheme of the model and by the exact solution (`PNG`) | `equiv`: bed 1.8e-11 (change 0.07 m), the face velocities 1.6e-10 (round-off of `alpha h U K_c/U` against `alpha h K_c`); `profile`: `L_a` 1.48 m (C 26.4, h 0.2 m), model = upwind scheme to 0.4 % of `I_e` max, = exact relaxation to 0.4 %; 63 % of the build-up 2.04 m after the bend entrance (model, scheme, exact; `I_e` itself rises over 1 m), of the decay 2.05 m after the exit (exact 1.95) |
| `t9_bend_blocks.py [steps]` | P7: the LFM bend as one block and as two blocks along the flow (block 1 the inflow reach with `G 54 1 1 Q`, block 0 the bend and the outflow reach), bed and water level per cell after 1500 steps (500 sediment steps, MF 100): `S 95 1`, `S 95 2`, `S 95 2` + `S 64 0.5`; `CASES`, `NA` | the flow (bed held) 5e-16 m in the water level; `S 95 1` 8e-9 m in the bed (bed change 0.054 m), `S 95 2` with `S 97 0` and `0.1` the same (8e-9), with `S 97 1` 4.6e-6, `S 95 2` + `S 64 0.5` 1.6e-6, 2-4 passes per step; the curvature `I_e` of the first cell at the inflow differs between the native inflow and the open face `G 54` of a block (their ghost velocities), the lag carries it downstream: with `I_e` set to 0 in the first metre (a test build) 9.5e-9 and 2.2e-8.  The bed volume change of both the same to 1e-9.  (Sediment leaves the domain only through the native outflow: an open face `G 54` of a block is closed for the sediment, three blocks with the outflow on block 2 fill its end with sand) |
| `t9_bend_blocks.py` P8a (`LAYOUT=three`: block 1 the inflow reach `G 54 1 1 Q`, block 0 the bend, block 2 the outflow reach `G 54 2 2 h`; `NOGUARD`, `MF`, `PERTURB`, the cases `eh*` with Engelund-Hansen, `*noneq` with `S 33`) | the sediment through the open face of block 2 | before P8a three blocks differed from one by 0.2 m (the end of block 2 silted up, the face closed for the sediment); now 1.7e-4 m with MF 100 and 1.6e-6 m with MF 10 (bed change 0.05 and 0.013 m), the flux through the outflow face and the bed of the outflow reach equal to one block to 1e-12 over the first 150 sediment steps; the remaining difference starts at the interface of blocks 0 and 2 next to the outer wall, where the bed changes fastest after the bend, and grows with the morphology (two blocks with the interface in the straight inflow reach: 8e-9).  `S 33 1` across the interfaces (`ehnoneq`, MF 10): two blocks 2.0e-9, three blocks 2.9e-7.  Without `S 71` the one-block run puts 1.2e-2 m of sand into the cell at its inflow corner next to the wall (j = 0 only), the block not: the native inflow corner |
| `t5_sed_mb.py` P8a (`XCTRL "S 33 1"`, `"S 11 4;S 33 4;S 40 2"`) | `S 33` across the blocks, with a refinement box, on 2 ranks | `S 33 1` (Ls 0.1 m): four blocks = one block to 7.7e-10 (P7: 1.9e-6), with refinement 6.4e-10, 2 ranks = 1 rank; `S 33 4` with Ls 2h (0.4 m, recirculating dam break, Engelund-Hansen): 5.5e-7 (P7: 2.2e-6), with refinement 1.4e-6; with the relaxation tolerance 1e-12 (a test build) 8e-12, with refinement 1.3e-11, 2 ranks 8e-12: the coupling is exact, the rest is the tolerance 1e-6 of the sweeps.  Bed volume conserved |
| `make_regression_cases.py` | writes the regression cases `*_sed_*` (P5b: `*_sed_*bend*`, P7: `*_sed_*lags*`, P8a: `sflow_2d_mb_sed_open_faces_noneq`, `sflow_2d_mb_sed_noneq_amr_mpi2`) | |

P5b with `XCTRL "S 86 1;S 95 1"` (NHFLOW `S 86 1`): t1 Cartesian = curvilinear = turned by 90 deg to round-off (SFLOW
and NHFLOW, `S 31 1`/`2`, `S 32 1`/`2`, 1 and 4 ranks); t2 bed volume conserved (1e-19) on all grids, also with the
steep mound and the sand slide; t3 (`SOLVER`, `XCTRL` added) 2 and 4 ranks bit for bit with one rank (NHFLOW `A 509 0`,
SFLOW `A 209 0`); t4 leaf volume conserved (1e-17) for 1-2 levels, adaptive, 2 ranks, curvilinear = turned; t5 four
blocks = one block to round-off on 1-3 ranks (NHFLOW with refinement on 2 ranks: as the flow alone).

P7 against Talmon & de Graaff (1989) run 3 (`t7_talmon_bend.py` with `XCTRL "S 95 2;S 64 0.75;S 65 4"`, 75 % of the
load suspended, `S 65 4` from the Rouse profile of the run (`w_s` 7.3 mm/s, `u*` 0.032 m/s, Z 0.57; `L_s` up to 0.74 m,
`L_a` up to 0.55 m); `t7_compare.py` for the table and the figures).  The two grids now agree: 10 and 20 cells across
differ by 0.034 a0 rms (P6c `S 95 1`: 0.12, `S 95 2` alone: 0.15, both with a second transverse mode on 20 cells).

| 10 / 20 cells across | P6c `S 95 1` | `S 95 2` | P7 `S 95 2` + `S 64 0.75` + `S 65 4` | measured |
|---|---|---|---|---|
| rms a/a0, all 432 points | 0.149 / 0.208 | 0.150 / 0.247 | 0.090 / 0.097 | |
| rms, sections 6-24 (overshoot) | 0.195 / 0.266 | 0.194 / 0.296 | 0.076 / 0.082 | |
| rms, sections 25-45 | 0.089 / 0.160 | 0.117 / 0.212 | 0.111 / 0.118 | |
| outer - inner, sections 12-18 | 0.53 / 0.28 | 0.78 / 0.24 | 1.25 / 1.16 | 1.25 |
| outer - inner, sections 25-45 | 0.59 / 0.27 | 0.47 / 0.24 | 0.81 / 0.80 | 0.56 |
| deepest outer point | 9: 1.48 / 10: 1.27 | 13: 1.55 / 12: 1.47 | 15: 1.65 / 15: 1.59 | 16: 1.73 |
| inner-bank minimum | | | 14: 0.23 / 14: 0.35 | 15: 0.24 |

The overshoot (position, amplitude, the inner bank) is now right on both grids.  Still off: after it the bed stays too
steep (outer - inner 0.8 against 0.56 in the second half), the measured trough at sections 25-33 is missing; `S 95 2`
alone has the trough on 10 cells but not on 20.  The `DS=0.025` run used a build that differed from the final one only in
the upwind rule at walls (ghost cells beyond a wall were upwind neighbours); on 10 cells the two give the same table.

P7 with `XCTRL "S 95 2;S 64 0.5;S 86 1"` (SFLOW): t1 turned by 90 deg = curvilinear to round-off (5e-16 of the change),
curvilinear = Cartesian to 1.2e-8 of the change (the tolerance of the relaxations: the two grids stop the sweeps at
other round-off).  t5 with `S 95 2`, `S 97 0.1`, `S 64 0.5`, 40 steps: four blocks = one block to 3e-9 m (change 1e-3
m), with a refinement box 5e-8 (both levels), 2 ranks = 1 rank; the bed volume conserved (1e-15).  The dam break
amplifies the differences of the relaxations (tolerance 1e-6): `S 95 2` with `S 97 1` (`L_a` up to 6 m in a 2 m box,
recirculating flow) 1e-5 after 60 steps, 1.3e-8 after 20.  Defaults: the 107 cases `sflow_*` and `*sed*` bit for bit
with P6c.

## Notes

- `A 509` (NHFLOW bed smoothing at the start, Gauss-Seidel): the result depends on the rank layout
  unless `G 52 1` (smoothing from the old values, set on multi-block grids) or `A 509 0`.
- NHFLOW with mesh refinement depends on the rank count already without sediment (P3/P4: four blocks with a
  refinement box over the vertex on 2 ranks, about 1e-3 in eta after 20 steps; the same with the P4
  binary); the sediment step itself is the same on any number of ranks (SFLOW with refinement on 1-3
  ranks bit for bit, NHFLOW multi-block without refinement on 1-3 ranks to round-off).
- SFLOW, Cartesian against curvilinear with a static refinement box (`t4 LEV=1`): the Cartesian bed shear
  of SFLOW takes the cell velocity from the face velocities `P`, `Q` of the flow, the curvilinear one the
  cell velocities; at a coarse-fine face the flow's face velocity is the matched fine one, so the bed load
  of the coarse cells next to it differs (1e-4 of the bed change in t4).
- `S 31 2` and the sand slide `S 90 3` (steepest descent: ties between equally steep neighbours) are
  sensitive to round-off on the Cartesian grid as well; `S 82 > 0` uses the flow angle (0 at rest).
- SFLOW smooths the bed of the native grid (`A 209`) but not that of the blocks 1..: use `A 209 0` to
  compare a multi-block grid with one block.
- Found with P5b (not changed, they do not affect the results without the bend terms): the cells around a fresh
  NHFLOW patch hold a flow bed that can lag behind the sediment bed of their source for one step after a regrid
  (the bed slope of qb next to the patch edge sees it); the inflow guard (`S 78`) of the cells around a block on
  another rank is 0; SFLOW `A 209` (bed smoothing) depends on the rank layout like NHFLOW `A 509`.
- `t6` with a morphological factor of 100 runs at a bed celerity of 50-100 times the flow speed (`c/U` in the log):
  far from quasi-steady; use 10 or less.
- On the uniform NHFLOW grid eta keeps the old beduntil the first stage after a bed update (the stage
  uses it in the kinematic condition and the convection); the mesh refinement keeps that and takes the
  water level of the cells around the patches from their sources.
