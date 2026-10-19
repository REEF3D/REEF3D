# Sediment transport on curvilinear grids, mesh refinement and multi-block grids (P5a, P5b): tests

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
| `t7_talmon_bend.py` | P6c: run 3 of Talmon & de Graaff (1989) in the curved flume of the LFM, TU Delft (`data/talmon1989_run3.json`: 180 deg, R 4.10 m, W 0.5 m, Q 5 l/s, a0 0.051 m, C 19.2, D50 90 um, 75 % suspended): SFLOW, Engelund-Hansen total load as bed load, `S 86`, `S 95` with `E_spir 1`, the measured a/a0 of 48 sections x 9 positions against the depth of the model averaged over the last 300 s (the bed forms and the migrating bars are averaged out in the measurements); `KS`, `DS`, `MF`, `ESPIR`, `TAVG`, `PNG` | `DS=0.05` (10 cells across, 900 s flow, MF 100): the second half of the bend (sections 25-45) within 0.09 a0 rms, outer - inner 0.59 a0 (measured 0.56); the overshoot after the entrance too weak and too early (0.78 a0 at 23 deg against 1.45 a0 at 50 deg: the spiral flow follows the curvature at once, no adaptation length).  `DS=0.025` (20 cells): a second transverse mode (a trough 0.1 m from the inner bank), outer - inner only 0.27 a0: the result depends on the grid; missing physics are the lag of the suspended load (adaptation length h u/ws about 1.3 m, 25 a0) and of the spiral flow |
| `make_regression_cases.py` | writes the regression cases `*_sed_*` (P5b: `*_sed_*bend*`) | |

P5b with `XCTRL "S 86 1;S 95 1"` (NHFLOW `S 86 1`): t1 Cartesian = curvilinear = turned by 90 deg to round-off (SFLOW
and NHFLOW, `S 31 1`/`2`, `S 32 1`/`2`, 1 and 4 ranks); t2 bed volume conserved (1e-19) on all grids, also with the
steep mound and the sand slide; t3 (`SOLVER`, `XCTRL` added) 2 and 4 ranks bit for bit with one rank (NHFLOW `A 509 0`,
SFLOW `A 209 0`); t4 leaf volume conserved (1e-17) for 1-2 levels, adaptive, 2 ranks, curvilinear = turned; t5 four
blocks = one block to round-off on 1-3 ranks (NHFLOW with refinement on 2 ranks: as the flow alone).

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
