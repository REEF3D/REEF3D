# REEF3D regression test suite

Small, fast cases that exercise the main code paths, plus tools to

- **A/B-compare two binaries bitwise** (for refactors and clean-up: "did anything change?"), and
- **check against stored references** with a tolerance (long-term regression testing).

Pure Python 3 (standard library only) + bash. CFD, NHFLOW, FNPF and SFLOW cases; the layout is
solver-independent (`"solver"` in `case.json`).

## How it works

`src/regression_dump.{h,cpp}` is a small output class called from the CFD, NHFLOW, FNPF and SFLOW driver loops. It is
**inactive unless `REEF3D_REGRESSION_DIR` is set**, so normal runs are unchanged. When active, each
rank writes into that directory

| file | content |
|---|---|
| `steps_r<rank>.txt` | one line per step: count, simtime, dt and rank-local Σu², Σv², Σw², Σp², Σφ², Σνt² in C99 hexfloat (exact) |
| `state_<count>_r<rank>.bin` | full double-precision state (u, v, w, p, φ, νt, k, ε/ω, ρ, ν, topo, conc, flag4) at the initial and final step (and every n steps with `REEF3D_REGRESSION_EVERY=n`) |

`regression.py` writes the case inputs, runs DIVEMesh and REEF3D, and compares:

1. the exact final state, field by field (max abs diff, number of differing cells),
2. the per-step record, giving the **first time step where the runs diverge**,
3. the normal text output (wave gauges, probes, forces, 6DOF, CPM sediment log, NHFLOW particle log, NHFLOW boom forces) with a tolerance.

Result per case: **identical** (bitwise) · **close** (within `--rtol/--atol`) · **different** · **failed**.

## Quick start

```bash
cd Regression

# build two binaries with identical, reproducible flags (-O2 -ffp-contract=off, no -march=native/LTO)
./build_reef3d.sh -r hans_dev -o ~/reef3d_reg/ref      # a git revision (exported, working tree untouched)
./build_reef3d.sh             -o ~/reef3d_reg/new      # the working tree (incremental rebuilds)
# no hypre by default (as the Makefile); -H <hypre prefix> builds with hypre

export REEF3D_MPIRUN="mpirun"          # e.g. "mpirun --oversubscribe" if fewer cores than ranks
# or cap the ranks of all cases for an A/B run on a small machine: --max-np 2 (env REEF3D_MAX_NP);
# both binaries then use the same reduced decomposition, so ab stays valid (not for bless/check)

# run both and compare (default: must be bitwise identical)
./regression.py ab --ref-bin ~/reef3d_reg/ref/REEF3D --new-bin ~/reef3d_reg/new/REEF3D \
                   --divemesh ~/Codelite/DIVEMesh/bin/DiveMESH --out ~/reef3d_reg/ab --tags quick
```

The report is printed and written to `<out>/new/compare.md` (+ `compare.json`).

A case can depend on earlier runs (hydrodynamic coupling): `"chain": [{"case": "<case>", "copy":
["REEF3D_FNPF_STATE"]}]` runs that case first in `<rundir>/_chain` and copies the listed folders
into the run directory before DIVEMesh; `"keep_keys": ["P 40", "P 41"]` keeps print keys the suite
normally removes (state files for such a stage).

Other commands:

```bash
./regression.py list [--tags quick]                     # cases, tags, what they cover
./regression.py run  --reef3d BIN --divemesh DM --out DIR [--cases 'cfd_2d_*'] [--tags ...] [--steps N] [--keep]
./regression.py compare REF_DIR NEW_DIR [--require identical|close|different] [--rtol 1e-6] [--atol 1e-12]
./regression.py write --out DIR --cases cfd_2d_nwt      # only write the resolved input files
./regression.py bless DIR [--note "why"]                # store compact references in cases/<case>/reference.json
./regression.py check DIR [--rtol 1e-6]                 # check a run against the stored references
```

## Workflow for code changes

- **Refactor / removal (must not change results):** `ab --require identical`. Any difference
  shows the first diverging step and the fields affected.
- **Bug fix / physics change:** `ab --require close` or look at the report; then re-bless the
  affected cases with a note: `./regression.py bless DIR --cases <...> --note "F1 diffusion weight"`.
  The note and date are stored in `reference.json`, so the reference history is in git.
- **Bitwise comparisons need identical build conditions**: same compiler, MPI, flags
  (`build_reef3d.sh` fixes them) and the same number of ranks. Stored references
  (`check`) are compared with a tolerance and are meant to work across machines.

## Cases

Each case is a directory in `cases/` with a `case.json` and either its own `control.txt`
(DIVEMesh) + `ctrl.txt` (REEF3D), or `"base": "<other case>"` plus key overrides:

```json
{
 "base": "cfd_2d_nwt",
 "ctrl_set": {"N 40": "13"},          // replace all lines with this key (null removes, list = several lines)
 "ctrl_add": ["P 61 5.0 0.025 0.4"],   // append lines
 "control_set": {"B 1": "0.025"},      // same for the DIVEMesh control.txt
 "np": 4, "steps": 40, "dump_every": 0,
 "files": ["floating.stl"],            // extra input files copied into the run directory
 "tags": ["quick", "2d"],
 "description": "...", "covers": ["code paths this case is meant to exercise"]
}
```

The runner always sets `M 10` (ranks) in both files, `N 45` (steps), `N 41 1e9`, and removes the
VTU/state print keys (`P 20/30/40/41/42`), so runs are short and output stays small.

| case | np | what it exercises |
|---|---|---|
| `cfd_2d_nwt` | 1 | default N40=3 (FC3), relaxation waves, level set, implicit diffusion, gauges/probes |
| `cfd_2d_nwt_rk3` / `_rk2` / `_fc2` / `_rkls3` | 1 | N40=13 / 12 / 2 / 44 momentum + level-set variants |
| `cfd_2d_nwt_wenoz` | 1 | WENO-Z weights for momentum (D 12 1) and level set + reini (F 37 1) |
| `cfd_2d_nwt_komega` / `_kepsilon` | 1 | RANS models with free surface |
| `cfd_2d_nwt_porous_komega` | 1 | VRANS porous box + k-ω porous sources (B295) |
| `cfd_2d_nug_breaking_komega` | 1 | stretched grid, cnoidal waves, k-ω T36=2, slope solid |
| `cfd_2d_dambreak` (+ `_fcc3`, `_fcls3`, `_mpi2`) | 1/2 | closed tank, walls, N40=33 / 4, 2D MPI |
| `cfd_2d_dambreak_teno` | 1 | TENO5 weights for momentum (D 12 2) and level set + reini (F 37 2) |
| `cfd_2d_dambreak_plic` (+ `_rk3`, `_fcc3`) | 1 | PLIC VOF (F80=4) with N40=3 / 13 / 33 |
| `cfd_2d_dambreak_picard`, `cfd_2d_dambreak_plic_picard_lsm` | 1 | Picard volume correction: F46=2 (once per step), F46=3 with PLIC redistancing (F92=1) |
| `cfd_2d_cylinder_singlephase` | 1 | single phase, explicit diffusion, inflow/outflow, forces |
| `cfd_3d_dambreak_obstacle` | 4 | 3D, MPI halos, solid box |
| `cfd_3d_pier_komega` (+ `_rkls3_sf`, `_t33`) | 2 | 3D inflow/outflow, k-ω wall functions, cylinder, N40=14 sf loop, T33 k-gradient source |
| `cfd_3d_heave_sphere_6dof` (+ `_rk3`) | 4 | floating body 6DOF (FCLS3), N40=13→14 df loop |
| `cfd_3d_heave_sphere_6dof_mooring` | 2 | 6DOF with spring mooring (X310=4), all DOFs free, rotational damping |
| `nhflow_2d_6dof_box` (+ `_rk3`) | 1 | NHFLOW floating box in waves, two-way 6DOF, A510=2 / 3 |
| `nhflow_2d_6dof_box_init` | 1 | + linear damping X25/X26, initial pitch X101, initial velocities X102/X103 |
| `nhflow_2d_6dof_box_towed` | 1 | towed in surge (X11 u=2, X210) with velocity ramp X206 |
| `nhflow_2d_6dof_box_oneway` | 1 | prescribed motion, one-way coupling (X10=2) |
| `nhflow_2d_6dof_membrane_collar` | 2 | rigid membrane bag (X330) on a floating collar |
| `nhflow_3d_membrane_cylcone_current` | 2 | fixed membrane bag `cylcone` (cylinder on a cone, sloped floor) in a current: B60=1 + potential start, HLLC (needs patch 0002) |
| `nhflow_3d_membrane_cylcone_drain` | 2 | flexible `cylcone` bag, filling + drain, staggered coupling, compression; hydrograph inflow B60=2 (needs patches 0001, 0002) |
| `nhflow_3d_membrane_sharp_static80` | 2 | fixed `cylcone` bag filled to 80 % in still water, `mobility sharp` (nhflow_thinbody: blocked links, wall fluxes, head below the floor, rigid lid) (needs patches 0004-0006) |
| `nhflow_3d_membrane_sharp_current` | 2 | fixed `cylcone` bag in a current with `mobility sharp`: B60=1 + potential start, guarded fluxes near the bag (needs patches 0004-0006) |
| `nhflow_3d_membrane_sharp_kw` | 2 | as `nhflow_3d_membrane_sharp_current` with k-ω (A560=2): fabric as wall of the diffusion and the turbulence model (needs patches 0004-0006) |
| `nhflow_3d_6dof_box` | 2 | NHFLOW 3D box, all six DOFs free, initial roll/yaw |
| `nhflow_3d_sphere_force` | 1 | NHFLOW force integration P 81 in still water: submerged sphere (A 586) on 6 sigma layers, marching-tetra quads split into two triangles: Fx = Fy = 0 to round-off, Fz = rho g V of the reconstructed surface (needs patch 0001 of NHFLOW Forces; reference with 0001) |
| `nhflow_3d_cylinder_force_mpi2` | 2 | P 81 on the vertical cylinder in the k-eps channel (base `nhflow_3d_cylinder_kepsilon_mpi2`), surface across the rank interface, globalsum of the forces (needs patch 0001 of NHFLOW Forces; reference with 0001) |
| `nhflow_3d_shipwave_box` | 1 | NHFLOW moving pressure patch, ship-wave mode (X10=3, X400=2) |
| `sflow_shipwave_box` (+ `sflow_6dof_box_oneway`) | 1 | SFLOW ship pressure patch (X10=3) / one-way direct forcing (X10=2) |
| `fnpf_2d_6dof_box` (+ `_rk4`) | 1 | FNPF resolved floating body, A310=3 / 4, added-mass coupling |
| `nhflow_3d_ship_box_thrust` | 2 | ship module (X350, ship.dat): constant thrust with offset, ITTC-1957 friction + form factor |
| `nhflow_3d_ship_box_6dof` | 2 | ship module, all six DOFs: cross-flow drag (strips), roll damping |
| `nhflow_3d_ship_box_propeller` | 2 | ship module propeller: KT/KQ, wake fraction, Hough-Ordway actuator disk as NHFLOW momentum source (all three components, swirl without net force), shaft torque reaction |
| `nhflow_3d_ship_box_zigzag` | 2 | ship module, surge/sway/heave/yaw: propeller with sampled inflow (NHFLOW velocity), MMG rudder, zig-zag steering with rudder rate |
| `nhflow_3d_ship_box_yawfree` | 2 | ship box with sway and yaw free, side faces through cell centres: level set = 0 on the surface, X 15 1 forcing ramp, no side force / yaw kick |
| `nhflow_3d_ship_kvlcc2_mmg` | 2 | ship module MMG model: KVLCC2 L7 hull derivatives, implicit added mass, MMG wake, rudder f_alpha / asymmetric gamma_R, fluid loads masked (mmg_fluid 0), zig-zag start |
| `nhflow_3d_ship_box_current` | 2 | ship held in a current: discharge inflow, hull turned by 180 deg (X 101), friction from the velocity relative to the current |
| `cfd_3d_ship_box_propeller` | 2 | ship module in CFD: box barge in a two-phase tank, actuator disk on the staggered velocity points (water part outside the hull, exact T and Q), velocity sampling, SSP-RK3 |
| `fnpf_3d_ship_box_hybrid` | 2 | ship module in FNPF: box in head waves, six DOFs, hybrid MMG (mmg_fluid 2, Munk moment out of N'v), approach phase from rest, propeller inflow sampled from the FNPF velocity, FNPF psi_0 = chi (body-following time derivative), X 17 relaxation of eta next to the footprint |
| `cfd_2d_fem_obstacle` (+ `_mpi2`) | 1/2 | FEM solid (Z 30, N10=1): elastic obstacle hit by the bore, coupling across a subdomain border |
| `cfd_2d_fem_wall_failure` | 1 | FEM concrete wall cracking, erosion, the broken wall becomes a rigid fragment (`fragments rigid`), debris, ground and part contact |
| `cfd_3d_fem_column` | 4 | FEM solid in 3D on 4 ranks: elastic column hit by the bore |
| `cfd_2d_cpm_bedload_layer` (+ `_susp`, `_mpi2`) | 1/2 | CPM MP-PIC sand bed in a periodic channel with the sub-grid bedload layer (Q 58 1 / 2): pickup, hops, deposition, release into the suspension, column bed level, one-sided solid forcing at the particle bed, CPM log |
| `cfd_2d_fem_simple_wall` | 1 | FEM simple input (concrete C30 preset, fix base, monitor auto, resolution), settling with the initial water, hybrid loads, structural damping |
| `cfd_2d_fem_floating` | 1 | FEM rigid body (floating box, `material rigid 500`) in the collapsing water column: rigid-body integration, pressure loads with the probe correction and the added-mass stabilisation |
| `cfd_2d_fem_debris_impact` | 1 | FEM rigid timber block (`stiffness 2e5`) on the dry bed pushed by the dam-break surge onto an elastic post: debris impact as one spring of the debris stiffness, ground contact, pressure beside the body for faces on the floor (reference with patch 0021: the floor of the domain is no second ground under `ground 0.0`) |
| `nhflow_2d_fem_wall` | 1 | FEM coupled with NHFLOW (Z 30, A 10 5): concrete wall holding a 0.3 m water column on the dry bed; cells inside the structure marked solid (p->DF -1, no continuity flux, d->solid_flux), water at rest, hydrostatic load (needs patches 0016, 0017, 0019; references with 0019) |
| `nhflow_2d_fem_floating` | 1 | FEM rigid floating box in NHFLOW still water: probed pressure (filtered non-hydrostatic part plus hydrostatic of the free surface), added-mass factor 5, heave towards the 5 cm draft (needs patches 0016, 0017, 0019; references with 0019) |
| `nhflow_2d_fem_debris_impact` (+ `_mpi2`) | 1/2 | FEM rigid timber block carried by an NHFLOW dam-break surge over the dry bed onto a concrete post that stands dry: wet/dry columns, debris contact, probes kept out of solid cells, subdomain border at the post (needs patches 0016–0021; references with 0021) |
| `nhflow_2d_fem_light_debris` | 1 | FEM empty 20 ft container (rigid, 59.3 kg/m3, 2D slice) dropped 5 m onto 4 m of still NHFLOW water: added-mass factor limited to 30 body masses, speed relative to the water at most the terminal speed (needs patches 0016–0021) |
| `cfd_2d_fem_rc_wall` | 1 | FEM reinforced concrete (`rebar`): weak concrete wall with 1 % model steel (f_y 30 MPa) holding the dam-break column; cracks at the base, the bars yield, the wall bends 7.5 mm and stands; unilateral damage, steel utilisation and first yield in the summary (needs patch 0022) |
| `nhflow_2d_nwt_stokes5` (+ `_mpi2`) | 1/2 | NHFLOW relaxation generation + beach (B98=2/B99=1), Stokes 5th |
| `nhflow_2d_nwt_stokes5_b95` | 1 | the same with the step-size consistent relaxation zones (B 95 0.01): targets at the stage output times, factor r^(dt/0.01) |
| `nhflow_2d_dirichlet` | 1 | NHFLOW Dirichlet wave generation (B98=3) |
| `nhflow_2d_awa` | 1 | NHFLOW active wave generation + active absorption (B98=4/B99=3) |
| `nhflow_2d_current` | 1 | NHFLOW waves on a current: discharge inflow (B60=1, W10), beach relaxes to the current (B97=1) |
| `nhflow_2d_irregular_decomp` | 1 | NHFLOW JONSWAP waves, decomposed relaxation precalc (B89=1) |
| `nhflow_3d_irregular_decomp` | 2 | NHFLOW 3D short-crested waves (B130=2), decomposed precalc, CFL 0.5 |
| `fnpf_2d_regular` (+ `_mpi2`) | 1/2 | FNPF relaxation generation + beach, regular waves |
| `fnpf_2d_irregular_decomp` | 1 | FNPF JONSWAP waves, decomposed relaxation precalc (B89=1) |
| `fnpf_3d_shortcrested` | 2 | FNPF 3D short-crested waves (B130=2) |
| `nhflow_2d_two_sources` (+ `_mpi2`) | 1/2 | wave_field: B 92 linear wave + linear source 2 with phase (B 500/501); P 50 prints the summed target |
| `nhflow_2d_irregular_two_sources` | 1 | wave_field: two JONSWAP sources with their own seeds (B 139, B 504) |
| `fnpf_2d_two_sources` | 1 | wave_field in FNPF (cached potential path) |
| `nhflow_3d_amr_relax` | 2 | NHFLOW 3D static AMR, refinement box over the domain; relaxation zones stay unrefined (bc_zone boxes) |
| `fnpf_3d_amr_relax` | 2 | FNPF 3D static AMR, same check |
| `nhflow_3d_amr_still` | 1 | NHFLOW AMR static patch in a closed tank, bed ramps under the patch edges: still water stays still (fill, flux matching, restriction, composite pressure BiCGStab + FAC) |
| `nhflow_3d_amr_waves_mpi2` | 2 | NHFLOW AMR regular waves through a static patch split at the partition edge: patches on two ranks, remote fills and face runs |
| `nhflow_3d_amr_heave_mpi2` | 2 | NHFLOW AMR heave decay of a floating box (X 10 1), body zone G 12 across the partition edge: forcing on every grid, loads from the finest |
| `nhflow_3d_amr_tow_zone` (+ `_vref_rk3`) | 1 | NHFLOW AMR towed box (X 10 2) with moving zone and wake wedge (G 12/279), regrid every step: from_old, prolongation, patch deletion; `_vref_rk3`: vertical refinement G 6, RK3, A 520 1 |
| `nhflow_3d_amr_adaptive` | 1 | NHFLOW AMR solution-adaptive (G 20, G 21), ring wave from an F 72 hump, regrid every step: tagging, F 72 boxes on fresh patches |
| `nhflow_3d_amr_wetdry` | 1 | NHFLOW AMR wetting and drying in the patches (G 30): beach and cone island, a box crossing the island shoreline, run-up at t = 0 (wet flags of the cells around a patch, lexer::wetfix, covered-cell flags) |
| `nhflow_3d_amr_wetdry_adaptive` (+ `_mpi2`) | 1/2 | the same with the shoreline flag (G 22) only, regrid every step, layout changes from step 116; `_mpi2`: partition edge through the island |
| `nhflow_3d_amr_wetdry_n14d` | 1 | `nhflow_3d_amr_wetdry` with the V-cycles of the AMR pressure preconditioner in double (N 14 64; the default N 14 32 keeps them in float) |
| `nhflow_3d_amr_waves_place` | 2 | `nhflow_3d_amr_waves_mpi2` with the box on rank 1 only and the patches placed for the load of the ranks (G 40 1): pieces on rank 0 reach their parents on rank 1 through the block plans |
| `nhflow_3d_amr_heave_place` | 2 | floating body with placed patches (G 40 1): hull triangles taken by the rank of the finest grid at their centroid |
| `nhflow_3d_amr_tow_place2` | 2 | towed box with moving zone, placement test mode (G 40 2): patches migrate to other ranks at every regrid |
| `nhflow_3d_amr_wetdry_place2` | 2 | adaptive wetting/drying, G 40 2: shoreline wet flags restricted and prolonged across ranks |
| `nhflow_3d_amr_breaking_place2` | 2 | adaptive breaking (A 550 + A 512 2), G 40 2 |
| `nhflow_3d_amr_breaking` | 1 | NHFLOW AMR wave breaking (A 550 1, A 512 2): a bore runs into a static patch; the composite implicit diffusion over all grids (G 31 1, default), filled cells take the source grid's breaking (lexer::amrvb) |
| `nhflow_3d_amr_breaking_adaptive` (+ `_mpi2`) | 1/2 | breaking with the adaptive patch following the bore (G 20, G 23), wetting and drying, regrid every step; `_mpi2`: different patch counts per rank (local reductions, gcparaxijk_single) |
| `nhflow_3d_amr_breaking_g31_0` | 1 | `nhflow_3d_amr_breaking_adaptive` with the per-grid implicit diffusion of step 7 (G 31 0): a patch-local bicgstab_ijk per patch |
| `nhflow_3d_amr_cdiff_vref_rk3` | 1 | `nhflow_3d_amr_tow_vref_rk3` with A 512 2: the composite implicit diffusion (G 31 1) with RK3 (phase_D), vertical refinement (G 6) and a regrid every step |
| `nhflow_3d_amr_still_sub` | 1 | still water over the bed ramps with subcycling (G 7 1): level pressure solve, synchronisation projection |
| `nhflow_3d_amr_adaptive_sub` | 1 | ring wave with G 7 1, a regrid every level-0 step: flux registers, refluxing, a level of several patches |
| `nhflow_3d_amr_waves_place_sub` | 2 | placed patches (G 40 1) with G 7 1 and B 95 -1: remote fills in time, remote refluxing, the level solve across ranks |
| `nhflow_3d_amr_breaking_sub` | 1 | G 7 1 with RK2, breaking (A 550 1, per-grid diffusion) and wetting and drying in the patches (G 30 1) |
| `nhflow_3d_amr_heave_sub` (+ `heave_place_sub`) | 2 | heave decay with subcycling (G 7 1): the finest level advances the body with its stages and loads, level 0 a predicted copy, the synchronisation leaves the forcing band alone; `heave_place_sub` with placed patches (G 40 1) |
| `nhflow_3d_amr_tow_zone_sub` | 1 | towed box (X 10 2) with G 7 1: prescribed motion on every level at its own step, the moving zone regridded every level-0 step |
| `sflow_2d_amr_dambreak` (+ `_mpi2`, `_place`, `_place2`) | 1/2 | SFLOW AMR dam break (tutorial 8_6): 2 levels following the surface jump (G 20), regrid every 4 steps; `_mpi2` patches cut at the rank boxes, `_place` G 40 1, `_place2` G 40 2 (remote parents, migration: restrict_levels, ini_patch_old/prolong) |
| `sflow_2d_amr_bar_nh` (+ `_place2`) | 1/2 | SFLOW AMR non-hydrostatic bar (tutorial 8_7, A 220 1): static boxes, composite pressure; `_place2` through the block plans |
| `sflow_2d_amr_bar_bous` (+ `_place2`) | 1/2 | the bar with Boussinesq (A 220 4): u_a on the leaf cells of all levels; `_place2` G 40 2 |
| `sflow_2d_amr_dambreak_sub` (+ `_place2`) | 1/2 | the dam break with subcycling (G 7 1): level l takes 2^l steps per level-0 step, fills interpolated in time, flux registers and refluxing, regrid every level-0 step; `_place2` G 40 2 (remote fills and refluxing) |
| `sflow_2d_amr_bar_sub` | 1 | hydrostatic waves over the bar (A 220 0), static boxes, G 7 1 |
| `sflow_2d_amr_bar_sub_b95` | 1 | the subcycled bar with B 95 -1 (reference step of the finest level) |
| `fnpf_3d_amr_basin` (+ `_mpi2`, `_place2`) | 1/2 | FNPF AMR waves through a static box (composite Laplace); `_place2` G 40 2: remote parents in the FAC preconditioner |
| `fnpf_3d_amr_tow` (+ `_place2`) | 1/2 | FNPF AMR towed floating body (X 10 1, X 210), the zone follows it (G 12, G 2 4); `_place2` G 40 2: zone patch placed whole, body grids and triangle ownership |
| `nhflow_3d_two_edges` | 2 | zones with own sources (B 520/521/524): x- zone generates the B 92 wave, y- zone source 2 at 90 deg; beach zone from B 520 |
| `fnpf_3d_two_edges` | 2 | the same in FNPF |
| `nhflow_2d_custom_zones` | 1 | old input: custom B 108 generation zone, two B 107 beach zones |
| `nhflow_2d_origin` | 1 | old input: generation origin B 105 |
| `fnpf_2d_custom_zones` | 1 | old input: custom B 107 / B 108 zones in FNPF |
| `cfd_2d_nwt_custom_zones` | 1 | old input: custom B 107 / B 108 zones in CFD |
| `nhflow_3d_amr_custom_zones` | 2 | old input with AMR: custom zones wider than B 96; only the B 96 ranges stay unrefined, as before |
| `nhflow_2d_dirichlet_custom_zones` | 1 | old input: Dirichlet generation with two B 107 beach zones |
| `nhflow_2d_awa_custom_zones` | 1 | old input: active generation (B 98 4) with a B 107 relaxation beach, waves reach the beach |
| `fnpf_2d_dirichlet` | 1 | old input: FNPF Dirichlet generation with a B 107 beach |
| `fnpf_2d_dirichlet_awa` | 1 | old input: FNPF Dirichlet generation with active absorption (B 99 3), waves reach the beach |
| `sflow_2d_stokes2` (+ `_mpi2`) | 1/2 | SFLOW relaxation generation + beach, Stokes 2nd (tutorial 8_1) |
| `sflow_2d_stokes2_wenojs` | 1 | SFLOW with the previous WENO-JS weights (A 213 0; default WENO-Z) |
| `sflow_2d_custom_zones` | 1 | old input: SFLOW with custom B 107 / B 108 zones |
| `sflow_2d_dirichlet` | 1 | old input: SFLOW Dirichlet generation |
| `fnpf_2d_hdc` (stage `fnpf_hdc_source`) | 1 | old input: hydrodynamic coupling FNPF -> FNPF (FNPF state P 44 3, DIVEMesh H 10 44, B 92 61), as in FNPF nesting |
| `fnpf_2d_hdc_single` (stage `fnpf_hdc_source_single`) | 1 | the same with one state file per step (P 45 1); same result as `fnpf_2d_hdc` |
| `nhflow_2d_hdc` (+ `_mpi2`, stage `fnpf_hdc_source_uvw`) | 1/2 | FNPF -> NHFLOW nesting: FNPF state with velocities (P 44 1), DIVEMesh H 10 4, NHFLOW reads eta, u, w (B 92 61) |
| `nhflow_2d_hdc_single` (stage `fnpf_hdc_source_uvw_single`) | 1 | the same with one state file per step (P 45 1); same result as `nhflow_2d_hdc` |
| `nhflow_2d_hdc_nhflow` (stage `nhflow_hdc_source`) | 1 | NHFLOW -> NHFLOW nesting (DIVEMesh H 10 5); target time step close to the source state interval |
| `nhflow_2d_tide_flather` (+ `_mpi2`) | 1/2 | tidal background (B 510/511): progressive tide in through a Riemann edge (B 520 method 3), out through a Flather edge (method 4) |
| `nhflow_2d_tide_riemann` | 1 | tidal channel with Riemann edges at both ends |
| `nhflow_2d_tide_basin` | 1 | closed basin: Riemann edge at x-, wall at x+; the reflected tide leaves through the Riemann edge |
| `nhflow_2d_tide_waves_beach` | 1 | tide + waves: generation zone with background + waves, beach relaxing to the tide (B 523) |
| `nhflow_2d_current_background` | 1 | constant current background (B 514): Riemann in, Flather out, beach relaxing to the current |
| `nhflow_2d_tide_progressive` | 1 | one progressive background (B 515) for the Riemann edges at both ends |
| `nhflow_2d_tide_timeseries` | 1 | time-series background (B 510 mode 2, `background-1.dat`): Riemann in, Flather out |
| `nhflow_2d_riemann_waves` | 1 | waves from a Riemann edge with a wave source (B 524), no relaxation zone; beach at x+ |
| `nhflow_2d_riemann_waves_wall` | 1 | waves from a Riemann edge in a channel closed by a wall; the reflected waves leave through the edge |
| `nhflow_2d_tide_waves_riemann` | 1 | tide + waves through a Riemann edge (background + source), no generation zone |
| `nhflow_3d_tide_ychannel` (+ `_mpi2`) | 1/2 | tidal channel along y: Riemann edge at y- (B 520 method 3), Flather edge at y+ (method 4), progressive background at 90 deg; `_mpi2` splits the ranks in y |
| `nhflow_3d_tide_oblique` (+ `_mpi2`) | 1/2 | oblique tide (30 deg) in a basin with open edges on all four sides: Riemann at x- and y-, Flather at x+ and y+; `_mpi2` splits the ranks in x |
| `nhflow_2d_waves_current_doppler` (+ `_mpi2`) | 1/2 | linear waves on a following current (B 514) with Doppler (B 530 2 10): k from omega = sigma + k U, blended during the current spin-up |
| `nhflow_2d_waves_setup_heff` | 1 | linear waves on a 1 m set-up (B 514 eta) with B 530 1 10: k on h_eff |
| `nhflow_2d_tide_waves_doppler` | 1 | `nhflow_2d_tide_waves_beach` with B 530 2 10: k follows the tide level and current |
| `cfd_2d_channel_kepsilon` (+ `cfd_2d_channel_komega_mpi2`) | 1/2 | open channel, discharge inflow (B60 1) with the equilibrium k/ε/ω inflow profile, k-ε / k-ω across a rank border in x |
| `cfd_2d_channel_komega_t36` | 1 | k-ω free-surface damping T36 3 (y' = T37 h from the local water depth, dimensionless weight) |
| `cfd_2d_stillwater_plic_t41` | 1 | PLIC VOF still water, k-ω with T41 1: no NaN from the limiter at S = 0 |
| `nhflow_2d_channel_kepsilon` | 1 | NHFLOW open channel, discharge inflow with the equilibrium turbulence profile, k-ε, bed roughness A519 |
| `nhflow_2d_channel_suspended` | 1 | NHFLOW open channel with suspended load only (S 11 0, S 12 1), k-ε: conservative D·C scheme with FEx face fluxes, sediment leaving through the outflow, bed exchange only in the erodible region S 71 (needs Suspended Sediment patches 0001-0009) |
| `nhflow_2d_beach_suspended_mpi2` | 2 | NHFLOW waves on a 1:12 sand beach, bedload + suspended load, swash wetting/drying (dry-column deposit, Exner with the c_be the water column used), S 38, FEx across the MPI interface (needs Suspended Sediment patches 0001-0009) |
| `nhflow_3d_cylinder_kepsilon_mpi2` | 2 | NHFLOW 3D channel with a cylinder (A580), k-ε, ranks split in y: k/ε and ν_t across the rank interface |
| `sflow_1d_channel_ke` (+ `_kw`) | 1 | SFLOW depth-averaged k-ε / k-ω (A260 1/2): k, ε/ω relax to the Rastogi–Rodi equilibrium |
| `sflow_2d_channel_walls_kw_mpi2` | 2 | SFLOW k-ω with side walls, ranks split in y: production at wall cells, uniform k/ω across the width |
| `cfd_3d_heave_sphere_6dof_kepsilon` | 4 | k-ε at a direct-forcing body: wall functions at gcdf4 cells, flagsf4 turn-off, solid forcing for k/ε |
| `cfd_2d_channel_veg_kepsilon` | 1 | vegetation box (B 310, B 308 0): drag ½ Cd a \|u\| u_i, Lopez & Garcia k/ε sources |
| `cfd_2d_freesurface_komega_t45` | 1 | k-ω buoyancy term T 45 1: implicit sink, k ≥ 0 at the interface |
| `cfd_2d_channel_les_t21_2` | 1 | LES Smagorinsky with the second-order high-pass filter T 21 2 |
| `nhflow_3d_channel_walls_kepsilon_mpi2` | 2 | NHFLOW side walls with A519 2: friction and turbulence wall functions per wall face (nhflow_wall.h), ranks split in y |
| `seastate_2d_init` (+ `_mpi2`) | 1/2 | REEF3D::SEASTATE (A 10 7): shoal + land block, JONSWAP initial spectrum with Mitsuyasu spreading (A 710 1), block-sparse storage (8 x 8 tiles), no boundary spectrum: the waves leave the domain through its open sides |
| `seastate_2d_shoaling` | 1 | REEF3D::SEASTATE stationary (A 700 2) shoaling on a plane slope, unidirectional JONSWAP from x- (A 711 1); validated against energy-flux conservation (Dropbox SEASTATE/validation/01_shoaling) |
| `seastate_2d_refraction` (+ `_mpi2`) | 1/2 | REEF3D::SEASTATE stationary refraction at 30 deg on parallel contours, boundary spectrum on x- and y-, 72 directions; 2 ranks give the 1-rank result (halo exchange `seastate_exchange`) |
| `seastate_2d_front` (+ `_mpi2`) | 1/2 | REEF3D::SEASTATE nonstationary (A 700 1) front in a 20 km channel, dt 60 s (Courant ~8): implicit M-matrix, N >= 0, energy conservation; 2 ranks with 2 iterations per step |
| `seastate_2d_current` | 1 | REEF3D::SEASTATE opposing current (A 720 1, A 721): frequency shift c_sigma and Doppler term |
| `seastate_2d_spc` | 1 | REEF3D::SEASTATE boundary spectrum from a SWAN 2D spectral file (A 711 2, `seastate-boundary.spc`, nautical directions) |
| `seastate_2d_wind` | 1 | REEF3D::SEASTATE stationary fetch-limited wind sea (A 730-735): Komen wind input and whitecapping, DIA, action density limiter, pseudo time step, zero-gradient lateral sides (A 712 2); validated in Dropbox SEASTATE/validation/05_fetch_limited |
| `seastate_2d_windsea_mpi2` | 2 | REEF3D::SEASTATE nonstationary wind-sea growth from rest, oblique wind, 2 ranks: source terms with the limiter and the halo exchange |
| `seastate_2d_surf` | 1 | REEF3D::SEASTATE stationary surf zone on a plane slope: Battjes-Janssen breaking (A 740), JONSWAP friction (A 742), LTA triads (A 744), oblique Mitsuyasu spectrum, zero-gradient sides; single terms validated in Dropbox SEASTATE/validation/06-08 |
| `seastate_2d_handover` | 1 | REEF3D::SEASTATE handover of 2D spectra for the FNPF/NHFLOW wave generation (A 760/A 761, `REEF3D_SEASTATE_Spectra/spectrum-file-2d_P<n>.dat`, B 85 11 format), incl. a point outside the domain; validated in Dropbox SEASTATE/validation/11_handover |
| `seastate_sflow_setup` | 1 | REEF3D::SEASTATE coupled with SFLOW (A 10 2, A 750 1), radiation stress (A 751 1): set-down and set-up on a plane beach in a closed basin, water level and currents fed back, drying cells; validated in Dropbox SEASTATE/validation/09_setup |
| `seastate_sflow_vf_mpi2` | 2 | REEF3D::SEASTATE coupled with SFLOW on 2 ranks, vortex force (A 751 2) with the Stokes transport in the SFLOW continuity, oblique Mitsuyasu spectrum, handover point in a coupled run |
| `seastate_sflow_surfbeat` | 1 | REEF3D::SEASTATE surfbeat (A 770 1) coupled with SFLOW: wave groups (JONSWAP Hs 0.2 m, Tp 2.25 s) on the GLOBEX profile, Roelvink breaking (A 740 2), roller (A 748), radiation stress with roller, bound long wave in and absorbing long-wave boundary at x- (A 774), Crank-Nicolson and second-order group transport (A 775 2); validated in Dropbox SEASTATE/validation/12-14 |
| `seastate_sflow_surfbeat_vf_mpi2` | 2 | REEF3D::SEASTATE surfbeat coupled with SFLOW on 2 ranks, vortex force with roller: short-crested groups at 20 deg in a 90 m x 20 m basin, boundary rows along y, halo exchange with 2 layers |
| `seastate_nhflow_surfbeat` | 1 | REEF3D::SEASTATE surfbeat coupled with NHFLOW (A 10 5, A 750 1), 2 sigma layers: wave force in every layer, long-wave ghost cells at the open x- edge; NHFLOW vs. SFLOW in Dropbox SEASTATE/validation/13_globex |
| `seastate_2d_surfbeat` | 1 | REEF3D::SEASTATE surfbeat stand-alone (A 10 7): short-crested groups on the GLOBEX profile, T_rep 2 s (A 771), Roelvink breaking and roller |
| `seastate_2d_boundary_series` | 1 | REEF3D::SEASTATE time series of boundary spectra at many locations (A 711 3, nonstationary SWAN file): 3 locations, 4 hourly records, spectra on the x- and y- sides from the two nearest locations, start 30 min after the first record (A 780); validated in Dropbox SEASTATE/validation/15-16 |
| `seastate_2d_wind_field_mpi2` | 2 | REEF3D::SEASTATE wind field in space and time (A 730 2, seastate-wind.dat) with Komen physics and DIA, boundary spectra from a stationary 2-location file, 2 ranks; validated in Dropbox SEASTATE/validation/17-18 |
| `seastate_2d_forcing_stationary` | 1 | REEF3D::SEASTATE quasi-stationary sequence (A 700 2, one stationary solve per A 706) with boundary spectra changing in time and Battjes-Janssen breaking |
| `seastate_2d_amr_static` | 1 | REEF3D::SEASTATE static mesh refinement (G 1 2, REEFAMR): island on a slope from the bathymetry raster (A 790 1, all grids), criteria depth/coastline/gradient (A 791-793), box G 10, no refinement next to the forced side (A 794), stationary composite solve with whitecapping, breaking with the SWAN maximum energy (A 737), triads with SWAN 41.31 parameters (A 736), SWAN-type convergence (A 709), handover on the patches; validated in Dropbox SEASTATE/validation/19-23 |
| `seastate_2d_amr_mpi2` | 2 | REEF3D::SEASTATE mesh refinement on 2 ranks (G 1 1): the island case nonstationary from rest with wind (Komen, DIA), 2 iterations per step, patches cut at the rank boxes, fills and restriction across ranks |
| `seastate_2d_bathy_raster` | 1 | REEF3D::SEASTATE bathymetry raster without refinement (A 790 1, majority of wet/dry raster nodes per cell) with A 736, A 737, A 709 |
| `seastate_2d_phase7_sparse` | 1 | REEF3D::SEASTATE spectral sparsity (A 795 1e-6, Phase 7) on the refraction case: per cell, frequency and quadrant only the window of bins that can hold energy is solved; validated in Dropbox SEASTATE/validation/25, 26 |
| `seastate_2d_phase7_sordup_mpi2` | 2 | REEF3D::SEASTATE second-order upwind geographic fluxes (A 796 2, SORDUP type, Phase 7) on the refraction case, two halo layers; validated in Dropbox SEASTATE/validation/27 |
| `seastate_2d_phase7_composite` | 1 | REEF3D::SEASTATE composite sweep across the refinement levels (A 797 1, Phase 7) on the island case with A 795 and A 796 2; validated in Dropbox SEASTATE/validation/28 |
| `seastate_2d_phase7_composite_mpi2` | 2 | REEF3D::SEASTATE composite sweep (A 797 1) on the nonstationary island case with wind and DIA, patches cut at the rank boxes |
| `sflow_2d_bank_ediff` (+ `_idiff`) | 1 | SFLOW sloping bank (T 62) with a shoreline, constant viscosity, A212 1 / 2: free slip at the dry neighbours |
| `cfd_2d_channel_kepsilon_stretched_ifou` | 1 | CFD open channel on a grid stretched in x (B 101 1, B 111 2.0), k-ε with implicit first-order upwind T 12 1: per-face conservative upwind |
| `nhflow_2d_bump_ediff` | 1 | NHFLOW 2D flow over a bump, constant viscosity, explicit momentum diffusion A 512 1: σ face metrics, A 513 wall rule (no no-slip bed from the A 518 2 ghost), viscous time-step limit |
| `nhflow_3d_bank_kepsilon` | 1 | NHFLOW 3D channel with a sloping bank (T 62) and a shoreline along x, k-ε: zero gradient of k and ε towards dry columns |
| `cfd_2d_patch_inflow_komega` | 1 | CFD channel with a discharge inflow patch (B 441, B 411) instead of B 60, k-ω: IO flags at the patch faces, turbulence zero gradient there |
| `nhflow_2d_bore_beach_kepsilon` | 1 | NHFLOW 2D bore running up a dry 1:20 beach, k-ε: newly wetted columns start from the wet neighbours, length-scale limit l ≤ κh |
| `nhflow_3d_boom_particles_mpi2` | 2 | NHFLOW river plastic interceptor: angled boom (B 350, skirt depth below the moving free surface, normal-only resistance), catamaran hulls and a permeable conveyor as surface-piercing drag boxes (B 351), floating and buoyant particles blocked above the draft and passing under it (L 71), capture zone at the conveyor (L 72), boom/hull force output (REEF3D_NHFLOW_Boom); particles cross the rank interface along the boom (needs River Interceptor patches 0001-0002) |

Tag `quick` selects a subset that runs in a few minutes. Adding a case: copy a directory, edit,
run `./regression.py run ... --cases <new>`, check it, then `bless`.

## Notes

- Cases tagged `legacy` cover old input that must keep giving the same results as before
  the iowave redesign; A/B them against a binary of the original code after every change.
- Irregular-wave cases must fix the random seeds (`B 139`, and `B 138` for directional
  spreading); otherwise the phases come from `srand(time(0))` and no two runs agree.

- The stored references were blessed on hans_dev bacfbc4a2 (Linux, gcc, MPICH, `build_reef3d.sh`
  flags). `check` uses a tolerance, but another compiler or platform can still differ in the last
  digits; then bless once on your machine and keep using `ab` for bitwise checks.
  The NHFLOW AMR cases of steps 5-7 (`nhflow_3d_amr_still` ... `_breaking_mpi2`) were blessed
  on hans_dev 925e271a7 + REEFAMR patches 0001/0002 (patch lexer periodic flags, Z11 default).
- The HDC cases (`*_hdc*`) need DIVEMesh 95ae3a4 or later (state file names in
  `hdc_filename_in.cpp`); `nhflow_2d_hdc_nhflow` also needs the pressure skip in
  `hdc_read_nhflow.cpp` (NHFLOW state files with the default P 44 1).

- Runs use few steps on coarse grids: they test that code paths give the same numbers, not that
  the physics is right. Cases with an analytical or reference solution are in `Validation/`.
- A run fails on NaN/Inf in the per-step norms, a non-zero exit code, a timeout or a missing dump.

## Unit tests

`unit/` holds standalone tests of solver-independent kernels (no MPI, no REEF3D binary); the build
line is at the top of each file, run from `unit/`:

| test | covers |
|---|---|
| `rigidbody_test.cpp` | 6DOF rigid-body core: quaternion/Euler, constant force, oscillator and torque-free top (orders of RK2/RK3/RKLS3/RK4), DOF modes, damping |
| `ship_test.cpp` | ship module kernels: waterline clipping, wetted surface, draft strips, cross-flow drag, ITTC-1957 line, roll damping, propeller KT/KQ, actuator disk (discrete force and torque, swirl sense), MMG rudder (signs, symmetry, slipstream, course stability, f_alpha, asymmetric gamma_R), MMG hull polynomials, MMG wake |
| `fem_test.cpp` | FEM solid solver: cantilever (Timoshenko, frequency), objectivity, J2, crack band energy, contact, collapse, STL snapping, patch test, presets, settling/check, structural damping, rigid bodies, walls and inclined bed (stick/slide), debris impact (peak u √(kM), duration, restitution, crushing, two bodies in series, deformable wall, added mass not in the impact, crushing set per contact point), collapse with a rigid fragment, reinforced concrete (bars in the skin, reinforced tie: ρσ_s and rupture, RC cantilever pushover against a section analysis, plain concrete brittle) |
| `lagoon_store_test.cpp` | LAGOON store writer (P 18, needs `-lz`): VTU header parsing, σ-level offsets, shard files read back (append, CRC-32C index, inner chunks, components), Cartesian (CFD) blocks split in z |
| `lagoon_bodies_test.cpp` | LAGOON body writer (P 18, needs `-lz`; built with `../../src/lagoon_store.cpp`): quaternion of REEF3D's rotation matrix, two rigid bodies over 8 outputs (set and body attributes, mesh once, motion arrays, time counted once every body has it), a body off its rigid motion refused |
| `lagoon_particles_test.cpp` | LAGOON particle writer (P 18, needs `-lz`; built with `../../src/lagoon_store.cpp`): 5 outputs of 0 to 700000 particles (inner chunks and shards crossed), field names and their arrays, int32 fields, values read back from the shard files; objects with cells (ice floes: a cell set stored only when the cells change, a value per cell), inconsistent cells refused (needs patch 0009) |
| `lagoon_amr_test.cpp` | LAGOON AMR writer (P 18, needs `-lz`; built with `../../src/lagoon_store.cpp`): 4 outputs of 3 and 4 grids (level 0 of two ranks, a moving patch, a second one from the second output), set attributes and an int32 field, format version 0.3, grid_end and the grid table, coordinates as float64, levels, coordinates and cell values read back from the shard files; grids packed for MPI and back, a cut buffer refused; a grid with too few values refused and the set stopped (needs patch 0013) |
| `seastate_test.cpp` | REEF3D::SEASTATE kernels: spectral grid (logarithmic frequencies, bin widths, directions, direction faces, quadrants), block-sparse action storage (land tiles not allocated, tiles clipped at the range end, ghost cells, no overlap), integrated parameters (Hs, Tp, Tm01, Tm-1,0, direction, spread), dispersion relation (deep/shallow limits, residual, cg), SWAN 2D spectral file reader (CDIR/NDIR, VaDens/EnDens), memory budget, source terms (Battjes-Janssen Q_b, breaking with Newton linearisation, JONSWAP friction, moments with the spectral tail, Komen whitecapping, wind input with the Wu drag, DIA energy and action conservation, symmetry, shallow-water scaling and diagonal derivative, LTA energy conservation and Ursell limit, P >= 0 and D >= 0), surfbeat (single-frequency grid, Roelvink breaking, Herbers bound-wave coefficient vs. the LHS62 shallow-water limit, bichromatic envelope and bound wave, random-phase generator: variance, T_rep, directions, envelope mean, determinism), forcing files (date-times, series of SWAN spectra with many locations and times, wind field interpolation), bathymetry raster (reading, bilinear interpolation, cell mean and wet/dry majority rule), SWAN maximum energy (cap to (gamma d)^2/4, shape kept) |

