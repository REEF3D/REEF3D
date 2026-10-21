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
3. the normal text output (wave gauges, probes, forces, 6DOF, CPM sediment log, NHFLOW particle log, NHFLOW boom forces, DEM state and fixed-particle loads) with a tolerance.

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
normally removes (state files for such a stage). A restart from checkpoints is such a chain: the stage
writes them (P 48), the case copies `REEF3D_Checkpoint` (and the output folders, which the restart cuts back to the
checkpoint) and sets I 45 1 or I 45 <count>. `"rebalance": n` starts a run
again after a rebalance stop (G 46: exit code 3 and the marker `REEF3D_Checkpoint/rebalance`), at most
n times, as `tools/reef3d_rebalance.sh` does. The load maps (`*_loadmap.dat`) and the reports of the
SEASTATE pipeline (`*_pipeline.dat`, Phase 19), measured times, are not compared.

Other commands:

```bash
./regression.py list [--tags quick]                     # cases, tags, what they cover
./regression.py run  --reef3d BIN --divemesh DM --out DIR [--cases 'cfd_2d_*'] [--tags ...] [--steps N] [--keep]
./regression.py compare REF_DIR NEW_DIR [--require identical|close|different] [--rtol 1e-6] [--atol 1e-12]
./regression.py write --out DIR --cases cfd_2d_nwt      # only write the resolved input files
./regression.py bless DIR [--note "why"]                # store compact references in cases/<case>/reference.json
./regression.py check DIR [--rtol 1e-6]                 # check a run against the stored references
./regression.py ranks DIR                               # SEASTATE: cases with "same_as" against that case of the run, cell by cell
                                                        # (bitwise; with "same_as_rtol": Hs and the spectra within that tolerance);
                                                        # other models (a restart against the run it continues): the final state
                                                        # rank by rank and the text outputs, bitwise or within "same_as_rtol"
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
 "same_as": "<case>",                  // SEASTATE: bitwise the result of that case (other ranks); checked by `ranks`, not inherited
 "same_as_rtol": 1e-5,                 // with same_as: converged runs whose iterations differ (lagged halos): Hs relative to
                                       // itself and the bins relative to the largest bin of the cell within this, not bitwise
 "same_as_skip": ["REEF3D_FNPF_AMR/*_log.dat"], // with same_as, a restart (not SEASTATE): text outputs not compared
 "rebalance": 3,                       // restarts after a rebalance stop (G 46), at most this many
 "description": "...", "covers": ["code paths this case is meant to exercise"]
}
```

The runner always sets `M 10` (ranks) in both files, `N 45` (steps), `N 41 1e9`, and removes the
VTU/state print keys (`P 20/30/40/41/42`), so runs are short and output stays small.

| case | np | what it exercises |
|---|---|---|
| `cfd_2d_nwt` | 1 | default N40=3 (FC3), relaxation waves, level set, implicit diffusion, gauges/probes |
| `cfd_2d_relax_ramp_ls` | 1 | `cfd_2d_nwt` with level set RK3 (N40=13) and the start-up ramp B 101 1: the ramp acts on the wave elevation only (phi_relax, fixed in fe03b66a) |
| `cfd_2d_relax_ramp_plic` | 1 | `cfd_2d_nwt` with PLIC VOF (F80=4) and the start-up ramp B 101 1: the ramp acts on the wave elevation only, not on the still water level of zone and beach (vof_relax, fe03b66a) |
| `cfd_2d_nwt_rk3` / `_rk2` / `_fc2` / `_rkls3` | 1 | N40=13 / 12 / 2 / 44 momentum + level-set variants |
| `cfd_2d_nwt_wenoz` | 1 | WENO-Z weights for momentum (D 12 1) and level set + reini (F 37 1) |
| `cfd_2d_nwt_komega` / `_kepsilon` | 1 | RANS models with free surface |
| `cfd_2d_nwt_porous_komega` | 1 | VRANS porous box + k-ω porous sources (B295) |
| `cfd_2d_nug_breaking_komega` | 1 | stretched grid, cnoidal waves, k-ω T36=2, slope solid |
| `cfd_2d_dambreak` (+ `_fcc3`, `_fcls3`, `_mpi2`) | 1/2 | closed tank, walls, N40=33 / 4, 2D MPI |
| `cfd_2d_dambreak_teno` | 1 | TENO5 weights for momentum (D 12 2) and level set + reini (F 37 2) |
| `cfd_2d_dambreak_plic` (+ `_rk3`, `_fcc3`) | 1 | PLIC VOF (F80=4) with N40=3 / 13 / 33 |
| `cfd_2d_dambreak_picard`, `cfd_2d_dambreak_plic_picard_lsm` | 1 | Picard volume correction: F46=2 (once per step), F46=3 with PLIC redistancing (F92=1) |
| `cfd_2d_compressible_box` | 1 | sealed box, piston pushes water in: isentropic air compressibility (W 32 2), closed domain without pressure pinning |
| `cfd_2d_cylinder_singlephase` | 1 | single phase, explicit diffusion, inflow/outflow, forces |
| `cfd_3d_dambreak_obstacle` | 4 | 3D, MPI halos, solid box |
| `cfd_3d_dambreak_periodic` | 4 | periodic x and y with obstacle: REEFMG periodic halo (across ranks and on one rank) |
| `cfd_3d_nwt_narrow` | 1 | 3D wave flume 4 cells wide: REEFMG coarsening of a narrow direction |
| `cfd_3d_pier_komega` (+ `_rkls3_sf`, `_t33`) | 2 | 3D inflow/outflow, k-ω wall functions, cylinder, N40=14 sf loop, T33 k-gradient source |
| `cfd_3d_heave_sphere_6dof` (+ `_rk3`) | 4 | floating body 6DOF (FCLS3), N40=13→14 df loop |
| `cfd_3d_heave_sphere_6dof_mooring` | 2 | 6DOF with spring mooring (X310=4), all DOFs free, rotational damping |
| `cfd_3d_net_wall` | 2 | fixed quasi-static wall net (X 320 3) in an open-channel inflow, no floating body (sixdof_void): net forces spread with the Peskin kernel (net_interface::dlm_cfd) (references with patch 0022, Ship Hydrodynamics) |
| `nhflow_2d_6dof_box` (+ `_rk3`) | 1 | NHFLOW floating box in waves, two-way 6DOF, A510=2 / 3 |
| `nhflow_2d_6dof_box_init` | 1 | + linear damping X25/X26, initial pitch X101, initial velocities X102/X103 |
| `nhflow_2d_6dof_box_towed` | 1 | towed in surge (X11 u=2, X210) with velocity ramp X206 |
| `nhflow_2d_6dof_box_oneway` | 1 | prescribed motion, one-way coupling (X10=2) |
| `nhflow_2d_6dof_box_x15_5` | 1 | X 15 5: X 15 1 forcing with the sharp switch at the hull surface instead of the shell ramp (needs patch 0008, NHFLOW 6DOF code check) |
| `nhflow_2d_6dof_box_x15_6` | 1 | X 15 6 + A 520 1: X 15 1 forcing, faces between body and fluid cells reflective walls in the HLL fluxes (needs patch 0021) |
| `nhflow_2d_6dof_membrane_collar` | 2 | rigid membrane bag (X330) on a floating collar |
| `nhflow_3d_membrane_cylcone_current` | 2 | fixed membrane bag `cylcone` (cylinder on a cone, sloped floor) in a current: B60=1 + potential start, HLLC (needs patch 0002) |
| `nhflow_3d_membrane_cylcone_drain` | 2 | flexible `cylcone` bag, filling + drain, staggered coupling, compression; hydrograph inflow B60=2 (needs patches 0001, 0002) |
| `nhflow_3d_membrane_sharp_static80` | 2 | fixed `cylcone` bag filled to 80 % in still water, `mobility sharp` (nhflow_thinbody: blocked links, wall fluxes, head below the floor, rigid lid) (needs patches 0004-0006; references with 0007, 0008) |
| `nhflow_3d_membrane_sharp_current` | 2 | fixed `cylcone` bag in a current with `mobility sharp`: B60=1 + potential start, guarded fluxes near the bag (needs patches 0004-0006; references with 0007, 0008) |
| `nhflow_3d_membrane_sharp_kw` | 2 | as `nhflow_3d_membrane_sharp_current` with k-ω (A560=2): fabric as wall of the diffusion and the turbulence model (needs patches 0004-0006; references with 0007, 0008) |
| `nhflow_3d_membrane_sharp_mpi` | 2 | fixed 80 % `cylcone` bag centred on the subdomain border in a current ramp, `mobility sharp`: lid potential (CG), ghost faces and node sides across subdomains (needs patch 0007) |
| `nhflow_3d_membrane_sharp_flex` | 2 | flexible `cylcone` bag with `mobility sharp`: drain from 98 %, staggered coupling with added mass, bending hinges, smoothed geometry, normal wall velocity, hysteresis, moving floor (needs patches 0007, 0008) |
| `nhflow_3d_net_wall` | 2 | NHFLOW current (B60=1) through a fixed quasi-static wall net (X 320 3), no floating body (sixdof_void): net_interface::dlm_nhflow (references with patch 0022, Ship Hydrodynamics) |
| `nhflow_3d_6dof_box` | 2 | NHFLOW 3D box, all six DOFs free, initial roll/yaw |
| `nhflow_3d_6dof_box_towed_x15_3` | 2 | NHFLOW 3D box towed at 0.5 m/s in waves, heave and pitch free, X 15 3: forcing inside the body only, momentum fluxes at body faces from the fluid side (needs patch 0020) |
| `nhflow_3d_6dof_box_towed_x15_7` | 2 | as `nhflow_3d_6dof_box_towed_x15_3` with X 15 7 and HLLC (A 511 2): forcing inside the body only, body faces reflective walls in the HLL/HLLC fluxes (needs patch 0021) |
| `nhflow_3d_sphere_force` | 1 | NHFLOW force integration P 81 in still water: submerged sphere (A 586) on 6 sigma layers, marching-tetra quads split into two triangles: Fx = Fy = 0 to round-off, Fz = rho g V of the reconstructed surface (needs patch 0001 of NHFLOW Forces; reference with 0001) |
| `nhflow_3d_cylinder_force_mpi2` | 2 | P 81 on the vertical cylinder in the k-eps channel (base `nhflow_3d_cylinder_kepsilon_mpi2`), surface across the rank interface, globalsum of the forces (needs patch 0001 of NHFLOW Forces; reference with 0001) |
| `nhflow_3d_shipwave_box` | 1 | NHFLOW moving pressure patch, ship-wave mode (X10=3, X400=2) |
| `sflow_shipwave_box` (+ `sflow_6dof_box_oneway`) | 1 | SFLOW ship pressure patch (X10=3) / one-way direct forcing (X10=2) |
| `fnpf_2d_6dof_box` (+ `_rk4`) | 1 | FNPF resolved floating body, A310=3 / 4, added-mass coupling |
| `fnpf_2d_6dof_box_pto` | 1 | the FNPF box in heave with the 6DOF power take-off (X 500 2): PTO damper (X 502) and stiff end stops (X 504) as the CORAL PTO of `fnpf_2d_coral_fluid_body`; `REEF3D_6DOF_pto_0.dat` compared |
| `fnpf_2d_6dof_box_pto_loss` | 1 | `fnpf_2d_6dof_box_pto` with drive-train losses (X 510: loss damper and loss Coulomb friction, not absorbed power; P_loss and E_loss in `REEF3D_6DOF_pto_0.dat`); CORAL counterpart `fnpf_2d_coral_fluid_body_loss` (needs patch CORAL/0044) |
| `nhflow_3d_ship_box_thrust` | 2 | ship module (X350, ship.dat): constant thrust with offset, ITTC-1957 friction + form factor |
| `nhflow_3d_ship_box_6dof` | 2 | ship module, all six DOFs: cross-flow drag (strips), roll damping |
| `nhflow_3d_ship_box_propeller` | 2 | ship module propeller: KT/KQ, wake fraction, Hough-Ordway actuator disk as NHFLOW momentum source (all three components, swirl without net force), shaft torque reaction |
| `nhflow_3d_ship_box_zigzag` | 2 | ship module, surge/sway/heave/yaw: propeller with sampled inflow (NHFLOW velocity), MMG rudder, zig-zag steering with rudder rate |
| `nhflow_3d_ship_box_yawfree` | 2 | ship box with sway and yaw free, side faces through cell centres: level set = 0 on the surface, X 15 1 forcing ramp, no side force / yaw kick |
| `nhflow_3d_ship_kvlcc2_mmg` | 2 | ship module MMG model: KVLCC2 L7 hull derivatives, implicit added mass, MMG wake, rudder f_alpha / asymmetric gamma_R, fluid loads masked (mmg_fluid 0), zig-zag start |
| `nhflow_3d_ship_box_current` | 2 | ship held in a current: discharge inflow, hull turned by 180 deg (X 101), friction from the velocity relative to the current |
| `cfd_3d_ship_box_propeller` | 2 | ship module in CFD: box barge in a two-phase tank, actuator disk on the staggered velocity points (water part outside the hull, exact T and Q), velocity sampling, SSP-RK3 |
| `fnpf_3d_ship_box_hybrid` | 2 | ship module in FNPF: box in head waves, six DOFs, hybrid MMG (mmg_fluid 2, Munk moment out of N'v), approach phase from rest, propeller inflow sampled from the FNPF velocity, FNPF psi_0 = chi (body-following time derivative), X 17 relaxation of eta next to the footprint |
| `fnpf_3d_ship_box_hybrid_x18` | 2 | as `fnpf_3d_ship_box_hybrid` with X 18 0.5: the running mean of chi at body-fixed points (second-order low-pass) removed from the psi_0 free-surface data next to the hull (patch 0015) |
| `fnpf_3d_ship_box_hybrid_opts` | 2 | as `fnpf_3d_ship_box_hybrid` with the hybrid options of the ship module: `hybrid_munk` (A22 - A11 of the Munk moment taken out of N'v), `hybrid_surge 1` (MMG surge velocity term) and a reduced R0' (patch 0017) |
| `fnpf_3d_6dof_tumblehome_x27` | 1 | towed floating frustum with tumblehome (wider below the waterline) and X 27 2: columns with a body node within two sigma layers below their free-surface node belong to the footprint (patch 0019) |
| `cfd_2d_fem_obstacle` (+ `_mpi2`) | 1/2 | FEM solid (Z 30, N10=1): elastic obstacle hit by the bore, coupling across a subdomain border |
| `cfd_2d_fem_wall_failure` | 1 | FEM concrete wall cracking, erosion, the broken wall becomes a rigid fragment (`fragments rigid`), debris, ground and part contact (reference with 0024, adaptive integration) |
| `cfd_3d_fem_column` | 4 | FEM solid in 3D on 4 ranks: elastic column hit by the bore |
| `cfd_2d_cpm_bedload_layer` (+ `_susp`, `_mpi2`, `_k`, `_hyb`) | 1/2 | CPM MP-PIC sand bed in a periodic channel with the sub-grid bedload layer (Q 58 1 / 2): pickup, hops, deposition, release into the suspension and deposition of the suspended parcels at the bed, column bed level, one-sided solid forcing at the particle bed, pore pressure next to the bed, CPM log; `_k`: bed shear stress from the turbulent kinetic energy (Q 68 1); `_hyb`: hybrid suspension with the Eulerian concentration over the parcel bed (Q 58 3, patch 0017); since patch 0027 the bedload layer defaults to Q 72 2 and Q 68 1, so `_k` equals the base case (the k term is not active in this equilibrium channel; `_q72` covers Q 68 2) (references of all CPM bedload-layer cases with patch 0027) |
| `cfd_2d_cpm_bedload_layer_graded` | 1 | CPM bedload-layer channel on a grid with a change of cell size in x (2.5 mm on the first half, 5 mm on the second, `x-spacing.dat`): solid fraction from the kernel volume of the nodes, conservative smoothing across the change, bed shear stress with the vertical cell size, column bed level ending at the first almost empty cell (patch 0019; references of all 2D CPM cases with patch 0019; reference re-blessed with patch 0038: periodic ghost nodes from the opposite end) |
| `cfd_2d_cpm_bedload_layer_q72` | 1 | As `cfd_2d_cpm_bedload_layer` with the bed level seen by the fluid smoothed over the columns (Q 72 2, at most half a cell per column) and the bed shear stress with two velocity scales, sqrt(log law * 0.3 rho k) (Q 68 2) (patches 0024, 0026) |
| `cfd_2d_cpm_bedload_layer_q74` | 1 | As `cfd_2d_cpm_bedload_layer` with Q 74 1: grains of the bed that were not released into the suspension see no resolved flow and no turbulent dispersion, moved by the layer, gravity and the packed bed only (patch 0028) |
| `cfd_2d_cpm_bedload_layer_rk2` | 1 | As `cfd_2d_cpm_bedload_layer` with the RK2 time scheme of the parcels (Q 10 2): both stages with the bedload layer, the resting grains of the jammed bed sleep in both (patch 0032) |
| `cfd_2d_cpm_bedload_layer_mf` | 1 | As `cfd_2d_cpm_bedload_layer` with the morphological factor Q 65 10: pickup and hops 10 times faster, relaxation of the bed level seen by the fluid Q 63 2 s of bed time (patch 0023) |
| `cfd_2d_cpm_bedload_layer_merge` (+ `_mpi2`) | 1/2 | As `cfd_2d_cpm_bedload_layer_graded` with parcel weights, merge and split (Q 78 9): the resting grains of the jammed bed in the 5 mm cells merge in pairs to the target weight of their cell, awake and too heavy parcels split again; weights in the solid fraction, the grid-limited step, the bedload layer and the log; the CPM load map every 50 steps (G 45 50, `REEF3D_CFD_CPM_loadmap.dat`, a map of the cells); `_mpi2`: Q 78 4 on the uniform channel of `cfd_2d_cpm_bedload_layer_mpi2` (target weight 2), the weights across the packed parcel exchange of 2 ranks (patches 0034-0036) |
| `cfd_2d_cpm_wedge_yield` | 1 | CPM submerged 40 deg sand wedge in still water: jammed bed within the Mohr-Coulomb yield of the column (Q 67 1), packed-bed stress, Coulomb/mu(I) friction, deep slope failure (reference with patch 0037: bed level across cavities) |
| `cfd_2d_cpm_wedge_slide` | 1 | The sand wedge of `cfd_2d_cpm_wedge_yield` without the column yield (Q 67 0): sand slide of the bed surface beyond a slope of 0.62 on the smoothed column levels, top grains moved to the lower column, slide count in the CPM log (Q 75, patch 0029); fixed parcel bed for x < 0.05 m (Q 76, patch 0030; reference with patch 0037: bed level across cavities) |
| `cfd_2d_cpm_seepage_upflow` | 1 | CPM sand bed with an excess pore pressure of 700 Pa at the bottom (Q 69 1, Q 70, 1.44 x the fluidisation pressure): Laplace solve for the pore pressure, effective stress with the seepage force, unloading of the jammed bed, heave (patch 0014; reference re-blessed with patches 0034-0038, the reference of 0017 no longer matched hans_dev before them) |
| `cfd_2d_fem_simple_wall` | 1 | FEM simple input (concrete C30 preset, fix base, monitor auto, resolution), settling with the initial water, hybrid loads, structural damping |
| `cfd_2d_fem_floating` | 1 | FEM rigid body (floating box, `material rigid 500`) in the collapsing water column: rigid-body integration, pressure loads with the probe correction and the added-mass stabilisation |
| `cfd_2d_fem_debris_impact` | 1 | FEM rigid timber block (`stiffness 2e5`) on the dry bed pushed by the dam-break surge onto an elastic post: debris impact as one spring of the debris stiffness, ground contact, pressure beside the body for faces on the floor (reference with patch 0021: the floor of the domain is no second ground under `ground 0.0`) |
| `nhflow_2d_fem_wall` | 1 | FEM coupled with NHFLOW (Z 30, A 10 5): concrete wall holding a 0.3 m water column on the dry bed; cells inside the structure marked solid (p->DF -1, no continuity flux, d->solid_flux), water at rest, hydrostatic load (needs patches 0016, 0017, 0019; references with 0019) |
| `nhflow_2d_fem_floating` | 1 | FEM rigid floating box in NHFLOW still water: probed pressure (filtered non-hydrostatic part plus hydrostatic of the free surface), added-mass factor 5, heave towards the 5 cm draft (needs patches 0016, 0017, 0019; references with 0019) |
| `nhflow_2d_fem_debris_impact` (+ `_mpi2`) | 1/2 | FEM rigid timber block carried by an NHFLOW dam-break surge over the dry bed onto a concrete post that stands dry: wet/dry columns, debris contact, probes kept out of solid cells, subdomain border at the post (needs patches 0016–0021; references with 0024, adaptive integration) |
| `nhflow_2d_fem_light_debris` | 1 | FEM empty 20 ft container (rigid, 59.3 kg/m3, 2D slice) dropped 5 m onto 4 m of still NHFLOW water: added-mass factor limited to 30 body masses, speed relative to the water at most the terminal speed (needs patches 0016–0021) |
| `nhflow_2d_fem_light_debris_strong` | 1 | As `nhflow_2d_fem_light_debris` with strong coupling (`strong on`): the momentum step of a fluid step is repeated (fluid and FEM state saved and restored) until the motion of the container converges; no added-mass clamp or speed limit (needs patch 0025) |
| `cfd_2d_fem_rc_wall` | 1 | FEM reinforced concrete (`rebar`): weak concrete wall with 1 % model steel (f_y 30 MPa) holding the dam-break column; cracks at the base, the bars yield, the wall bends about 9 mm and stands; unilateral damage, steel utilisation and first yield in the summary (needs patch 0022; reference with 0024, adaptive integration) |
| `nhflow_2d_nwt_stokes5` (+ `_mpi2`) | 1/2 | NHFLOW relaxation generation + beach (B98=2/B99=1), Stokes 5th |
| `nhflow_2d_relax_ramp` | 1 | NHFLOW relaxation generation with the start-up ramp B 101 1: the ramp acts on the wave only, not on the still water depth of the zone (WL_relax, fe03b66a) |
| `nhflow_2d_nwt_stokes5_b95` | 1 | the same with the step-size consistent relaxation zones (B 95 0.01): targets at the stage output times, factor r^(dt/0.01) |
| `nhflow_2d_nwt_stokes5_p64` | 1 | the same with P 64 pressure probes in one column (near the bed, mid-depth, near the still water level, above it): the non-hydrostatic pressure P read from the vertical cell faces (ccipol7V) plus the hydrostatic part of the free surface; 0 above the surface (patch 261020) |
| `nhflow_2d_dirichlet` | 1 | NHFLOW Dirichlet wave generation (B98=3) |
| `nhflow_2d_awa` | 1 | NHFLOW active wave generation + active absorption (B98=4/B99=3) |
| `nhflow_2d_current` | 1 | NHFLOW waves on a current: discharge inflow (B60=1, W10), beach relaxes to the current (B97=1) |
| `nhflow_2d_irregular_decomp` | 1 | NHFLOW JONSWAP waves, decomposed relaxation precalc (B89=1) |
| `nhflow_3d_irregular_decomp` | 2 | NHFLOW 3D short-crested waves (B130=2), decomposed precalc, CFL 0.5 |
| `fnpf_2d_regular` (+ `_mpi2`) | 1/2 | FNPF relaxation generation + beach, regular waves |
| `fnpf_2d_break_kennedy` | 1 | FNPF Kennedy breaking (A 350 2) on steep deep-water Stokes waves with a low steepness threshold (A 355 0.25): 2D deep-water flag (cleared in the step it was set before fnpf breaking 01, run blew up at t = 13.6 s), two-pass viscosity spreading |
| `fnpf_2d_break_romero` | 1 | the same waves with the Romero model (A 350 3): stage values, depth from the bed, energy matching with the dissipation of the damping terms (fnpf breaking 02) |
| `fnpf_2d_break_barthelemy` | 1 | FNPF A 350 4 on the Ting & Kirby plunging breaker (nightly grid of the benchmark): crest tracking, onset B = u_s/C > 0.85 (first wave at x = 13.9 m), bore dissipation, per-crest viscosity, coastline viscosity, crest log P 59 |
| `fnpf_2d_coral_point_absorber` | 1 | REEF3D::CORAL inside REEF3D (Z 40 1): point absorber decay (Cummins, tether, PTO) in step with the FNPF run, `REEF3D_CORAL/results.dat` compared (needs patches CORAL/0001-0006) |
| `fnpf_2d_coral_fluid_body` | 1 | REEF3D::CORAL drives the FNPF box of `fnpf_2d_6dof_box` (Z 40 1, `fluid_body`): heave through a CORAL prismatic joint with PTO damper and stiff end stops, FNPF loads and added mass to CORAL, motion back at the RK stages with the predictor; `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0019) |
| `fnpf_2d_coral_fluid_body_loss` | 1 | `fnpf_2d_coral_fluid_body` with drive-train losses in the CORAL PTO (`c_loss`, `coulomb_loss`): they act on the box, are not absorbed power, channels `pto_Ploss`, `pto_Eloss` in `REEF3D_CORAL/box.dat`; 6DOF kernel counterpart `fnpf_2d_6dof_box_pto_loss` (needs patch CORAL/0044) |
| `fnpf_2d_coral_fluid_body_planar` | 1 | the same box in surge, heave and pitch through a CORAL planar joint, fluid loads extrapolated linearly without the predictor (`fluid_coupling order 1 no_predictor`) (needs patches CORAL/0011-0019) |
| `fnpf_2d_coral_kinematics` | 1 | REEF3D::CORAL one-way kinematics (Z 40 1, `fluid_kinematics`): Morison loads on a rigid and a flexible pile in the FNPF regular waves of `fnpf_2d_regular`, velocity and surface sampled at the strip centres once per fluid step; `REEF3D_CORAL/piles.dat` compared (needs patches CORAL/0011-0021) |
| `fnpf_2d_coral_one_way` | 1 | REEF3D::CORAL one-way body loads: the FNPF box of `fnpf_2d_6dof_box` held fixed (X 11 0), its loads on a rigid CORAL body (`fluid_body ... one_way`, CORAL one fluid step behind) and in strips along a clamped beam spanning it (`fluid_member`), journal of the fluid input (`record`); `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0021) |
| `nhflow_2d_coral_one_way` | 1 | the same with the NHFLOW box of `nhflow_2d_6dof_box` (loads at the end of each NHFLOW step and of the initialisation, `force_calc_stl` strips); `REEF3D_NHFLOW_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0023) |
| `nhflow_2d_coral_fluid_body` | 1 | REEF3D::CORAL drives the NHFLOW box of `nhflow_2d_6dof_box` (Z 40 1, `fluid_body`, two-way): heave through a CORAL prismatic joint with a PTO damper; the loads of the first stage at the predicted pose go to CORAL as provisional, CORAL steps, the later stages take its motion, the loads of the step's end are committed; CORAL's artificial added mass (`added_mass`); `REEF3D_NHFLOW_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0023) |
| `nhflow_2d_coral_fluid_body_x13` | 1 | `nhflow_2d_coral_fluid_body` with X 13 1 (the body advanced at the start of the stage): CORAL's poses set in `advance_stage`; `REEF3D_NHFLOW_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0025) |
| `nhflow_3d_coral_amr_sub` | 1 | `nhflow_3d_amr_heave_sub` on one rank driven by REEF3D::CORAL (Z 40 1, `fluid_body`, two-way): with the subcycling of the mesh refinement (G 7 1) every finest step is a CORAL step (predicted pose in the first stage, provisional loads, CORAL step, committed loads), the coarser level takes poses predicted from CORAL; `REEF3D_NHFLOW_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0025) |
| `cfd_2d_coral_fluid_body` | 1 | REEF3D::CORAL drives a CFD body (Z 40 1, `fluid_body`, two-way): heave decay of a 2D box (0.2 m x 0.15 m, density 500) in a level-set tank, low-storage RK3; the loads of the first stage committed, those of the second at the predicted pose provisional (8/15 dt), CORAL steps, the third stage on CORAL's motion; `REEF3D_CFD_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0026) |
| `cfd_2d_coral_one_way` | 1 | the box of `cfd_2d_coral_fluid_body` held (X 11 0): its CFD loads one-way on a rigid CORAL body (`fluid_body ... one_way`) and in 10 strips along a clamped beam spanning it (`fluid_member`, STL load integration); `REEF3D_CFD_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0026) |
| `fnpf_2d_coral_two_body` | 1 | REEF3D::CORAL two-body absorber in FNPF 2D (X 20 2, STL): a float released 2 cm above its equilibrium and a deeper second float ("spar") 1.2 m away in regular waves, the spar on a heave joint to the ground, the float on a prismatic joint to the spar with a PTO damper; the FNPF loads and the added mass with the cross terms of the two bodies go to CORAL; `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/two.dat` compared (needs patches CORAL/0011-0024) |
| `fnpf_2d_coral_drive` | 1 | forced heave of the FNPF box of `fnpf_2d_6dof_box` in still water through a CORAL joint drive (`drive ... sine`, 5 mm at T = 1 s with a ramp): the hydrodynamic load of a forced oscillation as for a database; `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/box.dat` compared (needs patches CORAL/0011-0024) |
| `fnpf_3d_coral_amr` | 1 | REEF3D::CORAL drives the FNPF point absorber of `fnpf_3d_amr_decay_sub` (Z 40 1, `fluid_body`, heave joint, no mooring) on a refined mesh (G 1 1, zone around the body): heave decay, the loads and the added mass of all grids to CORAL, CORAL's poses on all grids; `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/pa.dat` compared (needs patches CORAL/0011-0027) |
| `fnpf_3d_coral_amr_sub` | 1 | the same with subcycling (G 7 1): every finest step is a CORAL step, level 0 moves its copy of the body on poses CORAL predicts (needs patches CORAL/0011-0027) |
| `fnpf_3d_coral_amr_sub_mpi2` | 2 | `fnpf_3d_coral_amr_sub` on 2 ranks (the body and its refinement zone across the partition edge): every rank advances the same CORAL model on the loads of rank 0 (broadcast), the state of the model compared over the ranks after each step (the log warns at a difference) (needs patches CORAL/0011-0028) |
| `fnpf_3d_coral_array` | 1 | three FNPF point absorbers 1 m apart driven by CORAL (`X 20 3`, three `fluid_body`, a heave joint and a PTO damper each) with a refined zone around each body (G 1 1), released 3, 2 and 1 cm above equilibrium: the loads and the added mass with the cross terms of all bodies to CORAL, the poses back to each body (20 steps: CORAL = the 6DOF kernel with X 500 to 0.01 mm) (needs patches CORAL/0011-0030) |
| `fnpf_3d_coral_array_sub` | 1 | `fnpf_3d_coral_array` with subcycling (G 7 1): a CORAL step every finest step, the coarser levels on the predicted poses of all bodies (needs patches CORAL/0011-0030) |
| `fnpf_3d_coral_amr_fullscale` | 1 | `fnpf_3d_coral_amr` at full scale (Froude x100: dx 8 m, depth 100 m, the point absorber 30 m wide): the composite psi solves (pressure, added mass) with a resolved body on the patches; before patch CORAL/0047 the FAC preconditioner restricted a stale increment at the body nodes of the patches with the residual (small against it at model scale, dominant at full scale) and the psi solves diverged to NaN at the first step (CORAL: time step failed at t = 0); A33 1e6 times that of `fnpf_3d_coral_amr` within 0.6 % (needs patch CORAL/0047) |
| `fnpf_3d_coral_amr_fullscale_sub` | 1 | `fnpf_3d_coral_amr_fullscale` with subcycling (G 7 1): the window preconditioner (lap_prec_win) at full scale (needs patch CORAL/0047) |
| `fnpf_3d_coral_array_am3` | 1 | `fnpf_3d_coral_array` with the added mass every 3 time steps (X 28 3): the unit-mode solves of the three bodies and their cross terms at every third step, held in between (needs patch CORAL/0049) |
| `fnpf_3d_coral_amr_fullscale_am4` | 1 | `fnpf_3d_coral_amr_fullscale` with the added mass every 4 steps (X 28 4), earlier after a motion by more than 0.5 m (X 29 0.5 0), a checkpoint every 7 steps (P 48 7) (needs patch CORAL/0049) |
| `fnpf_3d_coral_amr_fullscale_am4_restart` | 1 | the restart (I 45 1) of `fnpf_3d_coral_amr_fullscale_am4` from step 14: the added mass of the checkpoint (fnpf_6DOF::ckpt_state version 2) until the next one is due; the final state bitwise that of the full run. The first restart case with mesh refinement under the AMR contract check (needs patches CORAL/0048, 0049) |
| `fnpf_3d_coral_amr_ckpt` | 1 | `fnpf_3d_coral_amr` with a checkpoint every 5 steps (P 48 5), the newest three kept (P 49 3: steps 10, 15, 20, the last also the checkpoint at the end of the run); the uninterrupted run of the restart cases (needs patches Checkpoints/0001-0008) |
| `fnpf_3d_coral_amr_restart`, `_restart_mpi2` | 1, 2 | the restart (`chain`) of `_ckpt` from its checkpoint at step 10 (I 45 10, not the latest) with the outputs of the full run copied in (CORAL, 6DOF, wave gauges): they are cut back to step 10 and continued; `ranks`: on 1 rank the final state and every text output bitwise those of the uninterrupted run, on 2 ranks the text outputs within 1e-6 (needs patches Checkpoints/0001-0008) |
| `nhflow_2d_coral_fluid_body_mpi2` | 2 | `nhflow_2d_coral_fluid_body` on 2 ranks (needs patches CORAL/0011-0028) |
| `cfd_2d_coral_one_way_mpi2` | 2 | `cfd_2d_coral_one_way` on 2 ranks (the held box across the partition edge): the loads and the 10 strips summed over the ranks (needs patches CORAL/0011-0028) |
| `fnpf_2d_coral_drag` | 1 | the FNPF box of `fnpf_2d_coral_fluid_body` with CORAL's quadratic drag (`drag ... Cd area`) in place of the PTO: heave decay; `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/box.dat` (with the drag channels) compared (needs patches CORAL/0011-0029) |
| `fnpf_2d_coral_mooring` | 1 | the FNPF box of validation 15 (Z 40 1, `fluid_body`, planar joint) in regular waves, moored fore by a quasi-static line (elastic catenary with seabed) and aft by a dynamic line (lumped mass, 10 segments, `fluid_kinematics`); `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/box.dat` (line forces, `Ta`/`Tb`) compared; re-blessed with CORAL/0031-0035 (the dynamic line starts from the static equilibrium of its segments: before, it rattled between slack and taut) (needs patches CORAL/0011-0035) |
| `fnpf_2d_coral_mooring_system` | 1 | the moored box of `fnpf_2d_coral_mooring` with a mooring system: fore a chain-rope-chain line (`line_type`, `section`, joined at points; the rope with `BA_slack`), aft two chains through a clump weight (`point`, resting on the bed; chain and clump with stick-slip friction), and a cable with bending stiffness (`EI`, `clamp_a`) that lets go of the ground at 1.5 s (`break_time`); anchors 2 cm above the FNPF bed, on which lines and clump rest (`bathymetry fluid`); `REEF3D_FNPF_6DOF/*.dat` and `REEF3D_CORAL/box.dat` (line forces, point positions, `umb_broken`) compared (needs patches CORAL/0011-0043) |
| `fnpf_2d_coral_mooring_system_lagoon` | 1 | `fnpf_2d_coral_mooring_system` with the LAGOON store instead of the VTP files (P 18 1, bodies every 0.1 s): the body and CORAL's six mooring lines (particles/coral_mooring: lines, tension at the points, line number) at the same outputs; the solution that of the case without the store (needs patch CORAL/0050) |
| `fnpf_2d_irregular_decomp` | 1 | FNPF JONSWAP waves, decomposed relaxation precalc (B89=1) |
| `fnpf_3d_shortcrested` | 2 | FNPF 3D short-crested waves (B130=2) |
| `nhflow_2d_two_sources` (+ `_mpi2`) | 1/2 | wave_field: B 92 linear wave + linear source 2 with phase (B 500/501); P 50 prints the summed target |
| `nhflow_2d_irregular_two_sources` | 1 | wave_field: two JONSWAP sources with their own seeds (B 139, B 504) |
| `fnpf_2d_two_sources` | 1 | wave_field in FNPF (cached potential path) |
| `cfd_2d_dambreak_amr` (+ `_mpi2`) | 1/2 | CFD AMR (G 1 1, static box G 15): dam break with one refined level around the water column, the front leaves the patch (fills, level set per grid, composite projection with matched fluxes at the patch edges, FAC-preconditioned BiCGStab, restriction); `_mpi2`: the box ends at the partition edge (G 4 4), the refined region is extended across it, patches on both ranks, remote fills |
| `cfd_2d_dambreak_amr3` | 1 | CFD AMR with two refined levels (G 1 2): fills from a patch parent, FAC V-cycle over three levels |
| `cfd_2d_nwt_amr` | 1 | CFD AMR in the wave tank: refined band at the free surface between the relaxation zones (level 0 only), the free surface across the patch edges |
| `cfd_2d_nwt_amr_zone` | 1 | the refined band over the numerical beach (B 99 2, G 7 0): the patch relaxes towards still water as level 0; the tiles at the outflow edge stay on level 0 |
| `cfd_2d_nwt_amr_zone_sub` | 1 | the band over the wave generation zone (B 98 2) with subcycling: the patch relaxes towards the wave theory of level 0 at its cells and the time of its step |
| `cfd_2d_nwt_amr_zone_komega_sub` | 1 | k-omega in the wave tank, the band over the generation zone (G 7 1): the patch relaxes the eddy viscosity in the zone as level 0 (iowave::turb_relax) |
| `cfd_2d_nwt_amr_active_sub` | 1 | active wave generation and absorption at the boundaries (B 98 4, B 99 3), two levels over the whole tank (G 7 1): the tiles at the in- and outflow boundary, and those whose nesting would reach them, stay on the coarser level |
| `cfd_2d_nwt_amr_porous` | 1 | k-omega in the wave tank with the porous box (B 270, x 5-6 m) inside a refined box (G 7 0): vrans_f on the patch (the zones at its cells), resistance, porosity factors and turbulence sources there |
| `cfd_2d_nwt_amr_porous_sub` | 1 | the same with the band z 0.3-0.75 refined and subcycling: the patch edge across the porous box (CPOR n at the faces of the patch edges, in the predictors and the level projections) |
| `cfd_2d_still_amr` | 1 | CFD AMR still water, the free surface across a patch edge: the hydrostatic pressure of both grids balanced at the edge |
| `cfd_3d_dambreak_amr` | 2 | CFD AMR 3D dam break: 3D boxes on 2 ranks, implicit diffusion of a patch with Dirichlet values around it |
| `cfd_2d_dambreak_amr_sub` | 1 | CFD AMR with subcycling (G 7 1): the patch takes two steps per level-0 step, the cells around it linear in time; the level projection with the coarse pressure at the patch edges, the synchronisation projection over all levels at the end of the level-0 step |
| `cfd_2d_dambreak_amr3_sub` | 1 | CFD AMR, two levels with subcycling: a patch as the coarser level in time (its old and new state), level 2 four steps per level-0 step |
| `cfd_2d_dambreak_amr_mpi2_sub` | 2 | CFD AMR subcycling with patches on both ranks next to each other: one component of the level projection across the ranks (block-Jacobi patch V-cycles) |
| `cfd_2d_still_amr_sub` | 1 | CFD AMR still water across a patch edge with G 7 1 and the volume correction over all levels (F 46 3) |
| `cfd_2d_still_amr_porous_sub` | 1 | still water with a porous box across the patch edge and the free surface (G 7 1): the hydrostatic balance of the porous faces at the patch edge |
| `cfd_2d_dambreak_amr_vol` | 1 | CFD AMR (G 7 0) with the volume correction over the leaf cells of all grids (F 46 2) |
| `cfd_2d_dambreak_amr_weno3_sub` | 1 | CFD AMR with WENO3 (D 10 7, F 35 8) and G 7 1: the patch builds its own WENO3 coefficient tables (weno3_nug_func) |
| `cfd_2d_dambreak_amr_hires` | 1 | CFD AMR (G 7 0) with a TVD scheme on all grids (hires, D 10 13, F 35 13) |
| `cfd_2d_dambreak_amr_ref_sub` | 1 | CFD AMR with G 7 1 and the pressure fixed at a point in the air (B 30 1, B 33 2): the reference over all grids once at the end of the step |
| `cfd_3d_dambreak_amr_sub` | 2 | CFD AMR 3D dam break with subcycling on 2 ranks |
| `cfd_3d_dambreak_amr_cn_sub` | 2 | the same with the Crank-Nicolson diffusion (D 20 3) on all grids |
| `cfd_2d_dambreak_amr3_adapt_edges` | 1 | two levels adaptive (G 1 2, G 24 3, regrid every 4 steps): level-2 patch edges on the edges between level-1 patches (the fine flux into the patch that owns the coarse face and into the copy of the coarse cell) |
| `cfd_2d_dambreak_amr_adapt` | 1 | CFD AMR adaptive (G 24 3: cells with \|phi\| < 3 psi flagged), regrid every 2 steps without hysteresis (G 2 2, G 5 0): flags, rebuild, transfer from the old patches, prolongation into new cells, freeing of patches |
| `cfd_2d_dambreak_amr_adapt_sub` | 1 | the same with subcycling (G 7 1): the regrid after the level-0 step |
| `cfd_2d_dambreak_amr3_adapt` | 1 | two levels adaptive (G 1 2, G 24 3): new level-1 patches as the parents of new level-2 patches |
| `cfd_2d_dambreak_amr3_adapt_sub` | 1 | two levels adaptive with subcycling (G 7 1): touching patches of a level share one multigrid in the level projections and the synchronisation; the pressure of level 1 held at its end value for level 2 |
| `cfd_2d_dambreak_amr_adapt_mpi2` | 2 | adaptive on 2 ranks: patches moving across the partition edge |
| `cfd_2d_dambreak_amr_flags` | 1 | the criteria together: interface band G 24 1, vorticity G 26 20, a static box G 15 |
| `cfd_3d_dambreak_amr_adapt` | 2 | CFD AMR 3D adaptive (G 24 1, tiles G 4 2, no buffer G 3 0) on 2 ranks |
| `cfd_2d_komega_amr` | 1 | CFD AMR with k-omega (T 10 2, T 45 1): open channel, box over the depth: kin/eps/eddyv fills, the model on the patch, restriction of kin, eps, eddyv |
| `cfd_2d_komega_amr_sub` | 1 | the same with G 7 1 (kin, eps, eddyv around the patch linear in time) |
| `cfd_2d_dambreak_amr_komega_adapt` | 1 | dam break with k-omega, adaptive (G 24 3), regrid every 2 steps: the turbulence fields and old time levels of new patches |
| `cfd_2d_dambreak_amr_solid` | 1 | cylinder (S 32) in the dam-break tank, adaptive G 24 3 + G 25 3 (cells near the solid): solid level set prolonged, direct forcing on the patches |
| `cfd_2d_dambreak_amr_solid_sub` | 1 | the same with G 7 1 |
| `cfd_2d_dambreak_amr_heat` | 1 | heat on the CFD patches (H 10 3, G 7 0): the refinement transports T on all grids in flux form (SSP-RK3, first-order upwind faces H 15 1, explicit diffusion D 20 1), flux registers at the patch edges refluxed per step; the density from T (H 3 1) at the patch edges and the cells around the patch |
| `cfd_2d_dambreak_amr_heat_sub` | 1 | the same with subcycling (G 7 1), WENO face values (H 15 4) and implicit diffusion (D 20 2): refluxing after the two steps of the patch |
| `cfd_2d_dambreak_amr_conc_sub` | 1 | concentration on the CFD patches (C 10 3, C 15 4, G 7 1): the density from C at the patch edges (C 1, C 3), C zero in the air as concentration_RK3 |
| `cfd_2d_dambreak_amr_heat_adapt_mpi2` | 2 | heat with adaptive refinement (G 24 3) and G 7 1 on 2 ranks: T through the regrids (old patches, prolongation) and the patch edges across the rank boundary |
| `cfd_2d_heave_amr` | 1 | floating box (X 10 1) heave decay, static box around the body and the surface: the hull cast on the patch grid in every stage, its forcing on the patches |
| `cfd_2d_heave_amr_sub` | 1 | the same with G 7 1: the finest level advances the body with the loads on its grid, level 0 steps with a predicted copy |
| `cfd_2d_heave_amr_ref_sub` | 1 | the same with the gage at the free surface above x = 0.3 (B 30 2, B 33 1): the gage of the loads from the finest grid |
| `cfd_2d_heave_amr2_sub` | 1 | the same with two refined levels: level 2 advances the body, level 1 and level 0 step with predicted copies |
| `cfd_2d_heave_amr_adapt` | 1 | the cells within 3 cells of the body (G 25 3), regrid every 2 steps: the hull cast on new patches |
| `cfd_2d_heave_amr_adapt_sub` | 1 | the same with G 7 1: the patches that follow the body advance it |
| `cfd_3d_heave_amr_mpi2_sub` | 2 | 3D floating box with G 7 1 on 2 ranks, the rank cut through the body and the refined box: the loads on the patches up to the cut; a second box on one rank (different patch counts per rank: no MPI exchange in the body forcing of the patches) |
| `cfd_2d_dambreak_amr_place_mpi2` | 2 | CFD AMR with the placement of the patches (G 40 1): the refined box (over rank 0's level 0) cut vertically into two units filling the ranks to the mean work; the unit on rank 1 reads its coarser level from its shadow (the level-0 cells under it and around it), the restriction, the coarse faces and the fluxes of the patch edges go back to rank 0 (patch 0058-0060) |
| `cfd_2d_dambreak_amr_adapt_place_mpi3` | 3 | adaptive refinement with G 40 1 on 3 ranks, 48 steps: units fill the ranks in turn, a shadow over the level 0 of two ranks, the state of new patches from old patches on other ranks (regrid_old) |
| `cfd_2d_dambreak_amr3_adapt_sub_place2` | 2 | two adaptive levels with G 7 1 and the test of the placement (G 40 2), 20 steps: small units placed anew at every regrid on the least loaded rank, the coarser level and the old patches mostly on the other rank; the overlap of the patch multigrids (RAS) at the cuts |
| `cfd_2d_dambreak_amr_heat_adapt_place_mpi2` | 2 | heat, adaptive, G 7 1, G 40 1: the coarse fluxes of the flux registers from the shadows, the refluxing added on the rank of the coarse cell |
| `cfd_2d_nwt_amr_porous_place_mpi2` | 2 | waves through a porous box in a refined box (k-omega), G 40 1: the porosity of the shadows, the turbulence fields in the fills |
| `cfd_3d_dambreak_amr_place` | 2 | 3D dam break, G 40 1: the refined box cut into two units (load max/mean 1.53 without G 40, 1.16 with it: the cut at a tile edge) |
| `cfd_3d_heave_amr_place_mpi2_sub` | 2 | 3D floating box with G 7 1 and G 40 1: the hull triangles of the loads on the rank of the finest grid (amr_cfd_mine), the body forcing of the patches through the shadows |
| `cfd_3d_spin_amr_sub` | 1 | 3D paddle (X 110) small against its refined box, G 1 2 and G 7 1, rotation about x prescribed with a ramp: three steps at rest (no body cast while the triangles are unchanged), then spinning; the body cast and reinitialised in the box of the body on both levels (needs ship_propeller_cfd_amr 261018 patches 1, 2) |
| `cfd_3d_spin_amr_sub_g42` | 1 | `cfd_3d_spin_amr_sub` with G 42 0.1: a patch casts the paddle again only after a vertex has moved more than 0.1 cells since its last cast (needs ship_propeller_cfd_amr 261018 patch 3) |
| `cfd_2d_dambreak_amr_periodic` | 1 | adaptive CFD refinement in a tank periodic in x (C 21 1): the water column 0.05-0.20 collapses to both sides, the left front crosses the periodic edge; the refined region kept alike on both sides of the edge, the margin across it from the patches at the other end, the level-0 multigrid of the composite solve periodic (patch 0061) |
| `cfd_2d_dambreak_amr_periodic_mpi2` | 2 | the same with the periodic edge between the ranks: remote siblings across the edge, the overlap of the multigrids at the sides to patches of the other rank |
| `cfd_2d_dambreak_amr_periodic_band` | 1 | a static band over the whole periodic length: one patch periodic over itself (its margin across the edge from itself, its multigrid periodic) |
| `cfd_2d_dambreak_amr3_periodic` | 1 | two adaptive levels, periodic: the nesting wrapped across the edge, the parents of level-2 margin cells in the level-1 patch at the other end |
| `cfd_2d_dambreak_amr_periodic_sub` | 1 | `cfd_2d_dambreak_amr_periodic` with G 7 1: the level projections and the synchronisation with the joint multigrid across the edge |
| `cfd_2d_dambreak_amr_heat_periodic_mpi2` | 2 | heat, adaptive, G 7 1, periodic on 2 ranks: the scalar fills and the refluxing with patches across the edge, the joint multigrids with overlap cells |
| `cfd_3d_dambreak_amr_periodic` | 2 | 3D, periodic in x (across the ranks) and y (on each rank): the box at x = 0 gets its counterpart at the other end of x, the patches periodic over themselves in y |
| `cfd_2d_patch_pressure_amr` | 1 | `cfd_2d_patch_pressure` with a refined band over the whole channel: the tiles at the inlet and outlet patches (B 441) stay on level 0, the patches keep off the faces of B 440-442 (patch 0062) |
| `cfd_2d_dambreak_amr3_nested` | 1 | a static box with three nested levels (G 1 3, G 15), each level one patch inside the coarser one: the multiplicative FAC V-cycle through all of them (patch 0068) |
| `cfd_2d_dambreak_amr4_adapt` | 1 | four adaptive levels (G 1 4, G 24 3) in the dam break, 6 steps: the FAC V-cycle over five grids, every correction passed on into the finer levels, with the interface band and the covered couplings (patch 0068) |
| `cfd_2d_dambreak_amr4_adapt_sub` | 1 | `cfd_2d_dambreak_amr4_adapt` with G 7 1, 3 steps: the level projections of a level with patches in no joint multigrid, their overlap on one rank (patch 0069); the synchronisation with four levels |
| `cfd_2d_dambreak_amr_adapt_place_w2` | 2 | `cfd_2d_dambreak_amr_adapt` on 2 ranks with G 40 1 and G 41 2.0: a patch cell counted as two level-0 cells in the placement (patch 0070) |
| `cfd_2d_dambreak_amr_periodic_place_mpi2` | 2 | `cfd_2d_dambreak_amr_periodic_mpi2` with G 40 1: the units placed for the load with the periodic edge between the ranks, the shadows take the level-0 cells across the periodic edge from the other end (patch 0076) |
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
| `nhflow_3d_amr_waves3` (+ `_sub`) | 2 | waves with three levels (G 1 2), one step for all levels; `_sub` subcycled (G 7 1): a level between two levels, two synchronisations per level-0 step |
| `nhflow_3d_amr_heave3_sub` | 2 | heave box with three levels and G 7 1: level 1 steps with a predicted copy of the body, level 2 advances it |
| `nhflow_2d_amr_waves` (+ `_lev2`, `_mpi2`, `_sub`) | 1/1/2/1 | NHFLOW AMR on a 2D grid (one cell in y): regular waves in a 24 m flume entering a static patch over 5-10 m; x-only refinement (2x1 children), the patch and level-0 pressure solves exact (reefamr_xzlu, patch 0083); `_lev2` two levels, `_mpi2` 2 ranks (level-0 V-cycle), `_sub` G 7 1 |
| `sflow_2d_amr_dambreak` (+ `_mpi2`, `_place`, `_place2`) | 1/2 | SFLOW AMR dam break (tutorial 8_6): 2 levels following the surface jump (G 20), regrid every 4 steps; `_mpi2` patches cut at the rank boxes, `_place` G 40 1, `_place2` G 40 2 (remote parents, migration: restrict_levels, ini_patch_old/prolong) |
| `sflow_2d_amr_bar_nh` (+ `_place2`) | 1/2 | SFLOW AMR non-hydrostatic bar (tutorial 8_7, A 220 1): static boxes, composite pressure; `_place2` through the block plans |
| `sflow_2d_amr_bar_bous` (+ `_place2`) | 1/2 | the bar with Boussinesq (A 220 4): u_a on the leaf cells of all levels; `_place2` G 40 2 |
| `sflow_2d_amr_dambreak_sub` (+ `_place2`) | 1/2 | the dam break with subcycling (G 7 1): level l takes 2^l steps per level-0 step, fills interpolated in time, flux registers and refluxing, regrid every level-0 step; `_place2` G 40 2 (remote fills and refluxing) |
| `sflow_2d_amr_bar_sub` | 1 | hydrostatic waves over the bar (A 220 0), static boxes, G 7 1 |
| `sflow_2d_amr_bar_sub_b95` | 1 | the subcycled bar with B 95 -1 (reference step of the finest level) |
| `sflow_2d_amr_bar_nh_sub` (+ `_place2`) | 1/2 | the non-hydrostatic bar (A 220 1) subcycled with B 95 -1: pressure level by level (cells under the finer level fixed), synchronisation projection; `_place2` G 40 2 |
| `sflow_2d_amr_bar_bous_sub` | 1 | the Boussinesq bar (A 220 4) subcycled: u_a level by level, u_a rebuilt after refluxing |
| `sflow_2d_amr_ship` | 1 | `sflow_shipwave_box` (X 10 3, A 220 2) with 2 levels around the towed box (G 12 margin, G 13 wake wedge), regrid every 4 steps |
| `sflow_1d_amr_bar_nh` (+ `_mpi2`) | 1/2 | `sflow_2d_amr_bar_nh` on a grid with one cell in y, the patch over 5-13 m that the waves enter: x-only refinement, the patch line solves exact (reefamr_xzlu), on one rank the level-0 line too (reefmg2D::set_exact_line, patch 0084); `_mpi2` 2 ranks |
| `sflow_1d_amr_bar_nh_sub`, `sflow_1d_amr_bar_sub` | 1 | the non-hydrostatic (B 95 -1) and the hydrostatic bar with one cell in y subcycled (G 7 1): level solves, synchronisation projection, refluxing with one fine face per coarse face |
| `sflow_1d_amr_bar_bous` | 1 | the Boussinesq bar (A 220 4) with one cell in y: u_a, x-only restriction and prolongation of the line systems |
| `sflow_2d_amr_ship_sub` (+ `_df_sub`) | 1 | the same subcycled (G 7 1): the body on the patches linear in time within the level-0 step; `_df_sub` direct forcing (X 10 2) |
| `sflow_2d_curv_dambreak_rot90` | 1 | SFLOW on a curvilinear grid (G 50 2): dam break on the Cartesian grid rotated by 90 deg and shifted (F 72 in physical coordinates); rotated HLL fluxes, walls, CFL |
| `sflow_2d_curv_lake_island` | 1 | lake at rest with an emerged island on a wavy, skewed grid: well balancing, bed levels from the nodes |
| `sflow_2d_curv_channel_rot30` | 1 | river reach rotated by 30 deg at UTM coordinates: discharge inflow along the face normal, stage outflow, emerging island |
| `sflow_2d_curv_bend_mpi2` | 2 | 90 deg bend with Manning friction on two ranks |
| `sflow_2d_curv_bend_lagoon_mpi2` | 2 | the bend with its free surface and bed in the LAGOON store as well (P 18 2): `xc`, `yc` of the curvilinear grid (needs patch 0020) |
| `sflow_2d_curv_river_ini_bend_mpi2` | 2 | sloping river (slope 0.002) through a 90 deg bend, free surface and velocities from the steady 1D backwater solution of W 10 (F 65 2) on two ranks (needs patch 0001 of Curvilinear/patches) |
| `sflow_2d_curv_river_ini_linear` | 1 | straight sloping river reach, free surface linear from F 61 to F 62 along the river (F 65 1) (needs patch 0001 of Curvilinear/patches) |
| `sflow_2d_curv_river_inflow_conveyance_mpi2` | 2 | the sloping bend with a rising hydrograph (10 to 120 m3/s in 10 s) onto the floodplains, discharge across the inflow face by conveyance (B 62 1) on two ranks (needs patch 0003 of Curvilinear/patches) |
| `sflow_2d_curv_wavy_nh` | 1 | non-hydrostatic pressure (A 220 2) on the wavy, skewed grid: chain-rule gradients, compact Poisson operator, deferred non-orthogonal correction |
| `sflow_2d_curv_wavy_bous` | 1 | Boussinesq (A 220 4) on the wavy, skewed grid: inversion per grid-aligned component, matrix-free BiCGStab with line preconditioner |
| `sflow_2d_curv_amr_bous_sub` | 1 | AMR on the rotated grid: Boussinesq undular bore, two levels following the surface, regrid, subcycling (G 7 1): composite curvilinear inversion, fills of the grid-aligned components in time |
| `sflow_2d_curv_amr_nh_place2` | 2 | AMR on the wavy grid: non-hydrostatic pressure, regrid, patches placed on two ranks (G 40 1): curvilinear composite pressure, patch nodes from the global curvilinear nodes |
| `sflow_2d_mb_quadrants` (+ `_mpi2`) | 1/2 | SFLOW multi-block grid (G 50 1, DIVEMesh R 20 1, `blocks.dat`): four curvilinear blocks around a regular vertex, indices turned by 0/90/180/270 deg, dam break over the vertex: blocks as level-0 patches, interface transforms, ring of the native grid, remote fills |
| `sflow_2d_mb_amr` (+ `_mpi2`, `_place`) | 1/2 | the quadrants with two levels in a static box over the vertex: patches cut at the block faces, fills and turned flux matching across the interfaces, rank boxes / placement (G 40 1) |
| `sflow_2d_mb_amr_adaptive` | 1 | the quadrants with two levels following the surface jumps, regrid every 4 steps: tags of the blocks, refinement continued across the interfaces |
| `sflow_2d_mb_wetdry_amr` | 2 | two blocks (the second turned by 90 deg), dam break running up a slope: wetting and drying across the interface, shoreline refinement (G 22) |
| `sflow_2d_mb_singular` | 1 | a hexagon split into three blocks (singular vertex): degenerate ring corner cells solid, no refinement around the vertex |
| `nhflow_3d_curv_dambreak_rot90` | 1 | NHFLOW on a curvilinear grid (G 50 2): dam break on the Cartesian grid turned by 90 deg, A 520 2: normal velocities, HLL fluxes, walls, sigma metrics, curvilinear Poisson rows |
| `nhflow_3d_curv_dambreak_rot30_mpi2` | 2 | the basin turned by 30 deg at UTM coordinates on two ranks |
| `nhflow_3d_curv_dambreak_rot30_lagoon_mpi2` | 2 | the basin's volume, free surface and bed in the LAGOON store instead of the VTU/VTP files (P 18 1), with `xc`, `yc` of the curvilinear grid (needs patch 0020) |
| `nhflow_3d_curv_wavy_lake` | 1 | lake at rest over a bump on the wavy, skewed grid (A 520 2): well balancing, non-orthogonal Poisson rows |
| `nhflow_3d_curv_waves_rot30` | 1 | wave channel turned by 30 deg: Dirichlet generation (B 98 3) and active absorption (B 99 3) along the face normals, wave frame B 105 |
| `nhflow_3d_curv_waves_relax_map30` | 1 | the channel stretched along x (orthogonal mapped grid) and turned: relaxation zones by B 107 / B 108 in physical coordinates (P3 gate) |
| `nhflow_3d_curv_river_kepsilon_rot30` | 1 | river reach turned by 30 deg: discharge inflow along the face normal, stage outflow, k-epsilon with the implicit diffusion (A 512 2) |
| `nhflow_3d_curv_river_ini_bend_mpi2` | 2 | NHFLOW, sloping river through a 90 deg bend, free surface and velocities on all layers from the steady 1D solution with the log law of the bottom layer (F 65 2), two ranks (needs patch 0004 of Curvilinear/patches) |
| `nhflow_3d_curv_river_inflow_conveyance` | 1 | NHFLOW, the bend with a rising hydrograph (10 to 120 m3/s in 10 s), log inflow profile across the face by conveyance (B 62 1) (needs patch 0005) |
| `nhflow_3d_curv_river_inflow_conveyance_plain` | 1 | as above with the plain inflow profile (B 61 1): u = Uk h^(2/3) (needs patch 0005) |
| `nhflow_3d_curv_tide_rot30` | 1 | oblique tide through a basin turned by 30 deg: Riemann and Flather edges on the edge-face normals |
| `nhflow_3d_curv_amr_adaptive_sub` | 1 | NHFLOW AMR on the turned grid: patches following the surface, regrid every step, subcycling (G 7 1) |
| `nhflow_3d_curv_amr_wavy_lake` | 1 | AMR on the wavy grid: lake at rest under six patches (restriction per area, bed prolongation with the parent's area mean) |
| `nhflow_3d_curv_amr_shore_rot90` | 1 | AMR with wetting and drying in the patches (G 30 1) on the turned grid: run-up on a plane beach |
| `nhflow_3d_mb_quadrants` (+ `_mpi2`) | 1/2 | NHFLOW multi-block grid (G 50 1, DIVEMesh R 20 1, `blocks.dat`): four curvilinear blocks around a regular vertex, indices turned by 0/90/180/270 deg, wavy skewed grid, non-hydrostatic dam break over the vertex: blocks as level-0 patches, stage fills of the blocks and the ring, composite pressure over the blocks (block-local V-cycles, coarse space G 53), isotropic walls (G 52) |
| `nhflow_3d_mb_amr` | 1 | the quadrants with one level in a static box over the vertex: patches cut at the block faces, fills and turned flux matching across the interfaces |
| `nhflow_3d_mb_idiff` | 1 | the quadrants with the implicit diffusion (A 512 2): one composite diffusion solve over the blocks |
| `nhflow_3d_mb_potential` | 1 | channel with discharge inflow and stage outflow on block 0, blocked by a dry obstacle, bypass through three turned blocks: potential flow initialisation (I 11) as one composite solve over the blocks |
| `nhflow_3d_mb_open_faces_series_mpi2` | 2 | `nhflow_3d_mb_open_faces_mpi2` with the open faces from time series files: G 54 1 3 1 (discharge-1.dat through the -i face of block 1) and G 54 2 4 2 (waterlevel-2.dat at the +i face of block 2), linear in time (ioflow_amr_block); constant files give the base case bitwise (needs patch 0002, Tagus_NHFLOW_multiblock) |
| `nhflow_3d_mb_edges_mpi2` | 2 | edge zones on the faces of the blocks (B 520-523 with B 522 b): a Riemann edge with a harmonic tide on the -i face of block 1 and a Flather edge on the +i face of block 2, non-hydrostatic, two ranks (ioflow_amr_block::edges_nhflow, open_xm / open_xp of the block); one block with the iowave edges gives the same solution (Dropbox Tagus_NHFLOW_multiblock/validation/01_block_edges) (needs patches 0002, 0003) |
| `nhflow_3d_si_mb_edges` | 1 | the same three blocks with the semi-implicit NHFLOW (A 534 1): a clamped discharge edge (B 525) on block 1 and a clamped level edge with a time series background (background-2.dat) on block 2; nhflow_momentum_SI::bc_method takes the zones of its own block (lexer::amr_block) (needs patches 0002, 0003) |
| `nhflow_3d_curv_sed_dambreak_rot90` | 1 | sediment transport (S 10 1, bed load, S 31 1) with NHFLOW on the curvilinear grid turned by 90 deg: bed-load dam break, per cell index the Cartesian run (sediment_curv: face normals, contravariant cell velocity, face sign, chain-rule bed gradient; the symmetric filter of the curvilinear grid) |
| `sflow_2d_curv_sed_river_wavy` | 1 | SFLOW river reach over a steep sand mound on the wavy, skewed grid: discharge inflow, stage outflow, Manning friction, bed load and the sand slide (S 90 1) with the curvilinear distances and area ratios |
| `nhflow_3d_curv_sed_river_rot30_mpi2` | 2 | the river with NHFLOW (hydrostatic, fixed dt) turned by 30 deg at UTM coordinates on two ranks: P, Q and the bed gradient across the partition, the slide on several ranks |
| `nhflow_3d_curv_sed_amr_adaptive` | 1 | sediment on the NHFLOW mesh refinement (turned grid, patches following the surface, regrid every 4 steps, sand slide): flux matching of the Exner fluxes, one dtsed, qb restriction and fill, the filter and the slide over all grids (reefamr_sed.h), the water level of the cells around the patches after a bed update |
| `sflow_2d_curv_sed_amr_lev2_mpi2` | 2 | sediment on the SFLOW mesh refinement: two levels in a static box on the turned grid, two ranks (remote fills of qb and the bed, remote flux matching) |
| `nhflow_3d_mb_sed_amr` | 1 | sediment on the NHFLOW multi-block grid (the four turned quadrants, steep bump, filter and sand slide, one level over the vertex): interfaces of the blocks and of the patches open for the filter and the slide; the same bed as one block |
| `sflow_2d_mb_sed_mpi2` | 2 | the quadrants with SFLOW and the weighted multi-direction slide (S 90 4) on two ranks: remote fills and the reverse exchange of the slid material (xrun_rev) |
| `fnpf_3d_amr_basin` (+ `_mpi2`, `_place2`) | 1/2 | FNPF AMR waves through a static box (composite Laplace); `_place2` G 40 2: remote parents in the FAC preconditioner |
| `fnpf_3d_amr_tow` (+ `_place2`) | 1/2 | FNPF AMR towed floating body (X 10 1, X 210), the zone follows it (G 12, G 2 4); `_place2` G 40 2: zone patch placed whole, body grids and triangle ownership |
| `fnpf_3d_amr_basin_sub` (+ `basin3_sub`, `basin_sub_place`) | 1/1/2 | FNPF waves with subcycling (G 7 1): level 0 alone, the patch steps with fills in time, the Laplace solve on the level window with fixed parent columns, the restriction after the steps; `basin3_sub` three levels (G 1 2), `basin_sub_place` 2 ranks with placed patches (G 40 1) |
| `fnpf_3d_amr_tow_sub` | 1 | FNPF towed body with G 7 1: the finest level advances the body (psi_0 solve on its patches), level 0 with a predicted copy, the zone margin widened to 12 level-0 cells |
| `fnpf_3d_amr_tow_x18`, `fnpf_3d_amr_tow_sub_x18` | 1 | `fnpf_3d_amr_tow` and `fnpf_3d_amr_tow_sub` with X 18 0.5: the running mean of chi removed next to the body on the base grid and on the refinement patches, with subcycling once per finest step (patch 0016) |
| `fnpf_3d_amr_decay_sub` | 1 | FNPF moored cylinder, heave decay with G 7 1: the free body, its added mass from all grids at the start of the level-0 step |
| `fnpf_3d_amr_break` (+ `_sub`) | 1 | FNPF AMR wave breaking (A 350 1) inside a static patch on a 1:10 slope with the static coastline (A 343 1), the patch kept off the coast: detection on the patch (set_local), the viscosity of the cells around it from the source grid, the implicit damping on the patch box after a fill of the cells around it; `_sub` with G 7 1 |
| `fnpf_3d_amr_coast_wd` (+ `_sub`, `_mpi2`, `_sub_place`) | 1/1/2/2 | the patch over the static coastline (G 30 1): fnpf_fsfbc_wd on the patch (wet mask and coastline distance from the fine bed, clamp, coastline relaxation), interpolation and restriction from wet cells only, with breaking; `_sub` G 7 1, `_mpi2` 2 ranks, `_sub_place` G 7 1 on 2 ranks with placed patches (G 40 1) |
| `fnpf_3d_amr_break_adapt` (+ `_sub`), `fnpf_3d_amr_coast_wd_adapt` | 1 | FNPF solution-adaptive refinement: the patches follow the breaking cells (G 23 1) and the surface jump (G 20), with G 30 1 also the shoreline (G 22 2), regrid every 4 steps: tagging, the regrid without a body, patch_wd_ini of new patches; `_sub` with G 7 1 |
| `fnpf_2d_amr_regular` (+ `_place_mpi2`, `2`, `_sub`) | 1/2/1/1 | FNPF AMR on a 2D grid (one cell in y): fnpf_2d_regular with a static patch over 12-28 m; the levels refine x only, the patch is one cell wide in y (patches 0077, 0078); `_place_mpi2` G 40 1 on 2 ranks, `2` two levels, `_sub` G 7 1; the patch-local and (one rank) level-0 solves exact, no pre-smoothing without a body (patch 0082) |
| `fnpf_2d_amr_break` | 1 | fnpf_3d_amr_break on a 2D grid: breaking A 350 1 in the patch, the static coastline |
| `fnpf_3d_amr_break_kennedy`, `fnpf_3d_amr_break_romero` | 1 | fnpf_3d_amr_break with the breaking of Kennedy (A 350 2, A 351 1, A 365 1.2) and of Romero (A 350 3) on every grid (patch 0079) |
| `fnpf_3d_amr_basin_picard` (+ `_mpi2`), `fnpf_3d_amr_basin_mixed` | 1/2/1 | fnpf_3d_amr_basin with Picard iterations on the sigma cross-derivatives (A 324 3, A 325 0.01; `_mpi2` with A 328 1 on 2 ranks) and with the cross-derivatives in the rows (A 328 1, N 10 1) (patch 0080) |
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
| `fnpf_2d_active_gen_wall` | 1 | FNPF active generation B 98 4 with a wall at x+: the reflected waves leave through the wavemaker (shallow-water absorption against the target, ramp) |
| `fnpf_2d_active_beach_b99_4` | 1 | FNPF active beach with the linear-theory profile B 99 4 (sinh(k h), level-wise ghost potential) |
| `fnpf_2d_combined_beach` | 1 | FNPF combined beach B 99 6: eta and surface potential relaxed to the low-passed state over the beach, active beach at x+ |
| `fnpf_2d_combined_gen` | 1 | FNPF combined generation B 98 6 (zone to target + low-passed deviation, active generation edge at x-) with B 99 6 |
| `fnpf_2d_twoband` | 1 | FNPF two-band active absorption at the wavemaker (B 98 7) and the beach (B 99 7), reference period B 528 2 s |
| `fnpf_2d_irregular_2nd_b_twoband` | 1 | 2nd-order irregular waves (B 92 33) with two-band active generation (bound-wave target) and two-band beach |
| `sflow_2d_stokes2` (+ `_mpi2`) | 1/2 | SFLOW relaxation generation + beach, Stokes 2nd (tutorial 8_1) |
| `sflow_2d_stokes2_lagoon_mpi2` | 2 | the wave tank with its free surface and bed in the LAGOON store instead of the VTP files (P 18 1) (needs patch 0019) |
| `cfd_2d_dambreak_amr_lagoon` | 1 | the dam break with one refined level, its volume and its refined patches in the LAGOON store instead of the VTU, .vtr and .vtm files (P 18 1): the AMR set `cfd_amr` of 3D blocks; the solution is that of `cfd_2d_dambreak_amr` (needs patch 0023) |
| `nhflow_2d_channel_kepsilon_fsf_lagoon` | 1 | the k-epsilon channel with its free surface, now with eddyv, kin and epsilon of the top layer, in the LAGOON store instead of the VTP files (P 18 1); the solution is that of `nhflow_2d_channel_kepsilon` (needs patch 0022) |
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
| `nhflow_2d_waves_current_doppler_stokes5` | 1 | 5th-order Stokes waves (B 92 5) on the following current with Doppler (B 530 2 10): the Stokes library refines k with its own dispersion relation omega = k (c_E + U_n) and recomputes its coefficients (needs patch 0001, NHFLOW Wave Current Scour Protection) |
| `nhflow_2d_waves_setup_heff` | 1 | linear waves on a 1 m set-up (B 514 eta) with B 530 1 10: k on h_eff |
| `nhflow_2d_tide_waves_doppler` | 1 | `nhflow_2d_tide_waves_beach` with B 530 2 10: k follows the tide level and current |
| `nhflow_2d_irregular_current_doppler` | 1 | JONSWAP irregular waves (B 92 31, 20 components) on a following current with B 530 2 10: k of every component, cached phases rebuilt |
| `nhflow_2d_irregular_current_doppler_b89` | 1 | `nhflow_2d_irregular_current_doppler` with decomposed precalc (B 89 1): B 530 with B 89 1, the space parts of the zone decomposed again when B 530 sets a new k |
| `nhflow_2d_irregular_current_blocked` | 1 | the same against an opposing current of 1 m/s: components above omega = g/(4 abs(U)) blocked and faded out |
| `nhflow_2d_tide_clamp_level` | 1 | `nhflow_2d_tide_flather` with a clamped level edge (B 520 method 5) at x-; edge mass balance (B 529) |
| `nhflow_2d_tide_clamp_q` | 1 | `nhflow_2d_tide_flather` with a clamped discharge edge (method 6) with the tide's velocity at x- |
| `nhflow_2d_river_clamped` (+ `_mpi2`) | 1/2 | river 40 m3/s in through a clamped discharge edge (B 525, no background), out through a clamped level edge; B 529 sums across ranks |
| `nhflow_2d_tide_waves_b530_off` | 1 | `nhflow_2d_tide_waves_beach` with B 530 0: the B 530 default turned off |
| `nhflow_2d_tide_waves_b89` | 1 | `nhflow_2d_tide_waves_beach` with irregular waves (B 92 31) and decomposed precalc (B 89 1) on the tide; B 530 on by default, the space parts decomposed again every B 530 N steps |
| `nhflow_2d_current_log_profile` | 1 | constant current 1 m/s through Riemann edges with a log-law profile (B 513 1 1 0.001), 20 sigma layers |
| `nhflow_2d_current_power_profile` | 1 | the same with the 1/7 power law (B 513 1 2 7) |
| `nhflow_2d_waves_current_log_profile` | 1 | linear waves on the log-law current: generation zone and beach with the profiled background, B 530 default (mode 2) |
| `nhflow_2d_awa_b99_4` | 1 | `nhflow_2d_awa` with B 99 4: the NHFLOW active beach (B 99 3 / 4) is an absorbing Riemann edge at x+ with still water outside |
| `nhflow_3d_source_global_frame` | 1 | 3D basin, B 105 20 5 3, no B 92 wave; source 2 given in the global frame (B 505 2 1), generated in a B 520 zone at x- |
| `nhflow_2d_irregular_two_sources_decomp` | 1 | `nhflow_2d_irregular_two_sources` with B 89 1: the components of both sources decomposed (wave_field::dspace / dtime) |
| `fnpf_2d_irregular_two_sources_decomp` | 1 | `fnpf_2d_irregular_decomp` plus an irregular source 2 (B 500, B 504): B 89 1 with several sources in FNPF |
| `nhflow_2d_b89_linear_fallback` | 1 | `nhflow_2d_two_sources` (linear) with B 89 1: falls back to B 89 0 with a message, bitwise as `nhflow_2d_two_sources` |
| `nhflow_2d_irregular_2nd_a` | 1 | NHFLOW 2D, 2nd-order irregular waves A (B 92 32, 32 components) in the relaxation zone: cached-point evaluation of the pair terms (`wave_lib_irregular_2nd_cache.h`) |
| `nhflow_2d_irregular_2nd_b` | 1 | as above with 2nd-order irregular waves B (B 92 33) |
| `nhflow_3d_irregular_2nd_dir` | 2 | `nhflow_3d_irregular_decomp` with directional 2nd-order irregular waves A and B 89 0 |
| `fnpf_2d_irregular_2nd_b` | 1 | FNPF 2D, 2nd-order irregular waves B (B 92 33): cached potential (`wave_fi_c`) and eta of the pair terms |
| `fnpf_2d_irregular_2nd_b_setdown` | 1 | `fnpf_2d_irregular_2nd_b` with B 531 1: mean set-down and return current of zero mass flux in the 2nd-order target (FNPF only; ignored in the other solvers) |
| `fnpf_2d_irregular_2nd_b_deep` | 1 | `fnpf_2d_irregular_decomp` (600 m, Hs 5 m) with 2nd-order irregular waves B and B 89 0: the corrected second-order theory in deep water (the former coefficients stopped with a HYPRE breakdown) |
| `nhflow_3d_two_edges_irregular_b89` | 2 | `nhflow_3d_two_edges` with irregular waves and B 89 1: zone sources (B 524) with decomposed precalc |
| `fnpf_3d_two_edges_irregular_b89` | 2 | as above for FNPF |
| `cfd_2d_zone_sources` | 1 | CFD 2D flume, B 520 generation zone at x- with the B 92 wave and a linear source 2 (B 524 1 1, B 524 1 2) |
| `cfd_2d_zone_sources_b89` | 1 | CFD, two irregular sources with B 89 1; the zone selects source 2 only, so the B 92 components get space parts 0 there |
| `sflow_2d_zone_sources` | 1 | SFLOW 2D flume, Stokes 2nd wave plus a linear source 2 in a B 520 zone (B 524): precalc over the zone columns |
| `sflow_2d_irregular_zone_sources` | 1 | SFLOW, two irregular sources in a B 520 zone: cached-point evaluation (wave_eta_c / wave_uvw_c) |
| `nhflow_2d_combined_beach` | 1 | combined beach B 99 6 (default B 526 = T / 4): relaxation to the low-passed state plus the absorbing edge at x+ |
| `nhflow_2d_combined_gen` | 1 | NHFLOW combined generation B 98 6: relaxation zone to target + low-passed deviation (B 527), Riemann edge with the waves at x-; combined beach B 99 6 |
| `fnpf_3d_waverecon_dir` | 2 | wave reconstruction (B 92 51) from a waverecon.dat with directions (column 4), zones at x- and y- |
| `nhflow_3d_irregular_dir_cached` | 2 | `nhflow_3d_irregular_decomp` with B 89 0 and a frequency-direction grid (B 84 1, B 136 1): cached kinematics, one vertical profile per distinct k |
| `fnpf_3d_two_edges_irregular` | 2 | `fnpf_3d_two_edges_irregular_b89` with B 89 0: cached potential (wave_fi_c) with zone sources |
| `sflow_2d_tide_flather` | 1 | SFLOW channel 3 km: progressive tide in through a Riemann edge at x-, out through a Flather edge with a still-water background at x+ (iowave_sflow_tide.cpp, open-edge fluxes in sflow_HLL) |
| `sflow_2d_tide_riemann` | 1 | the same with Riemann edges at both ends |
| `sflow_2d_tide_waves_beach` | 1 | SFLOW tide through Riemann edges, linear waves from a zone with background + waves (B 523), beach relaxing to the tide |
| `sflow_2d_riemann_waves` | 1 | SFLOW waves generated by a Riemann edge with source 1 (B 524), beach at x+ |
| `sflow_2d_river_clamped` (+ `_mpi2`) | 1/2 | SFLOW river 40 m3/s through a clamped discharge edge (B 525), clamped level out; edge mass balance B 529 (RK-stage fluxes) to REEF3D_SFLOW_Log |
| `sflow_2d_tide_oblique` (+ `_mpi2`) | 1/2 | SFLOW basin, oblique tide with open edges on all four sides (y edges flagged as outflow in DIVEMesh) |
| `sflow_2d_waves_current_doppler` | 1 | SFLOW linear waves on a following current of 1 m/s with Doppler (B 530 2 10) |
| `sflow_2d_awa_edge` | 1 | SFLOW active beach B 99 3 as the absorbing Riemann edge at x+ (as NHFLOW) |
| `sflow_2d_combined_beach` | 1 | SFLOW combined beach B 99 6: low-pass relaxation beach in front of the absorbing edge |
| `cfd_2d_irregular_cached` | 1 | CFD irregular waves (B 92 31) with B 89 0: relaxation over the zone columns, cached-point evaluation (cell centre, u and v point per column) |
| `nhflow_3d_irregular_2nd_grid` | 2 | `nhflow_3d_irregular_decomp` with 2nd-order irregular waves A on a frequency-direction grid (8 x 4, B 136 1) and B 89 0: same-direction runs and pairs of other directions per distinct \|K\| |
| `cfd_2d_irregular_2nd` | 1 | CFD 2D flume, 2nd-order irregular waves A (B 92 32, 32 components) with B 89 0: pair terms at the cached zone-column points |
| `nhflow_2d_irregular_2nd_lowrank` | 1 | `nhflow_2d_irregular_2nd_a` with 128 components: the pair sums of the direction from low-rank expansions of their coefficients (`lr_build`, `lr_kin`, `lr_eta` in `wave_lib_irregular_2nd_cache.h`) |
| `fnpf_2d_irregular_2nd_b_lowrank` | 1 | `fnpf_2d_irregular_2nd_b` with 128 components: eta and the potential of the sum terms low-rank, the potential of the difference terms direct |
| `fnpf_2d_irregular_2nd_b_twoband_lowrank` | 1 | `fnpf_2d_irregular_2nd_b_twoband` with 128 components: B 98 7 / B 99 7 against the low-rank bound-wave target |
| `nhflow_2d_irregular_2nd_combined_lowrank` | 1 | `nhflow_2d_irregular_2nd_lowrank` with B 98 6 / B 99 6: the low-rank target in the zone-column relaxation and at the generation edge |
| `nhflow_3d_irregular_2nd_rows` | 2 | NHFLOW 3D, long-crested 2nd-order irregular waves A (128 components) along x: the cells across the crests share their phases, eta once per row, u, v, w from the Chebyshev tables in z of the row (`wave_lib_zgroup.h`) |
| `nhflow_3d_irregular_rows` | 2 | `nhflow_3d_irregular_2nd_rows` with 1st-order irregular waves (64 components): groups and tables in z in `wave_lib_irregular_1st` |
| `fnpf_3d_irregular_2nd_b_rows` | 2 | FNPF 3D, long-crested 2nd-order irregular waves B (128 components): eta and the surface potential once per row of zone cells; the difference-term potential of the low-rank block with the pair phase in the loop |
| `cfd_2d_channel_kepsilon` (+ `cfd_2d_channel_komega_mpi2`) | 1/2 | open channel, discharge inflow (B60 1) with the equilibrium k/ε/ω inflow profile, k-ε / k-ω across a rank border in x |
| `cfd_2d_channel_komega_t36` | 1 | k-ω free-surface damping T36 3 (y' = T37 h from the local water depth, dimensionless weight) |
| `cfd_2d_stillwater_plic_t41` | 1 | PLIC VOF still water, k-ω with T41 1: no NaN from the limiter at S = 0 |
| `nhflow_2d_channel_kepsilon` | 1 | NHFLOW open channel, discharge inflow with the equilibrium turbulence profile, k-ε, bed roughness A519 |
| `nhflow_2d_channel_suspended` | 1 | NHFLOW open channel with suspended load only (S 11 0, S 12 1), k-ε: conservative D·C scheme with FEx face fluxes, sediment leaving through the outflow, bed exchange only in the erodible region S 71 (needs Suspended Sediment patches 0001-0009) |
| `nhflow_2d_beach_suspended_mpi2` | 2 | NHFLOW waves on a 1:12 sand beach, bedload + suspended load, swash wetting/drying (dry-column deposit, Exner with the c_be the water column used), S 38, FEx across the MPI interface (needs Suspended Sediment patches 0001-0009) |
| `sflow_2d_sediment_basin` | 1 | SFLOW closed basin, corner dam break over an MPM sand bed, default Exner (S 31 2, FOU, bed filter): no bedload through the walls, SFLOW sediment start-up (needs Exner Review patches 0001-0016) |
| `sflow_2d_sediment_basin_noneq_mpi2` | 2 | as above with S 31 1, WENO flux Exner (S 32 4) and non-equilibrium bedload (S 33 1) across the subdomain boundary (needs Exner Review patches 0001-0016) |
| `nhflow_3d_sediment_pier_mpi2` | 2 | NHFLOW closed basin with a pier, default Exner S 31 2 with closed faces at the structure, NHFLOW bed shear with the thin-water clamp, sand slide S 90 1 (needs Exner Review patches 0001-0016) |
| `nhflow_3d_sediment_pier_rk2` | 1 | as above with the fully coupled RK2 sediment step S 10 12 (bed advanced in both stages of every flow step, dtsed = S 17 dt) (needs Exner Review patches 0001-0016) |
| `nhflow_3d_sediment_pier_s1610` | 2 | as `nhflow_3d_sediment_pier_mpi2` with k-omega (A 560 2) and the bed shear stress S 16 10 (shared with CFD S 16 10 and the CPM bedload layer): log law 2 cells above the bed on the sigma grid (S 63 2.0), at most half the water depth, max(log law, 0.3 rho k) with k of the first cell (patch 0021) |
| `cfd_2d_sediment_flume` | 1 | CFD 2D flume, MPM sand bed with a submerged block, k-ω, bed shear S 16 1 at the velocity sample height, CFD topo update, sand slide (needs Exner Review patches 0001-0016) |
| `cfd_2d_sediment_flume_s1610` | 1 | As `cfd_2d_sediment_flume` with the bed shear stress S 16 10 (`bedshear::column_tau`, shared with the CPM bedload layer): log law 2 cells above the bed level (S 63 2.0), max(log law, 0.3 rho k) with k of the first fluid cell, flow over and behind the block (patch 0020) |
| `cfd_2d_sediment_flume_s1611` | 1 | As `cfd_2d_sediment_flume` with the bed shear stress S 16 11 (bedshear::column_tau with two velocity scales, sqrt(log law * 0.3 rho k)), log law 2 cells above the bed level (S 63 2.0) (patch 0026) |
| `cfd_2d_sediment_flume_susp` | 1 | As `cfd_2d_sediment_flume` with the Eulerian suspended load (S 12 1, d = 0.15 mm, W 10 0.0015): erosion w_s c_be and deposition w_s c_1 as the only bed exchange, faces of the concentration to the bed, the block, the air and the walls closed; the concentration is checked through the bed (patch 0018; reference with patch 0018) |
| `cfd_2d_bedshear_wave_s16_9` | 1 | CFD 2D oscillatory flow in a periodic channel with a slip lid (B 180/181 body force, U = 0.3 sin(pi t) m/s), k-ω, MPM bed (S 35 50): wave friction bed shear S 16 9 with S 58 (needs Bed Shear Stress patch 0002; reference with 0001-0003) |
| `cfd_2d_bedshear_wave_current_s59` | 1 | `cfd_2d_bedshear_wave_s16_9` with a current (W 20) and S 59 1 20: S 16 9 as max(log law, Soulsby DATA13 wave-current stress), current from the running mean, 20 deg phase lead of the wave stress (needs Bed Shear Stress patch 0004) |
| `nhflow_2d_bedshear_wave_s59` | 1 | NHFLOW 2D Stokes wave flume, 10 σ-layers, k-ω, sand bed with Exner x10: S 16 9 with S 59 2 15 (DATA13 only, 15 deg phase lead) at the NHFLOW call site (needs Bed Shear Stress patch 0004) |
| `cfd_2d_bedshear_s16_7` | 1 | CFD 2D periodic channel driven by W 20, MPM bed: bed shear S 16 7 in 2D, lint*_2D keeps the caller's j (needs Bed Shear Stress patch 0001; reference with 0001-0003) |
| `cfd_2d_bedshear_bed_at_centre` | 1 | As `cfd_2d_bedshear_s16_7` with S 16 1 and the bed S 57 exactly at a cell-centre height: solid forcing with phi = 0 on faces (needs Bed Shear Stress patch 0003; reference with 0001-0003) |
| `nhflow_3d_cylinder_kepsilon_mpi2` | 2 | NHFLOW 3D channel with a cylinder (A580), k-ε, ranks split in y: k/ε and ν_t across the rank interface |
| `sflow_1d_channel_ke` (+ `_kw`) | 1 | SFLOW depth-averaged k-ε / k-ω (A260 1/2): k, ε/ω relax to the Rastogi–Rodi equilibrium |
| `sflow_2d_channel_walls_kw_mpi2` | 2 | SFLOW k-ω with side walls, ranks split in y: production at wall cells, uniform k/ω across the width |
| `cfd_3d_heave_sphere_6dof_kepsilon` | 4 | k-ε at a direct-forcing body: wall functions at gcdf4 cells, flagsf4 turn-off, solid forcing for k/ε |
| `cfd_2d_channel_veg_kepsilon` | 1 | vegetation box (B 310, B 308 0): drag ½ Cd a \|u\| u_i, Lopez & Garcia k/ε sources |
| `cfd_2d_freesurface_komega_t45` | 1 | k-ω buoyancy term T 45 1: implicit sink, k ≥ 0 at the interface |
| `cfd_2d_channel_les_t21_2` | 1 | LES Smagorinsky with the second-order high-pass filter T 21 2 |
| `nhflow_3d_channel_walls_kepsilon_mpi2` | 2 | NHFLOW side walls with A519 2: friction and turbulence wall functions per wall face (nhflow_wall.h), ranks split in y |
| `nhflow_3d_channel_walls_sst_mpi2` | 2 | as `nhflow_3d_channel_walls_kepsilon_mpi2` with k-ω and the SST limiter (A564 2): F2 wall distance to the bed and the side walls (walldist_h), ranks split in y |
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
| `seastate_2d_phase7_composite` | 1 | REEF3D::SEASTATE composite sweep across the refinement levels (A 797 1, Phase 7) on the island case with A 795 and A 796 2; validated in Dropbox SEASTATE/validation/28; Phase 17: the second ring layer of the patches for A 796 2 (second order where two patches of a level meet), re-blessed |
| `seastate_2d_phase7_composite_mpi2` | 2 | REEF3D::SEASTATE composite sweep (A 797 1) on the nonstationary island case with wind and DIA, patches cut at the rank boxes; Phase 17: the composite sweep pipelined across the ranks (A 781 1), `same_as` `seastate_2d_phase17_composite_ns` (bitwise) |
| `seastate_2d_phase7b_threads_mpi2` | 2 | REEF3D::SEASTATE 2 threads per rank (A 798 2, Phase 7b, wavefront order) on the refraction case with A 795 and A 796 2; the result of the serial sweep; validated in Dropbox SEASTATE/validation/30 |
| `seastate_2d_phase7b_sector` | 1 | REEF3D::SEASTATE fine direction sector (A 715 330 70 3, Phase 7b) on the handover case: bins of different width, handover spectra on uniform directions; validated in Dropbox SEASTATE/validation/32 |
| `seastate_2d_phase7b_windsea` | 1 | REEF3D::SEASTATE source iterations per cell (A 738 8, A 739 1e-4), distance-to-solution convergence test (A 799 1) and 2 threads on the fetch-limited wind-sea case; validated in Dropbox SEASTATE/validation/31 |
| `seastate_2d_phase6_diffraction` | 1 | REEF3D::SEASTATE phase-decoupled diffraction (A 718 1, Phase 6a, as SWAN) behind a semi-infinite breakwater (bathymetry raster); validated in Dropbox SEASTATE/validation/34 |
| `seastate_2d_phase6_diffraction_mpi2` | 2 | REEF3D::SEASTATE diffraction per frequency (A 718 2) on 2 ranks with 2 threads per rank: smoothed energy and Ca across the rank border |
| `seastate_2d_phase6_obstacles` | 1 | REEF3D::SEASTATE line obstacles (A 722, Phase 6a): oblique obstacle with transmission and specular reflection, Goda dam; validated in Dropbox SEASTATE/validation/35 |
| `seastate_2d_phase6b_coasts` | 1 | REEF3D::SEASTATE reflecting coasts (A 723, A 724, Phase 6b): oblique coast on the staircase of cells, diffuse reflection cos^4, a reflecting obstacle with diffuse reflection (A 726), spectral sparsity A 795; validated in Dropbox SEASTATE/validation/36 |
| `seastate_2d_phase6b_coasts_mpi2` | 2 | REEF3D::SEASTATE reflecting coasts on 2 ranks with 2 threads per rank and second-order fluxes (A 796 2, first order next to blocked faces): coastline directions and the active-cell mask exchanged across the rank border |
| `seastate_2d_phase6b_structures` | 1 | REEF3D::SEASTATE structures (A 725, Phase 6b): d'Angremond breakwater, Kt and Kr per frequency from a file, porous breakwater per frequency with diffuse reflection; validated in Dropbox SEASTATE/validation/37 |
| `seastate_2d_phase6b_vegetation` | 1 | REEF3D::SEASTATE vegetation (A 754 1, Phase 6b, as SWAN VEGETATION 1) in a flume; validated in Dropbox SEASTATE/validation/38 |
| `seastate_2d_phase6b_vegetation_field_mpi2` | 2 | REEF3D::SEASTATE vegetation per frequency (A 754 2, as SWAN VEGETATION 2) with the stems from seastate-vegetation.dat (A 756 1) on 2 ranks |
| `seastate_2d_phase6b_amr` | 1 | REEF3D::SEASTATE Phase 6 on the mesh refinement (Phase 6b): diffraction (A 718 1) on every level, reflecting coasts and an obstacle on the patch, composite sweep; validated in Dropbox SEASTATE/validation/39 |
| `seastate_2d_phase8_autosector` | 1 | REEF3D::SEASTATE automatic fine direction sector (A 716 2.5 15, Phase 8) around the Mitsuyasu boundary spectrum of the handover case, with spectral sparsity (A 795); validated in Dropbox SEASTATE/validation/41 |
| `seastate_2d_phase9_westh` | 1 | REEF3D::SEASTATE van der Westhuysen deep-water physics (A 732 2, Phase 9, as SWAN GEN3 WESTH): Yan wind input and saturation-based whitecapping on the fetch-limited wind-sea case; validated in Dropbox SEASTATE/validation/43 |
| `seastate_2d_phase9_width` (+ `_mpi2`) | 1 (2) | REEF3D::SEASTATE refinement by water width (A 762 4, Phase 9): open sea and a fjord 3 level-0 cells wide from the bathymetry raster, one level along the fjord only; on 2 ranks the wet cells gathered from both; validated in Dropbox SEASTATE/validation/44 |
| `seastate_2d_phase9_fas` (+ `_mpi2`) | 1 (2) | REEF3D::SEASTATE FAS coarse-grid correction (A 758 2 4, Phase 9) on the fjord case of `seastate_2d_phase9_width`, tight convergence for several corrections; validated in Dropbox SEASTATE/validation/45 |
| `seastate_2d_phase10_srcit_rows` | 1 | REEF3D::SEASTATE source iterations (A 738 16, A 739 1e-4) on the fetch-limited wind-sea case on 6 rows (Phase 10 fix: the action density limiter from the spectrum of the cell before its first solve, non-contracting cells keep the first solve); validated in Dropbox SEASTATE/validation/46 |
| `seastate_2d_phase10_sector` | 1 | REEF3D::SEASTATE computational direction sector (A 717 270 90, Phase 10, as SWAN CGRID SECTOR) on the fetch-limited wind-sea case: 13 of 24 directions, one quadrant with one direction and one without any; validated in Dropbox SEASTATE/validation/47 |
| `seastate_2d_phase10_fas2` (+ `_mpi2`) | 1 (2) | REEF3D::SEASTATE FAS (A 758 2 4) with two refinement levels (G 1 2, A 762 8) on the fjord case of `seastate_2d_phase9_width` (Phase 10); validated in Dropbox SEASTATE/validation/48 |
| `seastate_2d_phase10_fas_ns` | 1 | REEF3D::SEASTATE FAS (A 758 1 4) in a nonstationary run (5 steps of 120 s, 2 iterations per step) on the fjord case (Phase 10); validated in Dropbox SEASTATE/validation/48 |
| `seastate_2d_phase10_threads_amr`, `seastate_2d_phase10_threads_composite` | 1 | REEF3D::SEASTATE 2 threads (A 798 2, Phase 10) on the refinement, level by level (bitwise the reference of `seastate_2d_phase9_width`) and in the composite sweep (A 797 1; bitwise the serial composite sweep); validated in Dropbox SEASTATE/validation/48 |
| `seastate_2d_phase12_fas_uniform` (+ `_mpi2`) | 1 (2) | REEF3D::SEASTATE FAS V-cycle on the uniform grid (A 758 2 4 without G 1, Phase 12) on the fetch-limited wind-sea case: 6 coarse grids of 2 x 2 cells; on 2 ranks every rank coarsens its own cells, coarse halo exchange; validated in Dropbox SEASTATE/validation/52 |
| `seastate_2d_phase12_fas_uniform_threads` | 1 | as `seastate_2d_phase12_fas_uniform` with 2 threads (A 798 2, Phase 12 worker threads) on the fine, residual and coarse sweeps; bitwise the reference of `seastate_2d_phase12_fas_uniform` |
| `seastate_2d_phase12_fas_uniform_forcing_mpi2` | 2 | REEF3D::SEASTATE FAS on the uniform grid in a nonstationary run (A 758 1 4, Phase 12) with a wind field (A 730 2) and boundary spectra per side cell (A 711 3) restricted to the coarse grids |
| `seastate_2d_phase12_threads_patches` | 1 | REEF3D::SEASTATE `seastate_2d_amr_static` with four more boxes (4 patches on level 2) and 2 threads (A 798 2, Phase 12): the patches of a level shared out to the threads along the dependency graph of their ring reads; bitwise the same case with one thread; validated in Dropbox SEASTATE/validation/54 |
| `seastate_2d_phase13_tiles` | 1 | REEF3D::SEASTATE 3 threads (A 798 3, Phase 13) with spectral sparsity (A 795) and second-order fluxes (A 796 2) on the refraction case: the sweeps over tiles of cells along their dependency graph, the first visits and the stationary Hs test threaded; bitwise the same case with one thread; validated in Dropbox SEASTATE/validation/57 |
| `seastate_2d_phase14_refraction_mpi4` | 4 | REEF3D::SEASTATE sweeps pipelined across the MPI ranks (A 781 1, Phase 14): `seastate_2d_refraction` on 4 ranks, `same_as` it (`regression.py ranks`: bitwise, cell by cell); validated in Dropbox SEASTATE/validation/58 |
| `seastate_2d_phase14_tiles_mpi3` | 3 | `seastate_2d_phase13_tiles` (3 threads, sparsity, second order) on 3 ranks, pipelined: two halo layers per message, the second halo layer of the rank boundaries wet (seastate_halo2), the first visits exchanged; `same_as` the one-rank case |
| `seastate_2d_phase14_front_mpi4` | 4 | `seastate_2d_front` (nonstationary) on 4 ranks, pipelined: one iteration per step as on one rank; `same_as` it |
| `seastate_2d_phase14_current_mpi2` | 2 | `seastate_2d_current` on 2 ranks, pipelined: the current carries energy against the sweep (the downstream halo of the last exchange); `same_as` it |
| `seastate_2d_phase14_sordup`, `_windsea`, `_diffraction`, `_coasts`, `_vegetation_field`, `_wind_field` | 1 | the one-rank twins of `seastate_2d_phase7_sordup_mpi2`, `seastate_2d_windsea_mpi2`, `seastate_2d_phase6_diffraction_mpi2`, `seastate_2d_phase6b_coasts_mpi2`, `seastate_2d_phase6b_vegetation_field_mpi2`, `seastate_2d_wind_field_mpi2`, which have `same_as` them (Phase 14); also `refraction_mpi2`, `front_mpi2`, `init_mpi2`, `phase7b_threads_mpi2` (`same_as` `phase13_tiles`) |
| `seastate_2d_phase15_fas_coasts` | 1 | REEF3D::SEASTATE FAS on the uniform grid (A 758 2 4) with reflecting coasts (diffuse), a reflecting line obstacle and the spectral sparsity (A 795): coarse faces of their own, no correction under blocked fine faces, zero bins held at zero (Phase 15); validated in Dropbox SEASTATE/validation/61 |
| `seastate_2d_phase15_fas_structures` | 1 | REEF3D::SEASTATE FAS on the uniform grid (A 758 2 4) with the three structure types (d'Angremond, per frequency, porous) and diffuse reflection: the coarse coefficients from the restricted spectra (Phase 15); validated in Dropbox SEASTATE/validation/61 |
| `seastate_2d_phase16_faces` | 1 | REEF3D::SEASTATE refinement box (G 1 1) whose east and north edges lie on the rank borders of the `_mpi2` twin: swell over a shoal, stationary to A 708 1e-7, handover points downstream of the patch (Phase 16); validated in Dropbox SEASTATE/validation/65 |
| `seastate_2d_phase16_faces_mpi2` | 2 | `seastate_2d_phase16_faces` on 2 ranks: the level-0 cells next to the patch of the other rank take the faces of its covered cells (covered halo cells, faces sent after every halo exchange); `same_as` the one-rank case with `same_as_rtol` 1e-5 (before Phase 16: 0.7 % in Hs) |
| `seastate_2d_phase16_faces_levels` | 1 | two levels from the depth criterion (A 791) over the shoal, the levels cut into several patches on one rank: ring cells next to a covered cell of a neighbouring patch of the level take its face (Phase 16; results changed against Phase 15 by up to 0.17 % in Hs) |
| `seastate_2d_phase16_faces_levels_mpi4` | 4 | `seastate_2d_phase16_faces_levels` on 4 ranks: the ring faces of covered source cells on other ranks with the fills, the level-0 faces across the rank borders; `same_as` the one-rank case with `same_as_rtol` 1e-5 (before Phase 16: 0.25 %) |
| `seastate_2d_phase16_faces_composite`, `_composite_mpi2` | 1, 2 | `seastate_2d_phase16_faces_levels` with the composite sweep (A 797 1): the ring faces from the latest fine cells (ghost_fill); on 2 ranks the composite sweep pipelined (Phase 17), `same_as` the one-rank case, bitwise (Phase 16: within 1e-5 with the lagged sweeps) |
| `seastate_2d_phase17_composite_mpi4` | 4 | `seastate_2d_phase16_faces_composite` on 4 ranks: the composite sweep pipelined across the ranks (A 781 1, Phase 17), the boundary tiles carry the faces of the covered level-0 cells and the ring cells and ring faces of the patches downstream; `same_as` the one-rank case, bitwise; validated in Dropbox SEASTATE/validation/66 |
| `seastate_2d_phase17_sordup`, `_sordup_mpi2` | 1, 2 | the same with A 795 1e-4 and A 796 2, 30 iterations: the second ring layer of the patches (a cell two upwind beyond the box used where it is a cell of the same level not covered by a finer grid), the covered cells of the second halo layer marked; the `_mpi2` case pipelined, `same_as` the one-rank case, bitwise |
| `seastate_2d_phase17_windsea`, `_windsea_mpi4` | 1, 4 | the same with a wind sea (Komen, DIA), 30 iterations; on 4 ranks with 2 threads each (A 798 2), pipelined, `same_as` the one-rank case with one thread, bitwise |
| `seastate_2d_phase17_fas`, `_fas_mpi2` | 1, 2 | `seastate_2d_phase17_windsea` with FAS (A 758 2 4): on 2 ranks the composite sweeps pipelined, the coarse sweeps of level 0 pipelined as on a uniform grid; `same_as` the one-rank case, bitwise |
| `seastate_2d_phase17_composite_ns` | 1 | `seastate_2d_amr_mpi2` (nonstationary island case) on one rank with the composite sweep: the one-rank twin of `seastate_2d_phase7_composite_mpi2` |
| `seastate_2d_phase18_ckpt`, `_ckpt_mpi2` | 1, 2 | `seastate_2d_phase17_windsea` with A 795 1e-4, A 796 2, A 799 1: a checkpoint every 7 iterations of the stationary step (P 48 7, Phase 18: `REEF3D_Checkpoint/<count>_<iteration>`, the latest after iteration 28) and the load map of the measured work (G 45 1); the full runs of the restart cases; `_mpi2` `same_as` the one-rank case, bitwise; validated in Dropbox SEASTATE/validation/67 |
| `seastate_2d_phase18_restart_mpi2`, `_restart_mpi2to1` | 2, 1 | the restart (I 45 1, `chain`) from the latest checkpoint of `_ckpt` (1 rank) on 2 ranks and of `_ckpt_mpi2` on 1 rank: the spectra of all grids by global cell index, the active ranges of the sparsity, the last changes of Hs; the step continues at iteration 29; `same_as` `seastate_2d_phase18_ckpt`, bitwise |
| `seastate_2d_phase18_ckpt_uniform`, `_restart_uniform_mpi2` | 1, 2 | `seastate_2d_phase14_sordup` (uniform grid, second order) with A 795 1e-4, A 799 1, A 708 1e-8, P 48 2, and its restart on 2 ranks (two halo layers, pipelined); `same_as`, bitwise |
| `seastate_2d_phase18_ckpt_ns`, `_restart_ns_mpi2` | 1, 2 | `seastate_2d_phase17_composite_ns` (nonstationary island case) with P 48 4, and its restart after step 4 on 2 ranks; `same_as`, bitwise |
| `seastate_2d_phase18_rebalance_mpi2` | 2 | `seastate_2d_phase17_windsea` on 2 ranks with the rebalance (G 46 1.2, G 47 5, `"rebalance": 3`, DIVEMesh M 50 1): the measured load 1.45 stops the run after iteration 5 (checkpoint, load map, exit code 3), it continues on the weighted partition of the global grid (M 40 1, M 41 1), load 1.01; `same_as` `seastate_2d_phase17_windsea`, bitwise |
| `sflow_2d_bank_ediff` (+ `_idiff`) | 1 | SFLOW sloping bank (T 62) with a shoreline, constant viscosity, A212 1 / 2: free slip at the dry neighbours |
| `sflow_2d_bank_inflow_conveyance` | 1 | as `sflow_2d_bank_ediff`, discharge across the inflow face by conveyance (B 62 1), Cartesian (needs patch 0003 of Curvilinear/patches) |
| `cfd_2d_channel_kepsilon_stretched_ifou` | 1 | CFD open channel on a grid stretched in x (B 101 1, B 111 2.0), k-ε with implicit first-order upwind T 12 1: per-face conservative upwind |
| `nhflow_2d_bump_ediff` | 1 | NHFLOW 2D flow over a bump, constant viscosity, explicit momentum diffusion A 512 1: σ face metrics, A 513 wall rule (no no-slip bed from the A 518 2 ghost), viscous time-step limit |
| `nhflow_3d_bank_kepsilon` | 1 | NHFLOW 3D channel with a sloping bank (T 62) and a shoreline along x, k-ε: zero gradient of k and ε towards dry columns |
| `cfd_2d_patch_inflow_komega` | 1 | CFD channel with a discharge inflow patch (B 441, B 411) instead of B 60, k-ω: IO flags at the patch faces, turbulence zero gradient there |
| `cfd_2d_patch_pressure` | 1 | as `cfd_2d_patch_inflow_komega` with an outlet patch on the x+ wall (B 441, B 413): hydrostatic ghost pressure below the patch water level, level set from the water level, pressure correction Dirichlet |
| `sflow_1d_channel_patch` | 1 | `sflow_1d_channel_ke` with patches instead of B 60: discharge inlet (B 440, B 411) and water level outlet (B 440, B 413) on walls (patchBC_2D) |
| `nhflow_2d_bore_beach_kepsilon` | 1 | NHFLOW 2D bore running up a dry 1:20 beach, k-ε: newly wetted columns start from the wet neighbours, length-scale limit l ≤ κh |
| `nhflow_3d_boom_particles_mpi2` | 2 | NHFLOW river plastic interceptor: angled boom (B 350, skirt depth below the moving free surface, normal-only resistance), catamaran hulls and a permeable conveyor as surface-piercing drag boxes (B 351), floating and buoyant particles blocked above the draft and passing under it (L 71), capture zone at the conveyor (L 72), boom/hull force output (REEF3D_NHFLOW_Boom); particles cross the rank interface along the boom (needs River Interceptor patches 0001-0002) |
| `nhflow_3d_dem_boom_mat` | 1 | NHFLOW + DEM: floating debris parcels (light and heavy) driven against a boom that acts on the flow as B 350 and on the parcels as a fixed DEM box; contact loads on fixed particles (REEF3D_DEM/REEF3D-DEM-Fixed-Loads.dat) and heave/roll damping of floating unresolved particles (E 31) (needs River Interceptor patches 0001-0004) |

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
- The wave gauges (P 51) of FNPF, NHFLOW and SFLOW, and the SFLOW AMR gauge file, are interpolated at the
  gauge between the cell centres (`ccslipol4`, REEFAMR patch 0085); the FNPF and SFLOW references with gauges
  were re-blessed with it (the cell that holds a gauge on a node lies dx/2 downstream).
- The runner sets `REEF3D_AMR_CHECK=1` (unless set in the environment): the grids of a mesh refinement
  check the contract of their runtime state with level 0 (`src/lexer_amr.cpp`, patch 0050) on every step and stop
  the run (exit code 3) on a violation; `REEF3D_AMR_CHECK=0` prints the violation only.

## Unit tests

`unit/` holds standalone tests of solver-independent kernels (no MPI, no REEF3D binary); the build
line is at the top of each file, run from `unit/`:

| test | covers |
|---|---|
| `rigidbody_test.cpp` | 6DOF rigid-body core: quaternion/Euler, constant force, oscillator and torque-free top (orders of RK2/RK3/RKLS3/RK4), DOF modes, damping |
| `ship_test.cpp` | ship module kernels: waterline clipping, wetted surface, draft strips, cross-flow drag, ITTC-1957 line, roll damping, propeller KT/KQ, actuator disk (discrete force and torque, swirl sense), MMG rudder (signs, symmetry, slipstream, course stability, f_alpha, asymmetric gamma_R), MMG hull polynomials, MMG wake |
| `fem_test.cpp` | FEM solid solver: cantilever (Timoshenko, frequency), objectivity, J2, crack band energy, contact, collapse, STL snapping, patch test, presets, settling/check, structural damping, rigid bodies, walls and inclined bed (stick/slide), debris impact (peak u √(kM), duration, restitution, crushing, two bodies in series, deformable wall, added mass not in the impact, crushing set per contact point), collapse with a rigid fragment, reinforced concrete (bars in the skin, reinforced tie: ρσ_s and rupture, RC cantilever pushover against a section analysis, plain concrete brittle), adaptive integration (RC cantilever as with full integration, elastic top on one point), threads (1 and N OpenMP threads give bitwise the same results; build with -fopenmp), ranks (the solid split over 3 ranks, threads with a shared-memory communicator in place of MPI, gives bitwise the state of one rank; link with -lpthread) |
| `lagoon_store_test.cpp` | LAGOON store writer (P 18, needs `-lz -lpthread`): VTU header parsing, σ-level offsets, shard files read back (append, CRC-32C index, inner chunks, components), Cartesian (CFD) blocks split in z, noisy values compressed by several threads (the same shard with 1 and 4 threads, chunks read back, fill chunks left out), the case's georef.json in the store's root (on one line; left out when it is not a JSON object), format version 0.5, a curvilinear surface's `xc` and `yc` (shape, x fastest, float64; none for a rectilinear grid; needs patch 0020), a new run's first writer moves the store of the earlier run aside, named by its run id, and later writers of the run keep the new store (SFLOW: surfaces and AMR without a volume; needs patch 0021) |
| `lagoon_bodies_test.cpp` | LAGOON body writer (P 18, needs `-lz`; built with `../../src/lagoon_store.cpp`): quaternion of REEF3D's rotation matrix, two rigid bodies over 8 outputs (set and body attributes, mesh once, motion arrays, time counted once every body has it), a body off its rigid motion refused |
| `lagoon_particles_test.cpp` | LAGOON particle writer (P 18, needs `-lz`; built with `../../src/lagoon_store.cpp`): 5 outputs of 0 to 700000 particles (inner chunks and shards crossed), field names and their arrays, int32 fields, values read back from the shard files; objects with cells (ice floes: a cell set stored only when the cells change, a value per cell), inconsistent cells refused (needs patch 0009) |
| `lagoon_amr_test.cpp` | LAGOON AMR writer (P 18, needs `-lz`; built with `../../src/lagoon_store.cpp`): 4 outputs of 3 and 4 grids (level 0 of two ranks, a moving patch, a second one from the second output), set attributes and an int32 field, format version 0.5, grid_end and the grid table, coordinates as float64, levels, coordinates and cell values read back from the shard files; grids packed for MPI and back, a cut buffer refused; a grid with too few values refused and the set stopped (needs patch 0013); CFD's refined blocks as a set of 3D grids: `dims` 3, sizes nx, ny, nz, `z` and `grids/z_end`, every cell, z and block values read back, blocks packed and unpacked with z, a 2D grid in the set refused (needs patch 0023) |
| `runlog_test.cpp` | run log (`REEF3D_Case/REEF3D_<SOLVER>_run.jsonl`) of outputs that go only into the LAGOON store (P 18 1; built on its own, `runlog.cpp` with a small lexer stand-in): the run line written lazily with the first output, a `zarr` stream with store folder, group and pieces registered once, output lines with the time and iteration of the output (the store counts the volume at the next one), first_step, the end line (needs patch 0014) |
| `geometry_test.cpp` | horizontal geometry layer of the grid (`grid_geometry.cpp`): Cartesian grids (uniform, stretched, UTM coordinates, 2D) give 2D nodes equal to XN/YN and metrics equal to the 1D products bit for bit, closure exactly zero; curvilinear grids (annulus sector, sheared, meandering corridor): metrics = quadrilateral formulas, face vectors close every cell, areas sum to the boundary polygon area; bilinear subdivision (AMR patches) equals the 1D subdivision bit for bit on Cartesian nodes and tiles curvilinear cells exactly (area and face vectors) |
| `seastate_test.cpp` | REEF3D::SEASTATE kernels: spectral grid (logarithmic frequencies, bin widths, directions, direction faces, quadrants), block-sparse action storage (land tiles not allocated, tiles clipped at the range end, ghost cells, no overlap), integrated parameters (Hs, Tp, Tm01, Tm-1,0, direction, spread), dispersion relation (deep/shallow limits, residual, cg), SWAN 2D spectral file reader (CDIR/NDIR, VaDens/EnDens), memory budget, source terms (Battjes-Janssen Q_b, breaking with Newton linearisation, JONSWAP friction, moments with the spectral tail, Komen whitecapping, wind input with the Wu drag, DIA energy and action conservation, symmetry, shallow-water scaling and diagonal derivative, LTA energy conservation and Ursell limit, P >= 0 and D >= 0), surfbeat (single-frequency grid, Roelvink breaking, Herbers bound-wave coefficient vs. the LHS62 shallow-water limit, bichromatic envelope and bound wave, random-phase generator: variance, T_rep, directions, envelope mean, determinism), forcing files (date-times, series of SWAN spectra with many locations and times, wind field interpolation), bathymetry raster (reading, bilinear interpolation, cell mean and wet/dry majority rule), SWAN maximum energy (cap to (gamma d)^2/4, shape kept), compact direction sector (A 717: directions from the start of the sector, not periodic, lower faces, no neighbour beyond the ends, direction bins), the axis quadrant of a compact sector (Phase 15: two sweeps for 270 -> 90 and 180 -> 360), porous structure with strong resistance (Phase 15: Kt falls with H, the progressive mode). Build line with `seastate_test_threads.cpp` (a serial `run_threads` for `seastate_store`) |
| `coordinates_test.cpp` | world <-> model coordinate transform of rotated/shifted grids (`cms_flag` 1; built on its own, `coordinates.cpp` with a small lexer stand-in): `XYout(XYin(p))` and `XYin(XYout(p))` give `p` over a sweep of grid angles, origins (up to UTM size) and points, `XYin` agrees with `Xin`/`Yin` and `XYout` with `Xout`/`Yout`, a 90 deg rotation by hand, distances kept, segment directions and `Alpha_*_in/out` consistent, `cms_flag` 0 leaves points unchanged (needs patch CMS/0001) |
| `static_check.py` | static data in src and the REEF3D::CORAL folders below src/coral (no build, `python3 static_check.py`): every non-const static and thread_local variable (class members, function-local, file scope) against `static_allowlist.txt`, where each has the reason it holds no data of one grid (every grid of a mesh refinement has its own lexer, fdm and solver objects); a new static fails until it is moved into the object of the grid or listed with a reason, `--list` prints the current list (needs patch 0050 for the lexer_amr.cpp entry) |
| `coral_test.cpp` | REEF3D::CORAL core (`src/coral/core`, `structure`; built without REEF3D and MPI): name and version, rotations (exp/log, tangent operator, Euler angles), Lie group generalized-α (free fall exact, oscillator second order and without numerical damping for ρ∞ = 1, torque-free top), joints (pendulum on a revolute joint against the ODE, prismatic joints in series with the joint force, fixed joint, heavy top on a spherical joint against Euler's equations about the pivot, cylindrical, planar), tether through slack phases, linear 6x6 stiffness, modified Newton (needs patches CORAL/0001-0002) |
| `coral_hydro_test.cpp` | REEF3D::CORAL hydrodynamics (`src/coral/hydro`): waves (dispersion, JONSWAP Hs, kinematics against the elevation), retardation function against an analytic kernel, A∞ from Ogilvie's relation, state-space realization, database file round trip, WAMIT reader, heave in regular waves against the frequency-domain RAO (convolution and state-space), decay (damping ratio), hydrostatic pitch stiffness with the reference point away from the CoG, two coupled bodies against the 2-DOF frequency domain, Morison (inertia in deep-water waves, drag in a current, added mass of a free member); quadratic damping `drag` (heave decay against m z'' + Q |z'| z' + k z = 0, energy dissipated = drag work, keyword with Cd and area, a pitch coefficient) (needs patches CORAL/0003, 0029) |
| `coral_pto_test.cpp` | REEF3D::CORAL power take-off (`src/coral/pto`): linear damper (decrement, absorbed energy = loss of mechanical energy), Coulomb friction (stopping time and distance), end stop (penetration, rebound), latching and declutching, power and force limits, rotational PTO in a revolute joint, two-body point absorber (float, spar, PTO in the joint between them) in regular waves against the frequency domain (relative motion, mean power); the PTO element library shared with the REEF3D 6DOF bodies (X 500 - X 508): element Jacobians against finite differences, composite with force and power limits and the declutching window, zero-crossing interpolation, latching and declutching controllers, hydrodynamic table and B_opt, prismatic joint of one body (q, qd, J), latching with a stiff latch spring; joint drive (`drive ... sine`: a prismatic and a revolute joint follow the prescribed sine, the drive force against m q'' + k q, keyword); drive-train losses (loss damper and loss Coulomb friction, X 510 / `c_loss`, `coulomb_loss`: decay with c + c_loss, absorbed + dissipated energy = loss of mechanical energy, friction stop, not limited or declutched, B_opt with a loss, keywords, channels, checkpoint) (needs patches CORAL/0004, 0016-0017, 0024, 0044) |
| `coral_mooring_test.cpp` | REEF3D::CORAL mooring (`src/coral/mooring`): elastic catenary (hanging free, seabed and friction, profile, anchor forces, inextensible limit), quasi-static lines on a body, lumped-mass line relaxed against the catenary, surge decay with dynamic against quasi-static lines, drag of a taut line in a current (needs patch CORAL/0004); shared line between two bodies without a seabed and an anchor above the seabed (0031), wave kinematics once per step against every Runge-Kutta stage (0032), rigid motion of a line carried by its end nodes and the step independence of the fairlead tension (Hermite end motion, 0033), buoyancy of a line through the free surface, seabed friction along and across, the stop of an unstable line (0034), start from the static equilibrium of the segments (0035); mooring systems with free points: a clump between two lines, a buoy on a riser, a bridle, a pendulum period, a clump and a junction on the seabed (0036, 0041), the static start of a stiff chain (0037), line failure at a time and a tension and the restart after it (0038), a 10 % slope, a grid against a plane, the bed of a sampled fluid (0039), the bending of a cantilever, a taut line with EI and the clamping moment (0040), the damping law of slack segments (`BA_slack`: a compressed slack segment pushes by BA v / l, a taut line the same bit for bit under both laws, 0042), stick-slip seabed friction (a chain on a 1:5 side slope stays put where the regularized law creeps, a clump breaks away at mu_static mu W and slides at mu W, dragged lines under both laws, restart with the anchors bit for bit, 0043) (needs patches CORAL/0031-0043) |
| `coral_io_test.cpp` | REEF3D::CORAL input and output (`src/coral/io`): a model with every keyword read with the right numbers (continued lines, comments), errors that name the line, a pendulum from the input file bit for bit as the same model built in code, one output column per channel; structure keywords (`section`, `member`, `point_mass`, `static`, `modes`, `output vtp`) give bit for bit the same modes as the same frame built in code, `modes.dat` (kinetic energy split), `mode_shapes.dat` and the VTP file are written ; offshore keywords (`superelement`, `soil`, `soil_layer`, `pile`, `morison ... member` on frame and superelement members, `force`, `point_mass ... at`, `linear6 ... K` full matrix) build the expected loads and find the static equilibrium, superelement elements in the VTP; mooring `point`, `line_type` and `section`, `bathymetry plane|grid|fluid`, `BA_slack`, `mu_static`, `ks`, `cs` build the expected systems (needs patches CORAL/0005, 0010, 0015, 0036-0039, 0042, 0043) |
| `coral_beam_test.cpp` | REEF3D::CORAL co-rotational beams (`src/coral/structure`, static solver in `core`): cantilever under a tip load (Timoshenko and Euler-Bernoulli), axial and torsion; roll-up into a circle under the tip moment 2πEI/L (large rotations); Bathe's 45° bend against the published tip position; bending frequency with lumped mass in the time domain; free beam in rigid motion (no strain); beam end on a rigid body at an offset; analytic strain matrix B and element tangent (stiffness, damping) against finite differences in a deformed state; a lattice tower of 128 nodes with sparse LU gives the same static and dynamic response as dense LU (needs patches CORAL/0007-0012) |
| `coral_linear_test.cpp` | REEF3D::CORAL linear analysis (`src/coral/linear`): six uncoupled oscillators and a damped mode (damping ratio), pendulum on a revolute joint (geometric stiffness of the constraint force), cantilever bending, torsion and axial modes against Euler-Bernoulli, top mass against Rayleigh, NREL 5 MW tower with gravity (point mass and rigid rotor-nacelle assembly) against a planar reference, `Dropbox/Claude/CORAL/validation/02_nrel5mw_tower` (needs patches CORAL/0007-0010) |
| `coral_superelement_test.cpp` | REEF3D::CORAL Craig-Bampton superelement (`src/coral/structure/coral_superelement`) and generalized coordinates of the core: element matrices (rigid-body modes, mass), cantilever with a tip body: tip load (Timoshenko) and self-weight (q L^4/8EI) with all modes, fixed-interface frequency against the clamped-clamped beam, tip mass against the Euler-Bernoulli frequency equation with 3 modes, free vibration (period, Newton iterations per step as a check of the tangent), tripod superelement under a frame tower against the full frame (needs patches CORAL/0011-0012) |
| `coral_soil_test.cpp` | REEF3D::CORAL soil springs (`src/coral/soil`): API soft clay p-y (static, cyclic), API sand p-y (C1, C2, C3 against the API figure, initial slope), t-z, Q-z; spring tangent incl. shaft torsion against differences; elastic pile in a Winkler soil against Hetenyi (head load, head moment, axial), pile in API clay under a large lateral load (equilibrium) (needs patches CORAL/0011-0013) |
| `coral_offshore_test.cpp` | REEF3D::CORAL Morison loads on flexible members (`src/coral/hydro/coral_member_morison`, `core/coral_points`): buoyancy, inertia force on a stiff monopile against Cm rho A w^2 H/(2k), added mass of a submerged cantilever (frequency ratio), cantilever in a current as frame and superelement, point mass at an offset (compound pendulum) (needs patches CORAL/0011-0014) |
| `coral_fluid_test.cpp` | REEF3D::CORAL coupling to a resolved fluid (`src/coral/hydro/coral_fluid_load`, `core/coral_track`): Hermite interpolation of the recorded body motion (position, velocity, rotation, angular velocity) against an analytic motion; fluid load (polynomial in time through the samples, provisional sample, added-mass force); a heaving body driven by an emulated fluid solver against the damped oscillator: second order, quadratic extrapolation and the predictor at the error of CORAL's own time integration, stable at ten steps per period and for a body 100 times lighter than its added mass; keywords `fluid_body`, `fluid_coupling`; sampled kinematics (strip centres and line nodes, the `wave_field` outside the fluid and before the first sample, monopile Morison force from sampled linear waves second order in the sampling step); one-way loads (attachment point and moment, `fluid_member` strips with their total force and moment, cantilever tip under the strip load against Timoshenko, `excitation`, `order`); journal round trip (`record`, `coupling::replay` gives the same time series); a fluid without the instantaneous added mass (NHFLOW: the reaction of the body acceleration in its loads) with CORAL's artificial added mass: a body ten times lighter than its added mass unstable without, stable with 0.6-2 times the added mass, as accurate as the implicit added mass at 1, second order; a provisional sample inside the step (CFD, low-storage RK3 at 8/15 dt) is extrapolated at most linearly beyond it (needs patches CORAL/0011-0024) |
| `coral_tests.sh` | builds and runs all REEF3D::CORAL unit tests and the layering check (`./coral_tests.sh`, CXX and JOBS from the environment) |
| `coral_layering_check.py` | layering of REEF3D::CORAL (no build, `python3 coral_layering_check.py`): each folder of `src/coral/` includes only from the folders below it (table in `src/coral/README.md`), REEF3D headers and MPI only in `src/coral/reef3d/`, file names unique across the folders (needs patch CORAL/0001) |

