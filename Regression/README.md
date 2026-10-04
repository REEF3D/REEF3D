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
3. the normal text output (wave gauges, probes, forces, 6DOF) with a tolerance.

Result per case: **identical** (bitwise) · **close** (within `--rtol/--atol`) · **different** · **failed**.

## Quick start

```bash
cd Regression

# build two binaries with identical, reproducible flags (-O2 -ffp-contract=off, no -march=native/LTO)
./build_reef3d.sh -r hans_dev -o ~/reef3d_reg/ref      # a git revision (exported, working tree untouched)
./build_reef3d.sh             -o ~/reef3d_reg/new      # the working tree (incremental rebuilds)
# no hypre by default (as the Makefile); -H <hypre prefix> builds with hypre

export REEF3D_MPIRUN="mpirun"          # e.g. "mpirun --oversubscribe" if fewer cores than ranks

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
| `cfd_2d_nwt_komega` / `_kepsilon` | 1 | RANS models with free surface |
| `cfd_2d_nwt_porous_komega` | 1 | VRANS porous box + k-ω porous sources (B295) |
| `cfd_2d_nug_breaking_komega` | 1 | stretched grid, cnoidal waves, k-ω T36=2, slope solid |
| `cfd_2d_dambreak` (+ `_fcc3`, `_fcls3`, `_mpi2`) | 1/2 | closed tank, walls, N40=33 / 4, 2D MPI |
| `cfd_2d_dambreak_plic` (+ `_rk3`, `_fcc3`) | 1 | PLIC VOF (F80=4) with N40=3 / 13 / 33 |
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
| `nhflow_3d_6dof_box` | 2 | NHFLOW 3D box, all six DOFs free, initial roll/yaw |
| `nhflow_3d_shipwave_box` | 1 | NHFLOW moving pressure patch, ship-wave mode (X10=3, X400=2) |
| `sflow_shipwave_box` (+ `sflow_6dof_box_oneway`) | 1 | SFLOW ship pressure patch (X10=3) / one-way direct forcing (X10=2) |
| `fnpf_2d_6dof_box` (+ `_rk4`) | 1 | FNPF resolved floating body, A310=3 / 4, added-mass coupling |
| `nhflow_3d_ship_box_thrust` | 2 | ship module (X350, ship.dat): constant thrust with offset, ITTC-1957 friction + form factor |
| `nhflow_3d_ship_box_6dof` | 2 | ship module, all six DOFs: cross-flow drag (strips), roll damping |
| `cfd_2d_fem_obstacle` (+ `_mpi2`) | 1/2 | FEM solid (Z 30, N10=1): elastic obstacle hit by the bore, coupling across a subdomain border |
| `cfd_2d_fem_wall_failure` | 1 | FEM concrete wall cracking, erosion, debris, ground and part contact |
| `cfd_3d_fem_column` | 4 | FEM solid in 3D on 4 ranks: elastic column hit by the bore |
| `cfd_2d_fem_simple_wall` | 1 | FEM simple input (concrete C30 preset, fix base, monitor auto, resolution), settling with the initial water, hybrid loads, structural damping |
| `nhflow_2d_nwt_stokes5` (+ `_mpi2`) | 1/2 | NHFLOW relaxation generation + beach (B98=2/B99=1), Stokes 5th |
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
| `nhflow_3d_amr_heave_mpi2` | 2 | NHFLOW AMR heave decay of a floating box (X 10 1), body zone A 278 across the partition edge: forcing on every grid, loads from the finest |
| `nhflow_3d_amr_tow_zone` (+ `_vref_rk3`) | 1 | NHFLOW AMR towed box (X 10 2) with moving zone and wake wedge (A 278/279), regrid every step: from_old, prolongation, patch deletion; `_vref_rk3`: vertical refinement A 281, RK3, A 520 1 |
| `nhflow_3d_amr_adaptive` | 1 | NHFLOW AMR solution-adaptive (A 273, A 282), ring wave from an F 72 hump, regrid every step: tagging, F 72 boxes on fresh patches |
| `nhflow_3d_amr_wetdry` | 1 | NHFLOW AMR wetting and drying in the patches (A 283): beach and cone island, a box crossing the island shoreline, run-up at t = 0 (wet flags of the cells around a patch, lexer::wetfix, covered-cell flags) |
| `nhflow_3d_amr_wetdry_adaptive` (+ `_mpi2`) | 1/2 | the same with the shoreline flag (A 284) only, regrid every step, layout changes from step 116; `_mpi2`: partition edge through the island |
| `nhflow_3d_amr_breaking` | 1 | NHFLOW AMR wave breaking (A 550 1, A 512 2): a bore runs into a static patch; patch implicit diffusion with a patch-local BiCGStab, filled cells take the source grid's breaking (lexer::amrvb) |
| `nhflow_3d_amr_breaking_adaptive` (+ `_mpi2`) | 1/2 | breaking with the adaptive patch following the bore (A 273, A 285), wetting and drying, regrid every step; `_mpi2`: different patch counts per rank (local reductions, gcparaxijk_single) |
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
| `sflow_2d_custom_zones` | 1 | old input: SFLOW with custom B 107 / B 108 zones |
| `sflow_2d_dirichlet` | 1 | old input: SFLOW Dirichlet generation |
| `fnpf_2d_hdc` (stage `fnpf_hdc_source`) | 1 | old input: hydrodynamic coupling FNPF -> FNPF (FNPF state P 44 3, DIVEMesh H 10 44, B 92 61), as in FNPF nesting |
| `fnpf_2d_hdc_single` (stage `fnpf_hdc_source_single`) | 1 | the same with one state file per step (P 45 1); same result as `fnpf_2d_hdc` |
| `nhflow_2d_hdc` (+ `_mpi2`, stage `fnpf_hdc_source_uvw`) | 1/2 | FNPF -> NHFLOW nesting: FNPF state with velocities (P 44 1), DIVEMesh H 10 4, NHFLOW reads eta, u, w (B 92 61) |
| `nhflow_2d_hdc_single` (stage `fnpf_hdc_source_uvw_single`) | 1 | the same with one state file per step (P 45 1); same result as `nhflow_2d_hdc` |
| `nhflow_2d_hdc_nhflow` (stage `nhflow_hdc_source`) | 1 | NHFLOW -> NHFLOW nesting (DIVEMesh H 10 5); target time step close to the source state interval |

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
  the physics is right. Validation cases (long runs vs. measurements) are a separate set.
- A run fails on NaN/Inf in the per-step norms, a non-zero exit code, a timeout or a missing dump.

## Unit tests

`unit/` holds standalone tests of solver-independent kernels (no MPI, no REEF3D binary); the build
line is at the top of each file, run from `unit/`:

| test | covers |
|---|---|
| `rigidbody_test.cpp` | 6DOF rigid-body core: quaternion/Euler, constant force, oscillator and torque-free top (orders of RK2/RK3/RKLS3/RK4), DOF modes, damping |
| `ship_test.cpp` | ship module kernels: waterline clipping, wetted surface, draft strips, cross-flow drag, ITTC-1957 line, roll damping |

