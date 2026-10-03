# REEF3D regression test suite

Small, fast cases that exercise the main code paths, plus tools to

- **A/B-compare two binaries bitwise** (for refactors and clean-up: "did anything change?"), and
- **check against stored references** with a tolerance (long-term regression testing).

Pure Python 3 (standard library only) + bash. CFD, NHFLOW and FNPF cases; the layout is
solver-independent (`"solver"` in `case.json`), SFLOW cases can be added the same way once its
driver loop calls `regression_dump`.

## How it works

`src/regression_dump.{h,cpp}` is a small output class called from the CFD, NHFLOW and FNPF driver loops. It is
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
cd tests/regression

# build two binaries with identical, reproducible flags (-O2 -ffp-contract=off, no -march=native/LTO)
./build_reef3d.sh -r hans_dev -o ~/reef3d_reg/ref      # a git revision (exported, working tree untouched)
./build_reef3d.sh             -o ~/reef3d_reg/new      # the working tree (incremental rebuilds)

export REEF3D_MPIRUN="mpirun"          # e.g. "mpirun --oversubscribe" if fewer cores than ranks

# run both and compare (default: must be bitwise identical)
./regression.py ab --ref-bin ~/reef3d_reg/ref/REEF3D --new-bin ~/reef3d_reg/new/REEF3D \
                   --divemesh ~/Codelite/DIVEMesh/bin/DiveMESH --out ~/reef3d_reg/ab --tags quick
```

The report is printed and written to `<out>/new/compare.md` (+ `compare.json`).

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
| `nhflow_2d_nwt_stokes5` (+ `_mpi2`) | 1/2 | NHFLOW relaxation generation + beach (B98=2/B99=1), Stokes 5th |
| `nhflow_2d_dirichlet` | 1 | NHFLOW Dirichlet wave generation (B98=3) |
| `nhflow_2d_awa` | 1 | NHFLOW active wave generation + active absorption (B98=4/B99=3) |
| `nhflow_2d_current` | 1 | NHFLOW waves on a current: discharge inflow (B60=1, W10), beach relaxes to the current (B97=1) |
| `nhflow_2d_irregular_decomp` | 1 | NHFLOW JONSWAP waves, decomposed relaxation precalc (B89=1) |
| `nhflow_3d_irregular_decomp` | 2 | NHFLOW 3D short-crested waves (B130=2), decomposed precalc, CFL 0.5 |
| `fnpf_2d_regular` (+ `_mpi2`) | 1/2 | FNPF relaxation generation + beach, regular waves |
| `fnpf_2d_irregular_decomp` | 1 | FNPF JONSWAP waves, decomposed relaxation precalc (B89=1) |
| `fnpf_3d_shortcrested` | 2 | FNPF 3D short-crested waves (B130=2) |

Tag `quick` selects a subset that runs in a few minutes. Adding a case: copy a directory, edit,
run `./regression.py run ... --cases <new>`, check it, then `bless`.

## Notes

- Irregular-wave cases must fix the random seeds (`B 139`, and `B 138` for directional
  spreading); otherwise the phases come from `srand(time(0))` and no two runs agree.

- Runs use few steps on coarse grids: they test that code paths give the same numbers, not that
  the physics is right. Validation cases (long runs vs. measurements) are a separate set.
- A run fails on NaN/Inf in the per-step norms, a non-zero exit code, a timeout or a missing dump.
