# REEF3D validation suite

Cases with an analytical (or reference) solution: they check that the numbers are *right*, while
`tests/regression` checks that they stay the *same*. Same case format and runner as the regression
suite (`case.json` + `control.txt`/`ctrl.txt`, or `"base"` + key overrides), plus

| key | meaning |
|---|---|
| `"time"` | simulated time (sets `N 41`; instead of a number of steps) |
| `"check"` | `{"type": <checker>, ...}` – how to evaluate the run |
| `"check_set"` | override single check parameters of the base case |

```bash
cd tests/validation
./validation.py list
./validation.py run --reef3d BIN --divemesh DM --out DIR --check [--cases 'channel_*'] [--tags quick]
./validation.py check DIR            # evaluate an existing run again
```

The report (`DIR/validation.md`) lists PASS/FAIL per case with the error measure and the
tolerance, plus the probe-by-probe details.

## Checkers

| type | what | error measure |
|---|---|---|
| `conserved_integral` | volume integral of a dumped field (Σ field·vol, needs the `vol` field of the state dump) | \|I_end − I_start\| / \|I_start\| |
| `max_abs` | maximum \|field\| in the final state, e.g. `div` (divergence after the projection) | max \|field\| |
| `channel_startup` | body-force driven laminar channel flow from rest between no-slip walls (exact Fourier-series solution); probes `P 61` | max \|u − u_exact\| / u_max at the last time (`tol_final`) and the earlier times (`tol_transient`) |

## Cases

| case | what |
|---|---|
| `channel_2d_fc3` (+ `_rk3`, `_rk2`, `_fc2`, `_fcls3`, `_rkls3`, `_rkls3_sf`, `_fcc3`) | all momentum time schemes with implicit diffusion (D 20 2) |
| `channel_2d_fc3_plic`, `_rk3_plic`, `_fcc3_plic` | the PLIC (F 80 4) momentum classes |
| `channel_2d_fc3_explicit` | explicit diffusion (D 20 1), reference for the time integration |
| `channel_2d_fc3_cfl01`, `_cfl003` | time-step convergence (CFL 0.3 / 0.1 / 0.03) |
| `closed_dambreak_2d_div`, `closed_dambreak_3d_obstacle_div` | dam break in a closed tank (2D; 3D with obstacle, 4 ranks): velocity divergence-free in all cells after the projection, also next to walls and solids |
| `conc_diffusion_2d_uniform`, `_stretched` | passive concentration blob diffusing in still water (closed tank, 2D), uniform and stretched grid: total amount conserved |

Set-up notes: the channel cases use `B 20 3` (second-order Dirichlet wall ghost cells) and
`D 22 2` (implicit diffusion uses the wall ghost values). With the defaults `B 20 2` (ghost = 0, the
wall sits half a cell outside) and `D 22 1` (zero-gradient at walls in the implicit diffusion, i.e.
no wall friction without wall functions) the laminar channel is not reproduced.
