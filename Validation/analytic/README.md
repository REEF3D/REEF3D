# REEF3D analytical validation (CFD)

Cases with an analytical (or reference) solution: they check that the numbers are *right*, while
`Regression/` checks that they stay the *same*. Same case format and runner as the regression
suite (`case.json` + `control.txt`/`ctrl.txt`, or `"base"` + key overrides), plus

| key | meaning |
|---|---|
| `"time"` | simulated time (sets `N 41`; instead of a number of steps) |
| `"check"` | `{"type": <checker>, ...}` – how to evaluate the run |
| `"check_set"` | override single check parameters of the base case |

```bash
cd Validation/analytic
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
| `log_law` | fully developed rough-wall turbulent flow: velocity at the probes `P 61` against u = (u*/κ) ln(30 z/ks) (or ln(z/z0) with `z0`) at the last time | max \|u − u_log\| / u_log |
| `wave_height` | free-surface elevation at the wave gauges `P 51` (max − min in `t_start`..`t_end`) against the target height `H` | max \|H_gauge − H\| / H |

## Cases

| case | what |
|---|---|
| `channel_2d_fc3` (+ `_rk3`, `_rk2`, `_fc2`, `_fcls3`, `_rkls3`, `_rkls3_sf`, `_fcc3`) | all momentum time schemes with implicit diffusion (D 20 2), second order in time (D 23 2): errors 2.5e-3 = the spatial error of 20 cells |
| `channel_2d_fc3_mpi2` | as `channel_2d_fc3` on 2 ranks with 16 cells in x: the periodic faces are partition boundaries (periodic halo exchange, periodic neighbours in the solvers) |
| `channel_2d_fc3_plic`, `_rk3_plic`, `_fcc3_plic` | the PLIC (F 80 4) momentum classes |
| `channel_2d_singlephase_fc3` (+ `_fc2`, `_fcls3`, `_fcc3`, `_rk3`) | the same without level set (F 30 0) |
| `channel_2d_fc3_explicit` | explicit diffusion (D 20 1), reference for the time integration |
| `channel_2d_fc3_cfl01`, `_cfl003` | time-step convergence (CFL 0.3 / 0.1 / 0.03) |
| `channel_2d_fc3_d23_1` | first-order weights of the implicit diffusion (D 23 1), for comparison |
| `closed_dambreak_2d_div`, `closed_dambreak_3d_obstacle_div` | dam break in a closed tank (2D; 3D with obstacle, 4 ranks): velocity divergence-free in all cells after the projection, also next to walls and solids |
| `inflow_outflow_2d_div` | single-phase channel with constant-discharge inflow and outflow: divergence-free after the projection, incl. the inflow and outflow columns (~1e-6) |
| `conc_diffusion_2d_uniform`, `_stretched` | passive concentration blob diffusing in still water (closed tank, 2D), uniform and stretched grid: total amount conserved |
| `conc_diffusion_2d_walls` | as `_uniform` with the blob in the corner on bed and wall: no-flux walls for the scalar diffusion (was +1.4e-5 per 100 steps with the lagged ghost value) |
| `turbchannel_2d_komega`, `turbchannel_2d_kepsilon` | rough-wall closed channel (walls at z = 0 and 2 m), periodic in x, k-ω / k-ε with wall functions, 400 s to steady state: log law within 5 % at z = 0.025..0.275 m (k-ω 0.8 %, k-ε 2.7 %); the periodic faces must not be treated as walls by the wall laws, the k-ω boundary matrix and the initialisation |
| `nwt_2d_stokes2_height` | 2D wave tank (regression case `cfd_2d_nwt` for 24 s), Stokes 2nd order H = 0.02 m: wave height at x = 6, 8, 10 m within 5 % (0.7 % now, before and after the inflow/wave-generation pressure change) |

Set-up notes: the channel cases use `B 20 3` (second-order Dirichlet wall ghost cells) and
`D 22 2` (implicit diffusion uses the wall ghost values). With the defaults `B 20 2` (ghost = 0, the
wall sits half a cell outside) and `D 22 1` (zero-gradient at walls in the implicit diffusion, i.e.
no wall friction without wall functions) the laminar channel is not reproduced.
The channel is periodic in x with 4 cells and uses `N 10 3` for the pressure: the default
hypre PFMG Poisson solver (N 10 14, N 11 11) fails for this periodic set-up once the rhs is not
exactly zero (with 6 cells and with 4).

Time accuracy of the implicit diffusion (channel, error against a run at CFL 0.01, same scheme):

| N40 | D 23 1, CFL 0.3 / 0.1 / 0.03 | D 23 2, CFL 0.3 / 0.1 / 0.03 |
|---|---|---|
| 3 (SSP-RK3) | 6.6e-3 / 2.0e-3 / 4.4e-4 | 1.2e-4 / 2.7e-5 / 1.1e-5 (solver tolerance) |
| 2 (SSP-RK2) | 7.9e-3 / 2.4e-3 / 5.3e-4 | 2.2e-3 / 2.8e-4 / 3.5e-5 |
| 4 (low-storage RK3) | 3.3e-3 / 9.9e-4 / 2.2e-4 | 6.0e-4 / 8.3e-5 / 8.8e-6 |
