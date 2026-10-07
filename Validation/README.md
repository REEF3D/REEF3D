# REEF3D validation cases

Cases that check a model change against physics: analytical or equilibrium solutions, balances,
serial against MPI, before against after a patch. Each folder has the inputs, the result plots and a
README with the set-up, a results table, the builds that were compared and how to run it.

| set | question | runs |
|---|---|---|
| `Regression/` | does the code give the **same** numbers as before? | short runs, stored references, automatic |
| `Benchmark/` | does REEF3D reproduce the **published benchmarks**? | experiments and exact solutions, PASS/FAIL |
| **`Validation/` (this)** | is a **model change** right? | cases and evaluation scripts written for one review, run by hand |

## Layout

```
Validation/<topic>/tools/          shared run and evaluation scripts of a topic
Validation/<topic>/NN_<name>/      README.md, cases/<case>/{control.txt,ctrl.txt}, results/
```

Cases run as-is: copy `control.txt` and `ctrl.txt` into a run folder, run DiveMESH and REEF3D
(`<topic>/tools/run_case.sh` does this). The scripts need Python 3 with numpy (matplotlib for plots).

## Contents

| folder | what it checks |
|---|---|
| `turbulence/01_channels_mpi_freesurface` | SFLOW, NHFLOW and CFD k-ε / k-ω in open channels against the equilibrium (Rastogi–Rodi, log law), MPI against serial, free-surface damping T36 |
| `turbulence/02_nhflow_momentum_diffusion_rk` | NHFLOW implicit momentum diffusion with the RK stage weight: momentum balance (ν+ν_t)du/dz = gS(h−z), RK2 against RK3 |
| `turbulence/03_turbulence_options` | NHFLOW side-wall turbulence vs momentum walls, buoyancy term T45, vegetation sources, LES T21 2, SFLOW defaults (A212, A264) |
| `turbulence/04_sflow_shoreline_diffusion` | SFLOW horizontal momentum diffusion at the shoreline: free slip at dry neighbours (A212 1/2), sloping bank against the local Manning velocity |
| `turbulence/05_stretched_grids_diffusion_wetdry` | CFD `ifou` on stretched grids (k-ε, T 12 1), NHFLOW explicit vs implicit momentum diffusion (A 512 1/2) on a bump, NHFLOW k-ε at a wet/dry shoreline |
| `turbulence/06_inflow_wallshear_wetting` | CFD discharge inflow profile and wall shear per RK stage (k-ε channel), NHFLOW k-ε in columns wetted during the run (bank against local equilibrium) and a bore on a dry beach |
| `analytic/` | CFD against analytical solutions, automatic PASS/FAIL (`validation.py`, Regression runner): laminar channel start-up for all momentum schemes (also PLIC, single phase, 2 MPI ranks with periodic partition faces), time-step convergence, divergence after the projection, conservation of a diffusing scalar, NWT wave height, log law in the turbulent channel |
