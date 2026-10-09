# 14 KP505 resolved open water on Betzy (4 nodes, 512 ranks)

Fine version of case 13: KP505 model scale (D = 0.25 m, Z = 5), n = 10 rps, J = 0.6 (U = 1.5 m/s,
`W 10 1.5` over the 1 m² inlet), Re(0.7R) ≈ 5e5, LES WALE (`T 10 33`), 1 mm cells over the rotor
region (blade 0.7R ≈ 4 cells thick), growing by 1.08 to 20 mm; 237×372×372 = 32.8 M cells
(~64 k per rank). Fixed Δt = 3.5e-5 s (CFL 0.27 at the tip, 2857 steps per revolution), ramp over
the first revolution (`X 205 2`, `X 206 0 0.1`), 8 revolutions (`N 41 0.8`, ~23 k steps).
FW-H box [−0.05, 0.08]×[−0.15, 0.15]² (25 cells off the tips) with the mean flow `U 50 1.5 0 0`,
observers in the rotor plane at 1R and 1.5R off the tip (r = 2R, 3R), off the plane, upstream on
the axis, and at 40 D (lateral and downstream). Experiment at J = 0.6: KT = 0.219, 10 KQ = 0.357.

## Files to add on Betzy

- `floating.stl`: `KP505_model_250mm_5blades_1.5mm.stl` (made by `scripts/convert.py`, 30 MB).
- `REEF3D`: built from `ahmet_dev` (`make -j 16` with the same toolchain module as in the job).
- `DiveMESH`: built from the current DIVEMesh (grid format v2).
- In `job_*.slurm`: the project account (`--account=nnXXXXk`) and, if needed, the toolchain module.

## Steps

1. `sbatch job_test.slurm` — 300 steps in `./test`: grid, `FW-H:` line, time per step, no
   `EMERGENCY`. The blades first sweep cells during the ramp; case 13 diverged there with the
   default `X 41 0.6` at 4 mm. `X 41 1.0` is set here; if the test diverges use 1.5.
2. Time per step × 23 k steps gives the wall time; set `--time` of `job_run.slurm` (default 12 h).
3. `sbatch job_run.slurm`.
4. `python3 ../../scripts/kp505_analysis.py 0.5 0.8` in the case directory (last 3 revolutions).

Rough estimate from the 8-rank runs here (~10 µs per cell and step on one core): 0.6–1 s per step
on 512 ranks, 4–7 h for 23 k steps.

## Restart

State files every 0.1 s (`P 40 1`, `P 42 0.1`). The state holds the flow fields and the time, not the
6DOF orientation: give the propeller angle at the restart time with `X 101`. With n = 10 rps every
multiple of 0.1 s after the ramp is a full revolution, and the angle is

    phi = omega * 0.1 s * (1/2 - 2/pi^2) = -107.05 deg   (ramp X 205 2 over 0.1 s)

so for a restart from any state at t = 0.1, 0.2, ... add to `ctrl.txt`

    I 40 1
    I 41 0
    X 101 -107.05 0.0 0.0

(`P 40 1` overwrites state 0 each time, so `I 41 0` is the last state; `X 206` can stay: the
ramp is 1 after 0.1 s). The state can be up to 0.1 s older than the last FW-H samples: the
samples after the restart overlap those, keep the later ones (the analysis scripts use the last
value of each time). The FW-H files are appended to after a restart.
Check the angle against `REEF3D_CFD_6DOF/REEF3D_6DOF_position_0.dat` of the first run (in case 13,
n = 5, the formula gives −69.7° at t = 0.509 s, the output −70.2°).
