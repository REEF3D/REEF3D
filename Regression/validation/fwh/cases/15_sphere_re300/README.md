# 15 Sphere at Re = 300: vortex-shedding dipole, FW-H against Curle, near to far field

A fixed sphere (D = 1 m, 6DOF direct forcing, `Regression/cases/cfd_3d_heave_sphere_6dof/floating.stl`)
in a uniform flow U = 1 m/s (`B 60 1`, `W 10 100` over the 10 m × 10 m inlet), ν = 1/300 m²/s, single
phase, laminar. At Re = 300 the wake sheds hairpin vortices periodically in one plane; the lateral
force oscillates at St ≈ 0.137 and radiates as a dipole. One lateral kick of 0.05 D in the first
4 s (`X 11 0 2 0 0 0 0`, `6DOF_motion.dat` from `scripts/make_inputs.py shedding 0.05`) fixes the
plane of the shedding (x–y) and shortens the transient.

References
- Hydrodynamics: Johnson & Patel (1999), Re = 300: Cd 0.656, mean Cl 0.069, St 0.137.
- Acoustics: Curle's compact dipole of the force on the sphere, with the force from the 6DOF output
  (surface integral) and from the momentum balance of the FW-H box (the dipole FW-H radiates). In
  water the wavelength at f = 0.137 Hz is ~11 km: the observers go from the hydrodynamic near
  field (p ~ 1/r²) through r = λ/2π ≈ 1.7 km to the acoustic far field (p ~ 1/r, delay r/c).

Grid (`B 127-129` cell-based): uniform dx around the sphere and the near wake
(x ∈ [−1, 5], |y|, |z| < 1.25 for D/20 and D/40), growing to 0.25 m; domain
[−5, 15]×[−5, 5]².

| run | dx | cells | ranks | where | end time |
|---|---|---|---|---|---|
| d10 | D/10 | 0.25 M | 8 | local, ~1 h | 150 s |
| d20 | D/20 | 2.0 M | 128 | Betzy | 200 s |
| d40 | D/40 | 12.5 M | 512 | Betzy, VTU frames 180–200 s every 0.33 s | 200 s |

`make_case.py` writes `d<N>/control.txt`, `ctrl.txt`, `6DOF_motion.dat` and `job.slurm`. FW-H box
[−1.5, 2.5]×[−1.5, 1.5]² with the mean flow (`U 50 1 0 0`). Observers: ring r = 5 m (x–y plane),
y = 3 m … 100 km, z = 10 m and 20 km; for d20/d40 also a ring at 20 km and two observer planes z = 0
for the animation, ±40 m (5 m spacing) and ±40 km (5 km spacing), sampled every 0.05 s (`U 31`).

## Running on Betzy

In each run directory: `control.txt`, `ctrl.txt`, `6DOF_motion.dat`, `job.slurm` from `d20/` or
`d40/`, `floating.stl` (`Regression/cases/cfd_3d_heave_sphere_6dof/floating.stl`), and `REEF3D`,
`DiveMESH` built on Betzy; `job.slurm` stops with a message if one of them is missing. Set the account in `job.slurm` (and the toolchain module if it differs), then
`sbatch job.slurm`. d20 needs one node (128 ranks); if the partition asks for more nodes, use
`python3 make_case.py 20 512 200`.

## Analysis

From a run directory inside `d20/` or `d40/` (otherwise give the path to `Regression/validation/fwh/scripts`):

    python3 ../../../scripts/sphere_analysis.py 120 200      # Cd, Cl, St; FW-H vs Curle per observer
    python3 ../../../scripts/fwh_planes_vtk.py 180 200 0.33  # FW-H planes -> REEF3D_FWH_Planes/fwh_L5_*.vtk, fwh_L5000_*.vtk
    python ../../../scripts/sphere_plots.py --noshow          # sphere_fig1..4.png (forces, |p| vs r, directivity, time series)
    python ../../../scripts/sphere_slide_figs.py slides 89 177            # the same in the style of the slides
    python ../../../scripts/sphere_fwh_animation.py slides 150 3 --label "sphere Re = 300, D/20"  # slides/fwh_planes.mp4

The window has to start after the transient and after the delay r/c of the farthest observer
(100 km: 67 s), so 120–200 s. Expect FW-H = Curle (box force) from r ≈ 30 m on; on the ring r = 5 m
(2.5 m from the box, downstream points in the wake) the near field is not compact and the wake
contributes, so FW-H and Curle differ there. `d10/results_partial.txt`: the local D/10 run stopped
at t = 94.7 s (window 50–94 s): FW-H/Curle(box) = 1.10, 1.02, 0.996, 0.994, 0.994 at r = 30 m,
100 m, 1 km, 5 km, 20 km (phase ≤ 5°); Cd 0.59 (surface) / 0.80 (box), St 0.124 at D/10.

## Results

| | Johnson & Patel | D/10 (local, partial) | D/20 (`d20/results.txt`) |
|---|---|---|---|
| St | 0.137 | 0.125 | 0.1358 |
| Cd, box momentum balance | 0.656 | 0.80 | 0.670 |
| Cd, 6DOF surface integral | 0.656 | 0.59 | 0.608 |
| mean Cl, box | 0.069 | 0.060 | 0.0685 |
| FW-H / Curle (box) at r = 100 m | | 1.02 | 1.020 |
| FW-H / Curle (box) at r = 1–100 km | | 1.00 | 1.002 |

- The FW-H integral reproduces the compact dipole of the box force within 0.2 % from 1 km to
  100 km, through the change from 1/r² to 1/r at λ/2π ≈ 1.7 km and with delays up to 67 s; the
  directivity at 20 km is the dipole figure eight of Curle.
- Closer than ~50 m FW-H is larger than Curle (×1.07 at 40 m, ×4 at 3 m, largest downstream on the
  ring r = 5 m): the sphere and the box are not compact there and the wake pressure contributes;
  the same ratios at D/10 and D/20, so not a grid error.
- The direct-forcing surface force converges slowly: Cd 0.59/0.80 (surface/box) at D/10,
  0.61/0.67 at D/20; the box momentum balance is the better force at these resolutions.

## Animation (ParaView)

- Flow: the VTU frames of d40 (`REEF3D_CFD_VTU`), Gradient filter with Q-criterion, contour of Q
  coloured by the streamwise velocity (hairpin vortices), a slice z = 0 of the pressure, the
  sphere from `REEF3D_CFD_6DOF_VTP`.
- Sound: `REEF3D_FWH_Planes/fwh_L5_*.vtk` (near field, p' r² shows the oscillating dipole lobes) and
  `fwh_L5000_*.vtk` (acoustic field, p' r shows the wave fronts leaving at 1500 m/s, λ ≈ 11 km),
  time-synchronised with the flow frames (same t0, dt).
