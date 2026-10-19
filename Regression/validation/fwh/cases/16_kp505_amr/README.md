# 16 KP505 open water with CFD mesh refinement (A 10 6, G 1 2)

The KP505 model propeller (D = 0.25 m, 5 blades, `KP505_model_250mm_5blades_1.5mm.stl` as
`floating.stl`, not in the repository) at n = 10 1/s (X 211 -62.83 rad/s about x, prescribed),
U = 1.5 m/s (W 10 over the 1 m² inlet), J = 0.6, BPF = 50 Hz; k-omega (T 10 2), single phase.
Experiment at J = 0.6: KT = 0.219, 10 KQ = 0.357.

Why: the run on a uniform 1 cm grid (246 k cells, user, 2026-10-16) reproduced the CFD pressure with
FW-H at the near observers (BPF within 3-6 %, one point 19 %) and the box dipole in the far field
(+11 % at 10 m), but the propeller itself was not resolved: blades thinner than a cell, KT = 0.016,
10 KQ = 0.122, a spurious side force at BPF (0.9 N) and a non-axisymmetric field (34 against 17 Pa
at 2R on the y and z axes). The far-field BPF of that run is numerical: in a uniform inflow the
total force of a symmetric rotor has no BPF component and the steady-loading (Gutin) noise is
negligible at a tip Mach number of 0.005.

Grid: base dx = 8 mm, uniform in x in [-0.12, 0.40], |y|, |z| < 0.22, growing to 0.2 m; two refined
levels (G 1 2), 4 mm and 2 mm (about 40 cells along the chord), in two static boxes (G 15):
the rotor x in [-0.05, 0.045], |y|, |z| < 0.14 (the rotor: x in [-0.030, 0.023], r <= 0.125) and the
wake x in [0.045, 0.30] (1 D), |y|, |z| < 0.13 (the contracting slipstream and the tip vortices).
Every G 15 box is refined to the finest level (reefamr3d: the coarser levels only nest around it
with 3 cells), so the wake box is at 2 mm too. Estimate: level 2 0.94 M (rotor) + 2.2 M (wake),
level 1 ~0.5 M, base ~0.3 M, ~3.9 M cells. A longer wake costs 0.85 M cells per 0.1 m at 2 mm;
1 mm (G 1 3) is 8 times that. All levels take the time step of the finest one (G 7 0, CFL N 47 0.3; tip speed
7.9 m/s at 2 mm: dt ~ 7.5e-5 s, ~11 k steps for 0.8 s = 8 revolutions).

CFD AMR needs (cfd_amr::scope): N 40 2 or 3 (here 3; the uniform run used 4), no partition in z
(the DIVEMesh default M 13 0), no patch boundaries (B 440-442), X 10 0 or 1. If a condition is not
met the run prints "CFD AMR (G 1 2): needs ...; running without mesh refinement" at the start -
check the log for it.

FW-H box U 20 x in [-0.09, 0.15], |y|, |z| < 0.18, in the uniform part of the base grid: 4-5 base
cells around the rotor box; its downstream face lies in the refined wake. The acoustics reads the
level-0 fields, which under a patch hold the values restricted from it.

Output: LAGOON store (P 18 1, instead of the VTU files) `REEF3D_CFD.lagoon` every 0.01 s (P 30):
the level-0 volume and the propeller (body mesh and motion). The refined patches are not in the
store (cfd_amr_print writes them as .vtr per patch in REEF3D_CFD_AMR/ with a .vtm per output; the
LAGOON AMR writer is used by FNPF, NHFLOW and SFLOW only). P 18 2 writes the VTU files as well.

Observers 1-7 =
probes 1-7: the near points of the uniform run plus two more at 2R in the rotor plane on -y and -z
(axisymmetry: observers 1, 3, 6, 7 should give the same BPF amplitude); 8, 9 at 10 m.

Running: copy `floating.stl` (the KP505 STL), set M 10 in control.txt and ctrl.txt to the number
of ranks, `DiveMESH`, `mpirun -np <ranks> REEF3D`. Not run here (prepared 2026-10-16).

Analysis (in the run directory):

    python3 ../../../scripts/kp505_analysis.py 0.4 0.8

KT, KQ from the 6DOF forces and the box momentum balance against the open-water curve, and the BPF
amplitude and phase of FW-H and the probes at each point.
