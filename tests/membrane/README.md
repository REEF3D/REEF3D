# Membrane (X 330) test cases, REEF3D::NHFLOW

Run DIVEMesh, then REEF3D with 2 processes (`M 10 2`). The Poisson solver is reefmg (`N 10 1`).
Time series: `REEF3D_NHFLOW_Membrane/REEF3D_NHFLOW_Membrane_0.dat`; membrane geometry and loads:
`REEF3D_NHFLOW_Membrane_VTP/REEF3D-NHFLOW-Membrane-0-*.vtp` with the collection `REEF3D-NHFLOW-Membrane-0.pvd`
(open the .pvd in ParaView for the time series).

## Fixed membrane (static bag)

Closed bag at rest with the inner water level raised by `fill` (dh = 0.05 m) above the outside level.
The membrane is modelled as a porous jump with hydraulic resistance R_n = 1e4 m/s, so the only flow
through it should be the Darcy leakage u_n = g dh / R_n, and the floor should carry the hydrostatic
load -rho g dh A_floor.

    python3 ../../analyse.py REEF3D_NHFLOW_Membrane/REEF3D_NHFLOW_Membrane_0.dat 1e4 <A_exposed>

Both cases run with `A 520 1`; `A 520 2` works as well. The default static floor pressure of a fixed
membrane is `floorpressure 3` (shape of the discrete free-surface ramp, running average over tau = 2 s,
frozen at t = 2 tau).

* `nhflow_bag2d_static`: 2D tank 20 m x 5 m depth, bag x = 8 .. 12 m, floor at z = 2 m, dx = 0.1 m, 20 sigma layers.
  A_exposed = (2 walls x 3 m + 4 m floor) x 0.1 m cell width = 1.0 m^2.
  Reference (30 s): flux through the membrane / Darcy 0.97 (A 520 1) and 0.96 (A 520 2), floor load 0.91 of
  -rho g dh A (the rest ends up on the wall panels at the corners), max|U| 9e-4 m/s (A 520 1) and 1.5e-3 m/s
  (A 520 2), water volume conserved to round-off. For this axis-aligned 2D box `floorpressure 0` is quieter
  (max|U| 3e-4 m/s).
* `nhflow_bag3d_static`: cylinder bag R = 5 m, floor at z = 2.5 m, depth 6 m, dx = 0.4 m, 12 sigma layers,
  `projections 2`. A_exposed = 78.5 m^2 floor + 110 m^2 wall.
  Reference (20 s): flux / Darcy 0.94, floor load 0.89, max|U| 6.8e-3 m/s (A 520 1) and 1.2e-2 m/s bounded
  (A 520 2; with floorpressure 0 it is 2.8e-2 m/s and growing). With `projections 1` the leak rate is
  1.65 x Darcy at dx = 0.4 m (converges with the grid).

## Moving membrane

* `nhflow_bag2d_flexible`: the 2D bag of `nhflow_bag2d_static` as a flexible membrane (1 kg/m^2, E t = 5e5 N/m),
  top edge held in place, inner level +5 cm at the start. The floor sags by 2 - 5 cm and the bag "breathes":
  the head difference oscillates with a period of ~8.6 s and decays (30 s, max tension ~500 N/m, volume
  conserved). The staggered coupling (default) adds an inertia of about rho R_n dt per area, so this oscillation
  is slower than the physical one and depends on the time step; see `nhflow_bag2d_flexible_iterated`.
* `nhflow_bag2d_towed`: flexible bag (4 m wide, floor 3 m below the surface) hanging from a collar of two pontoons,
  collar towed at 0.2 m/s in still water (prescribed motion: `X 10 2`, `X 11 2 0 0 0 0 0`, `X 210 0.2 0 0`,
  ramp `X 205 2`, `X 206 0 3`). Without a sinker weight the front wall is pushed in, the floor lifts by up
  to ~0.6 m and the inner level rises by up to ~0.2 m (30 s, stable).
* `nhflow_bag2d_collar_rigid`: rigid bag carried by the freely floating collar (`X 10 1`, pitch locked) in regular
  waves H = 0.2 m, T = 4 s. Stable over 20 s with the added-mass stabilisation of the body coupling
  (`bodyaddedmass`, default 2 rho V_bag); the bag suppresses the surge drift of the collar.

## Strong coupling of a flexible membrane (`coupling iterated`)

The projection of every Runge-Kutta stage is repeated until the membrane velocities the fluid used and the ones
the structure returns agree (IQN-ILS, see `src/net_membrane_coupling.cpp`). The time series gets two more columns:
iterations per time step (all stages) and the largest relative residual at the end of a stage.

Defaults of the iterated coupling: membrane mesh h = 1.5 delta (elements as large as the coupling layer; the layer
does not see shorter structural modes, which are then nearly free and stall the iteration: 3D > 100 iterations per
step with h = max cell size), Robin preconditioner scaled by 16 (`couplingrobin`; the tolerances are divided by the
same factor, so the converged solution does not depend on it). Each iteration restores the predicted velocity and
the pressure of the stage start, so the fluid response is a deterministic function of the membrane velocity.

* `nhflow_bag2d_flexible_iterated`: `nhflow_bag2d_flexible` with `coupling iterated`. The bag breathes with the
  real inertia of the bag water: zero crossings of the head difference at 0.74, 2.29 and 4.94 s; on the finer
  membrane mesh h = 0.25 m 0.74, 2.27, 4.94 s, identical for CFL 0.3 and 0.15 (`N 47`). The staggered coupling
  depends on the time step (first zero crossing at 2.14 s for CFL 0.3 and 1.63 s for CFL 0.15). 4.1 iterations
  per time step (RK2, 2 per stage).
* `nhflow_bag3d_flexible_iterated`: the cylinder bag of `nhflow_bag3d_static` as a flexible membrane, top edge held,
  +5 cm (253 nodes). First zero crossing of the head difference at 2.70 s (2.64 s at CFL 0.15), staggered 5.61 s;
  floor sags to z = 2.20 m (staggered 2.43 m). 5.4 iterations per time step (4.2 at CFL 0.15), all stages
  converged; cheaper per step than the staggered coupling on its default (finer) membrane mesh.
* `nhflow_bag2d_collar_flexible`: flexible bag on a freely floating collar (`X 10 1`, pitch locked) in regular
  waves H = 0.2 m, T = 4 s. Pontoons 1.2 m high (half submerged), bag top 0.4 m above the still water level.
  Stable over 30 s, all stages converged, 4.2 iterations per time step (95 % below 6), tension below 1640 N/m,
  collar heave -0.13 .. +0.10 m, mean surge +0.2 .. +0.5 m without mooring (no drift away).
  With the 0.6 m pontoons of `nhflow_bag2d_collar_rigid` (reserve buoyancy 471 N per 0.1 m slice, bag top
  0.2 m above the still water level) the bag pulls the collar under after about 3 wave periods, identically at
  half the time step (not a coupling instability). The staggered coupling diverges at ~3 s in both set-ups.
  On the finer membrane mesh h = max cell size (0.25 m) the unmoored collar drifted downwave (2.8 m in 30 s; the same
  for dt / 2 and dx / 2, i.e. converged in the fluid discretisation, but not in the membrane mesh): the structural
  modes shorter than the coupling layer are free, the wall then hangs from the collar like a hinge and restrains
  it horizontally only once its top segment tilts. With h = 1.5 delta the drift is gone.

Settings: `coupling iterated [rtol [n]]` (default 1e-3, 50), `couplingtol rtol atol` (atol default 1e-5 m/s rms),
`couplingrobin f` (16), `couplingreuse n` (8), `couplingrelax w` (0.5), `couplingcolumns n` (100),
`couplingfilter eps` (1e-2), `couplingqn ils|imvj` (ils; imvj carries the inverse Jacobian over, n x n per stage,
no gain in these tests), `couplinglog 1` (residual of every iteration). With the coupling off (default
`coupling staggered`) all cases reproduce hans_dev 6769a1e bit for bit.
