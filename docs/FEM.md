# REEF3D FEM — deformable, failing and collapsing solid structures

Explicit finite element solver for solids, two-way coupled with REEF3D::CFD.
Covers elastic deformation and stresses, plasticity, concrete cracking and crushing,
element erosion, fragmentation and debris (contact with the ground and between fragments).
No external libraries (Eigen only).

## Activation (ctrl.txt)

| Flag | Type | Meaning | Default |
|---|---|---|---|
| `Z 30` | int | 0 off, 1 FEM solids coupled with CFD; input file `fem.dat` in the case directory | 0 |
| `Z 31` | double | print interval of `REEF3D_FEM/REEF3D-FEM-*.vtu` in seconds, 0 = no VTU (`print` in `fem.dat` overrides) | 0.0 |

Works with the CFD momentum schemes that call `momentum_forcing_start`
(FC2, FC3, FCLS3, FCC3, RK2, RK3, RK3CN, RKLS3, PLIC variants), like the rod trees.
Not hooked into `momentum_RKLS3_df/sf` and not into NHFLOW/SFLOW/FNPF.

Fluid density `W 1` and gravity `W 20/21/22` are taken from ctrl.txt (set `W 22 -9.81`).
In 2D runs (`j_dir 0`) the solid is computed in plane strain: give it one voxel layer in y.

## Method

**Mesh.** All solids live on one uniform voxel lattice (spacing `hx hy hz`). Boxes and closed
STL surfaces fill the voxels whose centres they contain, in input order (later shapes overwrite,
`remove` carves). Touching voxels share nodes, so touching bodies are glued; leave one empty voxel
between bodies that should be separate.

**Elements.** 8-node hexahedra in a Total-Lagrangian formulation (finite rotations, small to
moderate strains). `reduced` (default): one-point integration with Flanagan–Belytschko stiffness
hourglass control on the current positions (rotation invariant). `full`: 2×2×2 Gauss points,
no hourglass modes, somewhat stiff in bending (shear locking with few elements over the thickness).
Linear + quadratic bulk viscosity in compression damps impact ringing.

**Time integration.** Central differences with lumped mass. The solid subcycles inside every fluid
step with `dt_s <= cfl * h_min / c_p`; the number of substeps is printed in the log.

**Materials.**

* `elastic`: St. Venant–Kirchhoff.
* `plastic`: J2 plasticity, radial return in Green–Lagrange strain, linear isotropic hardening,
  element erosion at the equivalent plastic strain `eps_fail` (0 = never).
* `concrete`: isotropic damage driven by the principal effective stresses — Rankine tension
  (`ft`, `Gf`) and compression crushing (`fc`, `Gc`), exponential softening, crack-band
  regularisation with the element size (the dissipated energy per crack area is `Gf`,
  independent of the mesh). Erosion at damage `derode` (default 0.99).

All materials: erosion of strongly distorted elements (`J < erode_J` or `J > 1/erode_J`).

**Failure and debris.** An element fails when at least half of its integration points fail.
Its mass stays on the nodes. Nodes without any intact element become debris particles: they
fall, collide and are carried by the flow (drag + buoyancy).

**Contact.** Penalty contact of the nodes with a horizontal ground plane and between all surface
nodes / debris particles that do not share an intact element (fragments, separate bodies),
regularised Coulomb friction. The penalty stiffness is scaled with the element frequency
(`k = kfac m (c_p/h)^2`), stable with the solid time step.

**Coupling with CFD (per RK stage).**

1. Lagrangian points on all exposed faces of intact elements (spacing below the smallest fluid cell)
   carry the structural surface velocity. The direct forcing `(u_s - u)/(alpha dt)` is spread with the
   3-point Roma kernel onto the staggered forcing fields, as for the FSI strips.
2. Debris particles get a quadratic drag `0.5 rho Cd A |u-v| (u-v)`; its reaction is spread onto the fluid.
3. Final stage, loads on the solid (`loads reaction`, default): the fluid parcel in the forcing volume
   of every surface point (mass `rho dV`, velocity `u_f` sampled with the kernel) is attached to the face
   nodes for the step: node and parcels are merged inelastically at the start of the step and move together
   through the substeps, `(m - rho_f V + M) a = F + (m - rho_f V) g`. The momentum handed over by the
   parcels is the fluid load (pressure, impact and inertia of the surrounding fluid). The fluid enclosed by
   the immersed boundary is accounted for per node (`m - rho_f V`: buoyancy and its inertia, as in
   Uhlmann's method).
   `loads pressure`: the pressure is probed `pressure_offset` cells outside every surface point and
   integrated as `-(p - p_gage) n dA` (explicit, only for heavy, stiff structures, see below).
4. The solid is advanced over the fluid time step with subcycling (staggered coupling).

Why not plain pressure loads: for thin, flexible structures the explicit pressure coupling is unstable
whenever the fluid added mass exceeds the structural mass; the relevant ratio is not rho_s/rho_f but
roughly `rho_s t / (rho_f L)` (thickness t, span L). The elastic obstacle test (rho_s/rho_f = 2.5,
t/L = 0.15) blows up with `loads pressure` within 15 ms of the impact. Two other variants were tried
and dropped: a lagged added-mass correction (pumps energy into the higher modes) and an implicit
penalty tie to the fluid velocity (Peskin-type odd-even instability once the fluid step grows,
because the elastic force reaches the fluid only once per fluid step). Attaching the local fluid
mass to the nodes during the substeps removes both.

## Input: `fem.dat`

One keyword per line, `#` starts a comment, SI units.

```
lattice   hx hy hz [ox oy oz]          # voxel size, optional lattice origin

material  id elastic   rho E nu
material  id plastic   rho E nu sigma_y H eps_fail
material  id concrete  rho E nu ft Gf fc Gc [derode]

box       x0 x1 y0 y1 z0 z1 id         # fill voxels with material id
stl       file.stl id                  # closed STL surface (ASCII or binary)
remove    x0 x1 y0 y1 z0 z1            # carve voxels (openings, windows)
remove_stl file.stl

fix       x0 x1 y0 y1 z0 z1 [xyz]      # support: fix the dofs of the nodes in the box (default xyz)

element   reduced|full [hg]            # default reduced, hourglass coefficient 0.1
cfl       0.5                          # solid time step factor
damping   alpha                        # mass-proportional damping [1/s], default 0
relax     T alpha                      # extra damping alpha for t < T (settle under gravity)
bulk_viscosity q1 q2                   # default 0.06 1.2
erode_J   0.05                         # distortion erosion threshold
ground    z [kfac mu]                  # ground plane, default kfac 1, mu 0.5
contact   on|off [kfac mu dist]        # part/debris contact, default on 1 0.5 0.8 (dist in h_min)
plane_strain 0|1                       # automatic in 2D runs
gravity   gx gy gz                     # overrides W 20-22

monitor   name x y z                   # time series of the nearest node

# coupling
loads     reaction|pressure            # default reaction (attached fluid parcels), pressure: explicit pressure integration
pressure_offset 1.5                    # probe distance in fluid cells (loads pressure)
points_per_face 0                      # Lagrangian points per face edge, 0: automatic
forcing   1                            # 0: no velocity forcing (one-way loads)
debris_cd 1.0
debris_reaction 1
print     dt                           # VTU interval, overrides Z 31
```

Material parameters for orientation:

| | rho | E | nu | other |
|---|---|---|---|---|
| normal concrete C30 | 2400 | 3.0e10 | 0.2 | ft 2.9e6, Gf 120, fc 3.0e7, Gc 1.5e4 |
| structural steel S355 | 7850 | 2.1e11 | 0.3 | sigma_y 3.55e8, H 1.0e9, eps_fail 0.15 |
| timber (along grain, elastic) | 500 | 1.1e10 | 0.3 | |
| rubber gate (Antoci et al. 2007) | 1100 | 1.2e7 | 0.4 | |

The crack band needs `h < 2 E Gf / ft^2` (about 0.8 m for C30), otherwise the input is rejected.

## Output: `REEF3D_FEM/`

* `REEF3D-FEM-*.vtu`: intact elements (hexahedra) and debris particles (vertices) on the deformed
  geometry. Point data `displacement`, `velocity`, `load` (fluid load per node); cell data `vonMises`
  (Cauchy stress, Pa), `damage`, `plastic_strain`, `material` (-1 for debris).
* `REEF3D_FEM_log.dat`: per fluid step: time, substeps, intact elements, eroded elements, debris particles,
  total fluid force, force on the supports, max von Mises stress, max displacement, kinetic energy,
  dissipated energy (damage + plasticity).
* `REEF3D_FEM_monitor_<name>.dat`: displacement and velocity of the monitored node.

## Verification (`tests/fem`)

`tests/fem/fem_test.cpp` checks the solid solver alone (no MPI, no REEF3D):

| test | result |
|---|---|
| cantilever under self weight, L/t = 10, 4 elements over the thickness | reduced: 1.006 × Timoshenko, refined 1.001; full: 0.969 |
| cantilever first natural frequency | reduced 0.990, full 1.008 × Euler–Bernoulli |
| free spinning block, 1.6 revolutions | energy drift 2e-6, no spurious strain |
| J2 uniaxial tension | hardening branch within 0.4 % |
| concrete bar with a weak band, two meshes | dissipated energy 1.000 and 0.998 × Gf A, peak = ft |
| elastic block dropped on the ground | comes to rest, no energy gain, penetration 0.2 mm |
| weak concrete cantilever collapsing under self weight | fails at the root, fragments rest on the ground |

```
cd tests/fem
g++ -O2 -std=c++20 -I../../ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY -I../../src fem_test.cpp ../../src/fem_solid*.cpp -o fem_test
./fem_test            # or: ./fem_test cantilever|freq|rotation|j2|crackband|drop|collapse
```

## Limitations and next steps

* The coupling is staggered (one exchange per fluid step). The enclosed-fluid correction
  limits the effective mass to 10 % of the solid mass, so bodies lighter than about 1.1 rho_f
  (floating structures, timber) are not supported yet (would need sub-iterations).
* The parallel path (owner sampling, one `MPI_Allreduce`, spreading at subdomain borders) is written
  like the FSI strips and rod trees but was only run on one rank so far: compare `M 10 1` and
  `M 10 4` on the test cases before production runs.
* The solid is replicated on every rank and computed serially: fine up to some 10^5 elements,
  beyond that it needs OpenMP or a partition of the solid.
* The fluid inside thick bodies is only constrained near the surface (diffuse direct forcing); the
  level set can leak into closed bodies over long times. A solid level set from the deformed mesh
  for the ghost-cell method would remove this.
* `loads pressure` neglects shear stresses; `loads reaction` includes the no-slip reaction but its accuracy for skin friction is that of a diffuse immersed boundary.
* Voxel geometry: surfaces are stair-stepped at the lattice size.
* Debris particles are point masses with drag and buoyancy, they do not displace fluid.
* Contact: node-node and node-ground only (no node-face contact), plane ground only.
