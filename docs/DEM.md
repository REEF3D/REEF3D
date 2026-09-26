# REEF3D::DEM

Standalone discrete element module for rigid particles of arbitrary shape, coupled to
REEF3D::CFD and REEF3D::NHFLOW. It does not depend on the 6DOF module or on the
particle container (`part`) of the sediment/CPM code.

## Method

* **Shapes (level-set DEM).** Every shape template has a signed distance function in its
  body frame (analytic for sphere, box and cylinder; a grid SDF from the triangulation for
  ellipsoids and STL files), a set of surface nodes and a volume quadrature. Mass, centroid and
  principal inertia are computed from the triangulation (Eberly), and the template is stored in
  its principal frame. All particles of one shape share the template.
* **Contacts.** Surface nodes of one particle are tested against the level set of the other,
  both ways; sphere centres are tested directly. Walls are planes (`dem.txt`) and, optionally, the
  REEF3D level sets: `topo` and `solid` in CFD, the bed and `SOLID` in NHFLOW. Contact patches are
  reduced to at most `E 21` points per body pair and normal cluster (farthest-point sampling
  plus the deepest point).
* **Time stepping: non-smooth contact dynamics** (Moreau-Jean). Contacts are hard and resolved
  as velocity-level impulses. A contact is active if its midpoint-predicted gap is closed (Newton
  impact law with restitution). Otherwise it is speculative (no penetration within the step;
  the approach velocity is kept for the restitution in the next step). Friction uses a Coulomb
  disc projection with an isotropic tangential mass, so the converged solution satisfies maximum
  dissipation. Contacts are solved with symmetric projected Gauss-Seidel and warm-started from the
  previous step. Penetration is removed with split impulses, which move positions but add no
  energy. The gyroscopic term is implicit (Catto 2015). The step is limited by kinematics
  (`E 20`), not by contact stiffness, so the DEM normally runs at the fluid time step or with a
  few substeps.
* **Parallelisation (domain decomposition).** A particle whose bounding radius is at most `E 24`
  times the smallest subdomain extent is *distributed*. It is owned by the rank that holds its
  centroid, and it moves to a new owner after each DEM substep. Copies (ghosts) are sent to every
  partner rank whose subdomain lies within the particle's halo, which covers the contact range
  plus the reach of its fluid coupling. Larger particles are *replicated* on all ranks. A contact
  is solved by one rank: the owner of the lower global id, the owner of the distributed particle
  in a mixed pair, or rank 0 for two replicated particles.
  - **Solver:** contacts of particles solved on one rank only are swept locally. Contacts of
    particles shared by several ranks are swept in colour phases: ranks closer than twice the halo
    get different colours, and one colour sweeps at a time. After each phase, the velocity
    corrections of ghosts go back to their owners and the new velocities out to the ghosts. The
    warm-start impulses are synchronised once before the first sweep.
    - For distributed particles this is a Gauss-Seidel iteration with a different contact order,
      not a Jacobi iteration between ranks.
    - Replicated particles touched by several ranks are mass-split (Tonge et al. 2012), and their
      corrections are averaged.
    - Convergence is tested globally.
    - Summed ghost corrections (block Jacobi) and mass splitting for all particles were tried
      first; they did not converge for stacks across rank boundaries.
  - **Fluid coupling:** fluid data, forcing integrals, kernel sums and wall contacts are
    evaluated where the grid point lies and returned to the particle's owner by point-to-point
    exchange with the partner ranks. Only replicated particles need global reductions.
  - **Consequences:** memory and work scale with the local particles plus their halo. With
    `E 24 0` all particles are replicated. Because contact forces in piles are statically
    indeterminate, results depend slightly on the decomposition.

## Fluid coupling (`E 11`)

| E 11 | mode | description |
|---|---|---|
| 0 | dry | gravity and contacts only |
| 1 | unresolved | particles smaller than the grid |
| 2 | resolved | particles larger than the grid |
| 3 | hybrid | per particle: resolved if d_eq/dx >= `E 12`, else unresolved |

**Unresolved.** Drag follows Haider & Levenspiel (1989), with the particle's sphericity computed
from the template and the Di Felice (1994) voidage correction. Drag is implicit in the DEM step
through the effective mass `m + m_a + dt K`, so particles with a response time shorter than the
fluid step stay stable. Buoyancy is integrated over the volume quadrature with the local fluid
density (CFD: two-phase level set) or the local free surface (NHFLOW), which also gives the
righting moment of floating particles. Added mass is `C_a = E 17`, and the fluid acceleration
force is optional (`E 22`). The reaction force (drag and added mass) is spread to the fluid momentum
with a compact quartic kernel of radius `max(d_eq, E 23 dx)`, normalised globally so the momentum
exchange is exact. The solid volume fraction is spread with the same kernel and enters the voidage
correction.

**Resolved.** The fluid velocity inside the smoothed particle indicator (Heaviside half width
`E 18 dx`) is forced to the rigid-body velocity in every RK stage. The hydrodynamic force is the
volume integral of the forcing, averaged with the weights of the RK stages (SSP-RK3, RK2 and
low-storage RK3 are recognised). In NHFLOW the forcing after the projection is included as well;
it carries the non-hydrostatic pressure. To that is added the
rate of change of the fluid momentum inside the particle, taken from the fluid field (Kempe &
Froehlich 2012), which is stable at lower density ratios than the rigid-body estimate. In CFD the
forced fluid is subject to gravity, so the buoyancy correction is `-m_f g` with
`m_f = rho_mean V`. Fixed resolved particles act as obstacles.

**NHFLOW specifics.** A forced particle has no free surface inside it, so resolved forcing is used
only while a particle is fully submerged (below the ambient water level). A surface-piercing
particle larger than the grid switches automatically to unresolved loads. Its hydrostatic buoyancy
uses the ambient water level, sampled on a ring around the particle, and it is coupled one-way
(fluid to particle). Radiation damping of large floating bodies is therefore not represented.

The fluid forcing is applied in the momentum schemes that use `momentum_forcing` (for example
`N 40 3`) and in NHFLOW. With `N 40 14` the particles feel the fluid, but they do not force it back.

## Control file (`ctrl.txt`)

| key | type | default | meaning |
|---|---|---|---|
| E 10 | int | 0 | DEM on (1), input in `dem.txt` |
| E 11 | int | 1 | coupling: 0 dry, 1 unresolved, 2 resolved, 3 hybrid |
| E 12 | double | 5.0 | hybrid threshold d_eq/dx |
| E 13 | int | 1000 | max contact solver iterations |
| E 14 | double | 1e-5 | solver tolerance (velocity change per sweep relative to max(v_max, 1 mm/s)) |
| E 15 | double | 0.1 | print interval [s] for VTP and state file, 0 off |
| E 16 | int | 1 | minimum DEM substeps per fluid step |
| E 17 | double | 0.5 | added mass coefficient (unresolved) |
| E 18 | double | 1.5 | resolved: Heaviside half width in cells |
| E 19 | double | 0.2 | penetration correction factor (split impulse) |
| E 20 | double | 0.25 | max travel per DEM step, fraction of the smallest bounding radius |
| E 21 | int | 6 | contact points per pair and normal cluster, 0 keeps all |
| E 22 | int | 1 | unresolved: fluid acceleration force on/off |
| E 23 | double | 2.0 | unresolved: kernel radius in cells |
| E 24 | double | 0.25 | distributed particles: max bounding radius as fraction of the smallest subdomain extent, larger ones are replicated; 0 replicates all |

Gravity is taken from `W 20-22`, so `W 22 -9.81` has to be set. For CFD cases in still water,
initialise the pressure (`I 10 1`), otherwise the start-up flow disturbs the particles.

## Particle input (`dem.txt`)

Lines starting with `#` are comments. Angles are in degrees (rotation about x, then y, then z).

```
material <id> <density> <friction> <restitution>
wall_material <friction> <restitution>
shape <id> sphere <r>
shape <id> box <lx> <ly> <lz>
shape <id> cylinder <r> <length>                 # axis along body z
shape <id> ellipsoid <a> <b> <c>
shape <id> stl <file> [scale]                    # ASCII or binary, closed surface
node_spacing <f>           # surface node spacing as fraction of d_eq (default 0.15)
sdf_resolution <n>         # grid SDF cells along the longest axis (default 32)
quadrature <n>             # volume quadrature points along the longest axis (default 6)
grid_walls <0|1>           # contacts with the REEF3D topo/solid/bed level sets (default 1)
wall <nx> <ny> <nz> <px> <py> <pz>               # plane, normal pointing into the domain
particle <shape> <mat> <x> <y> <z> [<rx> <ry> <rz> [<u> <v> <w>]]
fixed <shape> <mat> <x> <y> <z> [<rx> <ry> <rz>]  # immovable, acts as obstacle when resolved
block <shape> <mat> <nx> <ny> <nz> <x0> <y0> <z0> <dx> <dy> <dz> [<random rotation 0|1> [<seed>]]
```

In 2D cases (one cell in y) the particles move in the x-z plane, and the y offset is ignored in
all particle-grid distances. The particles keep their 3D mass, and the fluid feedback is
distributed over the one-cell-wide slab. The resolved forcing integral covers only the slab, so in
2D, resolved particles should be extruded shapes whose y extent equals the cell width.

## Output

* `REEF3D_DEM_VTP/REEF3D-DEM-*.vtp`: particle surfaces with velocity, id and coupling mode
  (ParaView).
* `REEF3D_DEM/REEF3D-DEM-state.dat`: time, id, position, quaternion, velocity, angular velocity,
  active flag, mode, hydrodynamic force and fluid velocity at the particle.
* Screen: substeps, contacts, solver iterations, residual, maximum penetration, kinetic energy
  and DEM time per step.

## Verification (hans_dev, September 2026)

* Core: the rebound velocity of a sphere with e = 0.5 is 0.3 %, 4 % and 11 % low for dt = 1, 10
  and 30 ms. A box sliding down a 30 degree incline with mu = 0.3 matches g(sin - mu cos) exactly,
  and it sticks for mu = 0.7. A sphere launched on a plane rolls at 5/7 v0. A stack of 5 boxes stays
  at rest at dt = 10 ms; a stack of 10 is quasi-static (KE ~ 1e-5 J, the solver runs at its
  iteration limit). 100 mixed particles (spheres, boxes, cylinders, ellipsoids) settle into a pile
  without escapes at dt = 5 ms.
* CFD unresolved: a 3 mm glass sphere reaches a terminal velocity of 0.3515 m/s
  (Haider-Levenspiel: 0.3513 m/s), and its hydrodynamic force balances its weight to 4 digits.
* CFD resolved (d/dx = 5, Re ~ 8000, narrow tank): the sphere is stable at density ratio 1.25. Its
  terminal velocity is about half the unbounded value, so resolved coupling needs finer grids
  (d/dx >= 10) to be quantitative.
* NHFLOW: a submerged resolved sphere settles onto the bed and rests with the correct buoyancy.
  Floating unresolved and surface-piercing boxes reach the correct draft. In regular waves,
  floating blocks follow the orbital motion (with `E 22 1`).

**Domain decomposition** (September 2026, 1/2/4 ranks):
- Settling sphere on a rank boundary: same terminal velocity to 10 digits on 1 and 2 ranks.
- NHFLOW floating/submerged case: same mean positions on 1 and 2 ranks.
- 7-box stack across the boundary with alternating owners: at rest on 1, 2 and 4 ranks, with the
  similar iteration counts (about 60-600 once settled) and less than 0.01 mm drift in 1.5 s.
- Replicated box resting on distributed spheres across two ranks: same resting height as on one
  rank.
- 250-particle pile: same final statistics on 1, 2 and 4 ranks.

## Limitations and next steps

* Projected Gauss-Seidel converges slowly for large piles; accelerated projected gradient
  solvers, or sleeping of resting particles, would reduce the cost.
* Replicated (large) particles: their contacts with each other are solved on rank 0, and their
  velocity corrections are summed globally at every synchronisation, so keep their number small.
* Distributed solver: every sweep needs one synchronisation per colour. That is typically 2-8
  colours for a Cartesian decomposition, and more when the halo is wider than a subdomain.
  Warm-start caches stay on the rank that solved a contact, so a contact that changes rank starts
  cold once. Particles whose centroid leaves all subdomains, or that jump beyond the partner ranks
  in one substep, stay with their owner.
* The unresolved coupling has no correction for the self-induced velocity. It is small for
  d/dx < 0.5 with the default kernel.
* NHFLOW surface-piercing particles are coupled one-way, without radiation damping.
* With several DEM substeps per fluid step, the reaction force returned to the fluid uses the
  end-of-step particle velocity, not the impulse accumulated over the substeps.
* When an NHFLOW particle switches from unresolved to resolved, it has no forcing integral for
  one step.
* No particle-6DOF contacts yet.
