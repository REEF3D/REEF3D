# REEF3D RodTree — flexible branching rods (soft corals, gorgonians, vegetation)

Two-way coupled FSI of flexible, branching slender structures in REEF3D::NHFLOW and REEF3D::CFD.
Phase 1 (structural solver) and Phase 2 (unresolved coupling) of the soft-coral FSI roadmap.

## Activation (ctrl.txt)

| Flag | Type | Meaning | Default |
|---|---|---|---|
| `Z 20` | int | 0 off, 1 unresolved (Morison / actuator-line) two-way coupling; input file `rodtree.dat` in the case directory | 0 |
| `Z 21` | double | print interval of `REEF3D_RodTree/REEF3D-RodTree-*.vtp` in seconds, 0 = no VTP | 0.0 |

Works with NHFLOW (all RK schemes, via `nhflow_forcing::forcing`) and with the CFD momentum schemes that call
`momentum_forcing_start` (FC2, FC3, FCLS3, FCC3, RK2, RK3, RK3CN, RKLS3, PLIC variants).
Not hooked into `momentum_RKLS3_df/sf`.

Density `W 1` and gravity `W 20/21/22` are taken from ctrl.txt (set `W 22 -9.81`).

## Input: `rodtree.dat`

```
# global options (optional, outside colony blocks)
integrator implicit          # implicit (default) | explicit (verification only, needs damping 0)
substeps   2                 # implicit sub-steps per fluid time step (default 1)
reaction   full              # full (drag+inertia, default) | drag | none (one-way)

colony fan1
  material  5.0e6 0.4 1100   # E [Pa], Poisson ratio, structural density [kg/m3]
  hydro     1.2 0.01 1.0     # Cdn, Cdt (skin friction), Ca  [Cm, default 1+Ca]
  damping   0.01             # Kelvin-Voigt time constant beta [s]
  polyps    0.004 0.5 1.0    # optional: polyp height h_p [m], frontal solidity phi_p, Cd_p  [k_t, default 1]
  polyp_response 0.15 0.1 1.0  # optional: retract above U_r, fully retracted at U_r+dU [m/s], time constant tau [s]
  refine    6                # every edge is split into 6 rod elements
  node 1  50.0 0.0 0.00 0.010   # id x y z radius
  node 2  50.0 0.0 0.40 0.008
  node 3  50.2 0.0 0.70 0.005
  node 4  49.8 0.0 0.70 0.005
  edge 1 2                   # parent node, child node (parent is closer to the root)
  edge 2 3
  edge 2 4
  clamp 1                    # root(s), clamped to the bed
  instance 0 0 0             # optional copies (dx dy dz), e.g. canopies
  instance 1.5 0 0
end
```

Each colony must be a tree (every node at most one parent edge), every root must be clamped.
The radius is interpolated linearly along an edge (taper).

## Model

**Structure.** Discrete Cosserat rod on a tree. Every element is a rigid cylinder with centre, velocity,
orientation quaternion and angular velocity. Neighbouring elements (and a root element and the bed) are
connected at the shared node by an elastic 6-DOF joint:

- translational penalty spring `k = E A / l_j` on the gap of the attachment points (extension, shear),
- rotational spring `K = EI/l_j (I - t t^T) + GJ/l_j t t^T` acting on the rotation vector of
  `q_rel q_rel0^-1` (bending and torsion, rest shape = input geometry, stress free at t = 0),
- Kelvin-Voigt damping `beta` on both,

with `l_j` the mean length of the two elements (half an element at a clamp). Branch points are ordinary joints
between the parent element and each child, so junction stiffness and torsional coupling come from the same model.
Converges with second order in the element length (see tests).

**Hydrodynamics (per element, relative velocity `u_r = u_f - v`, submerged fraction `chi`):**

- normal drag `0.5 rho Cdn D l |u_n| u_n`, tangential drag `0.5 rho Cdt pi D l |u_t| u_t`,
- inertia `rho V (Cm a_f - Ca a_s)` in the normal plane (Froude-Krylov + added mass),
- added mass `Ca rho V (I - t t^T)` and added rotational inertia in the element mass matrix,
- buoyancy and gravity.

**Polyps.** Extended polyps are a porous layer of height `h_p` and frontal solidity `phi_p` around the branch.
With extension `ext` (0 retracted .. 1 extended) they add a frontal width `b_p = 2 ext h_p phi_p` with drag
coefficient `Cd_p`:

- normal drag coefficient per length: `Cdn D + Cd_p b_p`
- tangential drag coefficient per length: `Cdt pi D + k_t Cd_p b_p` (`k_t` = axial/normal ratio of polyp drag, default 1)

Without `polyp_response` the polyps stay extended. With it, each element follows its own local normal relative
speed `|u_n|`: target extension 1 below `U_r`, 0 above `U_r + dU`, linear in between, and
`d ext/dt = (target - ext)/tau` (implicit, once per time step). Emerged elements retract. Sheltered branches in
the wake therefore keep their polyps out longer than exposed ones. The polyp drag enters the drag Jacobian, the
fluid reaction and the point-implicit fluid drag in the same way as the branch drag. Not modelled: polyp added
mass and volume, feeding currents, diel or other biological schedules.

`chi` is the fraction of the element's Lagrangian points that lie in water (NHFLOW: between bed and bed+WL,
CFD: level set `phi >= 0`), so emergent colonies and drying are handled.

**Time integration.** Linearly implicit Euler on the whole tree: joint Jacobians by central finite differences,
drag and gyroscopic Jacobians analytic, sparse LU (Eigen, COLAMD). Stable for the stiff joints; first order,
numerical damping ratio about `omega h / 2` per mode (use `substeps` for oscillatory problems with large fluid
time steps). Because only one Newton step is taken per sub-step, a sub-step is rejected and halved (roll-back)
if any element rotates by more than 0.25 rad or moves by more than 0.25 element lengths; this keeps very soft
colonies in strong currents stable at NHFLOW/CFD time steps (test 8). The explicit symplectic integrator sub-cycles automatically and is meant for verification: with
`beta > 0` its stable step collapses (stiffness-proportional damping) and the code stops with a message.

**Coupling.** Per RK stage: fluid velocity at the Lagrangian points (element split into
`ceil(l / min(dx,dy))` points, max 64) by the owner rank, one `MPI_Allreduce`; Morison loads; reaction spread with
the 4-point Peskin kernel (NHFLOW: horizontal kernel, vertical kernel renormalised over the wet column).
The fluid-side drag is point-implicit: the spread force is divided by `1 + alpha dt c S / (rho V_cell)` with `c`
the linearised drag slope and `S` the kernel self-weight, so dense canopies do not reverse the flow within a stage.
The structure is advanced once per time step in the final stage (staggered explicit coupling). The structure is
replicated on all ranks.

NHFLOW velocities are interpolated with column z-coordinates rebuilt from `ZP*WL+bed` rather than with
`lexer::ccipol4V`: `ccipol4V` reads `p->ZSP` in the neighbouring column, whose halo is stale at subdomain faces
(and it gave biased values for points on cell faces). With the local interpolation, 1-rank and 2-rank runs with a
colony straddling the partition agree to 1e-12 (sampled velocity, one-way) and 8e-5 of the response range
(two-way, 9 s of Stokes waves).

## Output (`REEF3D_RodTree/`)

- `REEF3D_RodTree_<colony>.dat`: time, tip displacement (xyz), base force = force of the colony on the bed (xyz),
  total hydrodynamic force (xyz), number of submerged elements, mean polyp extension. Written every time step.
- `REEF3D-RodTree-XXXXXXXX.vtp`: elements as lines with radius, velocity, hydrodynamic force, submergence,
  polyp extension
  (ParaView: *Tube* filter with *Vary radius → By scalar*).

## Verification (`tests/rodtree`, no MPI/hypre needed)

```
make -C tests/rodtree run
```

1. cantilever under self-weight vs `qL^4/8EI`: 0.07 % at 40 elements, second-order convergence
2. first bending frequency: 0.09 % (implicit), 0.12 % (explicit)
3. rigid rotation of a branching colony is stress free
4. explicit energy drift 0.1 % (undamped branching colony)
5. stem in steady current: implicit = explicit; drag reconfiguration factor 0.475 at Cauchy number ~19
6. branching colony in oscillatory flow: implicit vs explicit within 0.12 % of amplitude
7. large-deflection cantilever vs elastica: 0.08 %
8. very soft colony (Cauchy number ~50) in 0.5 m/s current at dt = 0.03 s: same steady state as dt = 0.001 s
9. polyps on a rigid stem: drag ratio exactly `(Cdn D + Cd_p 2 h_p phi_p)/(Cdn D)`; retraction follows the
   time constant (ext(1 s) = 0.6073 for tau = 2 s, dt = 0.01 s) and the drag returns to the bare value; half
   extension halves the polyp drag

Coupled smoke cases (mesh with DIVEMesh, then run REEF3D) in `tests/rodtree/cases/`:

- `nhflow_fan_waves`: 2D NHFLOW flume, 5th-order Stokes waves H = 1 m, T = 4.5 s, 4 m depth, 0.8 m fan at
  x = 40 m. Periodic response after the waves arrive: tip excursion about -0.36..+0.28 m, base force about +-1.9 N.
- `cfd_bush_current`: 3D CFD channel, 0.2 m depth, about 0.5 m/s, 0.15 m soft colony. Reconfigures to a steady
  tip deflection of 0.127 m, base force 0.056 N; one-way vs two-way differ by about 1 % (single thin colony).

## Limitations / next phases

- Unresolved coupling only: element diameters should be well below the cell size. Resolved IB coupling in CFD
  (Phase 3, `Z 20 2`, strongly coupled because rho_s/rho_f ~ 1) is not implemented yet.
- Fluid acceleration `a_f` is the local time derivative at the element (no convective part).
- No contact (branch-branch, branch-bed): Phase 5, via the DEM.
- Hydroskeleton inflation and active pulsing are not modelled (constant material); polyps only as drag.
- The buoyancy reaction is not applied to the fluid (hydrostatic part stays in the pressure field).
