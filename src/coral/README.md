# REEF3D::CORAL

**C**oupled **O**ffshore **R**otor **A**nalysis and **L**oads: the structural solver of REEF3D for
offshore wind (fixed and floating) and wave energy, built natively in REEF3D.

CORAL holds everything that moves, deforms or carries a load: rigid and flexible bodies, joints,
mooring lines, power take-offs, controllers, foundations, and the aerodynamic and engineering
hydrodynamic loads on them. The fluid comes from REEF3D's own solvers (CFD, NHFLOW, FNPF) or,
without a fluid grid, from CORAL's own wave and wind kinematics.

Architect: Hans Bihs

## Folders

| folder | content |
|---|---|
| `core/` | types, rotations, the multibody model (bodies, nodes, DOF numbering), constraints and joints, assembly, time integration (generalized-α), linear solvers |
| `structure/` | structural elements: rigid body, point mass, spring/damper, nonlinear beam, linear frame, modal body, Craig-Bampton superelement, drivetrain |
| `linear/` | eigenanalysis, linearisation about an operating point, Campbell diagrams, multi-blade coordinate transform |
| `wind/` | wind fields: steady and sheared inflow, IEC gusts and direction changes, turbulence generation (Kaimal, Mann), field interpolation |
| `aero/` | airfoil polars, blade element momentum with corrections, dynamic inflow, dynamic stall, tower and nacelle influence; later vortex wake |
| `control/` | wind turbine controller (torque, pitch, peak shaving, floating feedback, individual pitch), pitch and yaw actuators, brakes, faults, tuned mass dampers, wave energy controllers |
| `hydro/` | engineering hydrodynamics: hydrostatics, Cummins equation and radiation state-space models, wave excitation, drift forces, Morison, MacCamy–Fuchs, slamming, nonlinear Froude–Krylov loads, the hydrodynamic database |
| `panel/` | linear radiation/diffraction panel solver that writes the hydrodynamic database (alternative to computing it with REEF3D's wave solvers) |
| `mooring/` | quasi-static and dynamic lines, seabed contact and friction, buoys, clump weights, anchors, winches |
| `pto/` | power take-offs: linear, Coulomb, hydraulic, electrical generator, end stops |
| `soil/` | foundations: p-y, t-z and Q-z springs, apparent fixity, 6x6 foundation matrices |
| `ice/` | ice loads on structures |
| `io/` | CORAL input file, model description, readers of the data formats CORAL uses, time-series and VTP output |
| `loadcase/` | load case definitions and batch runs (IEC 61400-1/-3 design load cases, wave energy power matrices) |
| `post/` | statistics, rainflow counting, damage-equivalent loads, extremes, power matrices |
| `reef3d/` | the link to REEF3D: driver hook, fluid interfaces (CFD, NHFLOW, FNPF), loads from the resolved fluid, body motion back to the fluid, adapters to existing REEF3D code (6DOF, beam, mooring, PTO) |
| `app/` | standalone CORAL executable (no fluid grid); not part of the REEF3D build |

Each folder has a README with its scope and planned content.

## Layers

```
            reef3d/                app/
               |                     |
   io/   loadcase/   post/  <--------+
               |
  structure/ linear/ wind/ aero/ control/ hydro/ panel/ mooring/ pto/ soil/ ice/
               |
             core/
```

A folder includes only from itself and the folders below it:

| folder | may include from |
|---|---|
| `core` | `core` |
| `structure`, `wind`, `control`, `hydro` | `core` |
| `linear`, `pto`, `soil`, `ice` | `core`, `structure` |
| `mooring` | `core`, `structure`, `hydro` (wave kinematics) |
| `aero` | `core`, `structure`, `wind` |
| `panel` | `core`, `hydro` |
| `io` | all physics folders |
| `loadcase`, `post` | all physics folders, `io` |
| `app` | everything except `reef3d` |
| `reef3d` | everything, and REEF3D's own headers |

Only `reef3d/` includes REEF3D headers (`lexer.h`, `fdm.h`, `ghostcell.h`, ...) and MPI.
Everything else builds without REEF3D and without MPI, so it can be unit tested and run on
its own. `Regression/unit/coral_layering_check.py` checks these rules.

## Conventions

* Namespace `coral`; REEF3D-facing classes in `reef3d/` are named `coral_...` like the rest of REEF3D.
* File names `coral_<topic>.h/.cpp`, unique across all folders (the build puts objects of several
  folders side by side).
* Includes with the path from `src/`: `#include"coral/core/coral_types.h"`.
* No non-constant static data (`Regression/unit/static_check.py` covers `src/coral/` as well).
* SI units, right-handed coordinates, z up, origin at the still water level unless stated.
* Every new kernel gets a unit test in `Regression/unit/` and, where a reference exists, a
  validation case.

## Reuse of existing REEF3D code

CORAL does not copy existing code. The geometrically exact beam (`beam.h`), the dynamic mooring
lines (`mooring_*`), the PTO elements and controllers (`6DOF_pto_*`) and the 6DOF force
integration on the body surface (`6DOF_obj_*`) are used through adapters in `reef3d/` until
they are moved into CORAL. Without CORAL switched on, REEF3D runs as before.

## Roadmap

1. Core, rigid bodies, engineering hydrodynamics, hydrodynamic coefficients from REEF3D, PTO and
   mooring. First check: wave energy converters (RM3, OSWEC).
2. Structure: co-rotational frames with analytic tangents and sparse solvers, Craig-Bampton
   superelements, static equilibrium, eigenanalysis and linearisation, soil (API p-y, t-z, Q-z),
   Morison loads on flexible members; floating-frame superelements and modal bodies to follow.
3. Wind, aerodynamics and the native controller. Checks: OC3, OC4, IEA 15 MW (fixed and floating).
4. Strong coupling with the resolved fluid (CFD, NHFLOW, FNPF), flexible immersed geometry.
5. Load cases and post-processing; vortex wake, wind farms, marine operations.
