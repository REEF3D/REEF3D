# CORAL: app/

Standalone CORAL executable: runs a CORAL model without a fluid grid (engineering hydrodynamics, analytical waves, wind fields).

## Content

* `coral_main.cpp`: `CORAL [input file]` reads the model (`io/`), finds the static equilibrium (`static`), runs the linear analysis (`modes`: table on the screen, `modes.dat`, `mode_shapes.dat`), runs it with the generalized-α integrator (a failed step is repeated with sub-steps) and writes the time series and the VTP files (`coral_vtp/`); with `end` not after the start only the analyses run; `CORAL <input file> replay <journal>` runs the case on the fluid input that a coupled REEF3D run recorded (`fluid_coupling record`), without the fluid
* build: `make coral` (bin/CORAL) or the CMake target `CORAL`; not part of the REEF3D executable and not in the CodeLite project (it has its own `main`)

## Planned

* load case batches (`loadcase/`), restart

## May include from

everything except `reef3d`
