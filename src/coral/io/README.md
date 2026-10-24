# CORAL: io/

Input and output.

## Content

* `coral_input.h/.cpp`: the CORAL input file (keywords documented in the header): bodies, joints, waves, hydrodynamic databases (CORAL format or WAMIT), PTOs, prescribed joint motion (`drive`) and force (`jforce`), springs, fluid bodies of REEF3D (`fluid_body` two-way or `one_way`, `fluid_member`, `fluid_kinematics`, `fluid_coupling` with `origin` and `record`), linear 6x6 elements, quasi-static and dynamic mooring lines, Morison members, time and output settings; errors name the line
  and the structure: cross-sections (`section`, tube or properties), beam members (`member`), superelements (`superelement`), point masses (on the body or at a point), soil profiles and pile springs (`soil`, `soil_layer`, `pile`), Morison loads on frame and superelement members (`morison ... member`), point forces (`force`), full 6x6 foundation matrices (`linear6 ... K`), static equilibrium (`static`), linear analysis (`modes`), VTP output (`output vtp`)
* `coral_coupling.h/.cpp`: the coupling of a CORAL case with a resolved fluid without REEF3D headers or MPI (`coral_reef3d` is a thin wrapper around it): the mapping of the 6DOF bodies to `fluid_body` and `fluid_member` loads, the fluid frame (`origin`), the kinematics points and samples, the one-way step (CORAL trails the fluid to the time of its last loads), the two-way advance with the predictor and the stage poses; every input from the fluid can go into a binary journal (`record`), and `coupling::replay` runs the same CORAL case, or a changed one with the same fluid points, from the journal alone (`CORAL coral.txt replay <journal>`)
* `coral_output.h/.cpp`: time series of the bodies and of the output channels of all loads
* `coral_vtp.h/.cpp`: the structure for ParaView: bodies as vertices, beam elements, superelement elements (recovered nodes) and mooring lines as lines
* `coral_modes_out.h/.cpp`: `modes.dat` (frequency, period, damping and the split of the modal kinetic energy over the six global directions, so fore-aft, side-side and torsion modes can be told apart) and `mode_shapes.dat`

## Planned

* readers of the data formats CORAL uses (airfoil polars, cross-section properties, hydrodynamic databases, wind fields, STL)
* mode shapes as VTP animations
* restart files

## May include from

all physics folders
