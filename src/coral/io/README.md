# CORAL: io/

Input and output.

## Content

* `coral_input.h/.cpp`: the CORAL input file (keywords documented in the header): bodies, joints, waves, hydrodynamic databases (CORAL format or WAMIT), PTOs, springs, linear 6x6 elements, quasi-static and dynamic mooring lines, Morison members, time and output settings; errors name the line
  and the structure: cross-sections (`section`, tube or properties), beam members (`member`), point masses, static equilibrium (`static`), linear analysis (`modes`), VTP output (`output vtp`)
* `coral_output.h/.cpp`: time series of the bodies and of the output channels of all loads
* `coral_vtp.h/.cpp`: the structure for ParaView: bodies as vertices, beam elements and mooring lines as lines
* `coral_modes_out.h/.cpp`: `modes.dat` (frequency, period, damping and the split of the modal kinetic energy over the six global directions, so fore-aft, side-side and torsion modes can be told apart) and `mode_shapes.dat`

## Planned

* readers of the data formats CORAL uses (airfoil polars, cross-section properties, hydrodynamic databases, wind fields, STL)
* mode shapes as VTP animations
* restart files

## May include from

all physics folders
