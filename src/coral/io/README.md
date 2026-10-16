# CORAL: io/

Input and output.

## Content

* `coral_input.h/.cpp`: the CORAL input file (keywords documented in the header): bodies, joints, waves, hydrodynamic databases (CORAL format or WAMIT), PTOs, springs, linear 6x6 elements, quasi-static and dynamic mooring lines, Morison members, time and output settings; errors name the line
* `coral_output.h/.cpp`: time series of the bodies and of the output channels of all loads

## Planned

* readers of the data formats CORAL uses (airfoil polars, cross-section properties, hydrodynamic databases, wind fields, STL)
* time series output, VTP output of the structure for ParaView
* restart files

## May include from

all physics folders
