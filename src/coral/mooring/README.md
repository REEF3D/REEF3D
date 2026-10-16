# CORAL: mooring/

Mooring lines and cables.

## Content

* `coral_catenary.h/.cpp`: elastic catenary with a seabed portion and friction (fairlead forces from the fairlead position, profile, anchor forces)
* `coral_mooring.h/.cpp`: quasi-static line (`mooring_qs`) and dynamic lumped-mass line (`mooring_line`: axial stiffness and damping, normal and tangential drag and added mass, waves and current, seabed contact and friction, explicit RK4 sub-steps inside the multibody step, start from the catenary, settling)

## Planned

* buoys and clump weights on the lines, multi-segment lines with different properties, line-to-line connections
* winches and line length changes, line breaking
* cables and risers with bending stiffness: geometrically exact rod in `structure/`

## May include from

`core`, `structure`, `hydro` (wave kinematics)
