# CORAL: soil/

Foundations.

## Content

* `coral_soil.h/.cpp`: soil profiles (sand, clay and linear Winkler layers below a mudline, effective stress), API RP 2GEO / 2A curves: p-y sand (Reese wedge and flow resistance, tanh) and soft clay (Matlock, static and cyclic), t-z sand and clay, Q-z; `soil_spring` at a pile node (lateral p-y radial about the pile axis, axial t-z, tip Q-z, torsion from the shaft friction), nonlinear elastic with analytic tangent; `add_pile_springs` along the nodes of a frame member with tributary lengths below the mudline

## Planned

* hysteresis and cyclic degradation (unloading-reloading, Masing), soil damping
* DNV and PISA curves, distributed moment and base shear springs for large-diameter monopiles
* stiff clay p-y (Reese, Welch), sand p-y after O'Neill and Murchison
* suction buckets and gravity bases as foundation matrices

## May include from

`core`, `structure`
