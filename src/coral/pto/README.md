# CORAL: pto/

Power take-offs of wave energy converters.

## Content

* `coral_pto.h/.cpp`: PTO in a joint (axial force or torque about the axis): linear spring-damper with preload, Coulomb friction, force and power limits, end stops, latching and declutching; absorbed power and energy

## Planned

* hydraulic PTO: cylinders, accumulators, valves, motor
* electrical generator model (efficiency, delivered power)
* PTO between two points (linear actuator not in a joint)

## May include from

`core`, `structure`
