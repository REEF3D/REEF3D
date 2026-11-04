# CORAL: control/

Controllers and actuators.

## Content

* `coral_turbine_controller.h/.cpp`: wind turbine controller of the NREL 5 MW baseline type (Jonkman et al. 2009, the DISCON routine): first-order filter of the generator speed; variable-speed generator torque in regions 1, 1.5, 2 (k omega^2), 2.5 (induction-machine line through rated torque) and 3 (constant power, also with the pitch above a threshold), with a torque limit and a rate limit; gain-scheduled PI collective pitch (GK = 1/(1 + theta/theta_K)) with anti-windup, pitch limits and a rate limit; the NREL 5 MW's values as defaults; PI gains from the drivetrain inertia, rated speed and power sensitivity (`pi_gains`); `controller_load` runs it once per converged step in the model and sets the generator torque and the pitch command through functions (zero-order hold, as a controller DLL), so the controller does not depend on the rotor's implementation; its state for checkpoints; for floating turbines (`floating`) the part of the generator speed that the nacelle's fore-aft motion causes is taken out of the pitch loop's speed error (parallel compensation: a first-order drivetrain model driven by the hub's velocity along the shaft, from the rotor's dQ/dU and dQ/dOmega), so a pitch loop tuned as onshore does not give negative damping of platform pitch

Keyword `controller` (`io/coral_input.h`), with `tune`: the gains for any rotor from CORAL's steady BEM (`aero/coral_rotor_tune.h`), with `floating` also the drivetrain model of the compensation (`rotor_sensitivities`).

## Planned

* set-point smoothing and peak shaving (minimum pitch from a wind speed estimate), as ROSCO
* individual pitch control; ROSCO's floating feedback (nacelle velocity directly into the pitch) for comparison
* a Bladed-style DISCON interface to run ROSCO and other controller DLLs
* yaw control, brakes, start-up, shut-down and fault events
* structural control: tuned mass dampers
* wave energy controllers: damping, reactive, latching, declutching, model predictive control

## May include from

`core`
