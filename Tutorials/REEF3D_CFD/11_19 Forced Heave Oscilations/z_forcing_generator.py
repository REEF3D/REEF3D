import math

frequency = 1.0   # Hz
amplitude = 5000.0
duration = 20.0   # s
dt = 0.001         # s

with open("6DOF_forces.dat", "w") as f:
    t = 0.0

    while t <= duration:
        x = 0.0
        y = 0.0
        z = amplitude * math.sin(
            2.0 * math.pi * frequency * t
        )

        f.write(
            f"{t:.10e} "
            f"{x:.10e} "
            f"{y:.10e} "
            f"{z:.10e}\n"
        )

        t += dt
