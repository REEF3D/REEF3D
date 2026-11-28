# 6DOF_motion-0.dat (X 240 1): t x y z of the body centre, lifted at speed v from t0 until the lift L (then at rest),
# rows every 5 ms
# python3 motion_file.py x y z v T t0 L 6DOF_motion-0.dat
# Keep the lifted body below the top of the domain: a body through a closed top (C 16 21 / 3) splits the air into two
# compartments that exchange air only through the material below the body, and the trapped air holds the material back.
import sys
x, y, z, v, T, t0, L = map(float, sys.argv[1:8])
with open(sys.argv[8], "w") as f:
    n = int(round(T/0.005))
    for k in range(n + 1):
        t = k*0.005
        f.write("%.6f %.6f %.6f %.8f\n" % (t, x, y, z + min(L, v*max(0.0, t - t0))))
