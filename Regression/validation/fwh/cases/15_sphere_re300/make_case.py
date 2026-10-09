# Writes the inputs of validation case 15 (sphere, Re = 300, FW-H) for one resolution
# into the directory d<D/dx>/ : control.txt (DIVEMesh), ctrl.txt (REEF3D), 6DOF_motion.dat (the
# lateral kick of 0.05 D, as scripts/make_inputs.py shedding 0.05), job.slurm (Betzy):
#   python3 make_case.py <D/dx> <ranks> <end time> [animation: t0 t1 dt]
# e.g.  python3 make_case.py 10 8 150            (the local run)
#       python3 make_case.py 20 128 200
#       python3 make_case.py 40 512 200 180 200 0.33
# Without arguments (e.g. run from Spyder) it writes these three.
# The sphere (D = 1, Regression/cases/cfd_3d_heave_sphere_6dof/floating.stl) is copied into the
# run directory separately.
import math, os, sys

CONTROL = """C 11 1
C 12 3
C 13 3
C 14 2
C 15 3
C 16 3

B 1 {dx}
B 10 -5.0 15.0 -5.0 5.0 -5.0 5.0

B 101 11
B 102 11
B 103 11
B 127 {dx} {dxmax} {xfocus} {xwidth} {growth}
B 128 {dx} {dxmax} 0.0 {ywidth} {growth}
B 129 {dx} {dxmax} 0.0 {ywidth} {growth}

M 10 {ranks}
"""

CTRL_HEAD = """B 60 1
B 61 1
D 10 4
D 20 2
D 30 1

F 30 0
F 40 0

I 11 1

N 40 4
N 41 {tend}
N 47 0.3

M 10 {ranks}

P 10 {print3d}
"""

CTRL_BODY = """
T 10 0

W 1 1000.0
W 2 0.0033333333
W 10 100.0

X 10 1
X 11 0 2 0 0 0 0
X 21 1000.0
X 180 1
X 240 1

"""

JOB = """#!/bin/bash
# Sphere Re = 300, FW-H validation, D/{res}: {ranks} ranks. Run directory: control.txt, ctrl.txt,
# 6DOF_motion.dat (from this folder), floating.stl (Regression/cases/cfd_3d_heave_sphere_6dof),
# REEF3D and DiveMESH built on Betzy. sbatch job.slurm
#SBATCH --job-name=sphere_d{res}
#SBATCH --account=nnXXXXk
#SBATCH --partition=normal
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks}
#SBATCH --time={walltime}
#SBATCH --output=run_%j.out

set -e
for f in control.txt ctrl.txt floating.stl 6DOF_motion.dat REEF3D DiveMESH; do
    [ -e "$f" ] || {{ echo "missing $f in $(pwd)"; exit 1; }}
done

module purge
module load foss/2023b          # GCC + OpenMPI; the toolchain REEF3D and DIVEMesh were built with

./DiveMESH > divemesh.log 2>&1
srun ./REEF3D > reef3d.log 2>&1
"""


def make(argv):
    res, ranks, tend = int(argv[0]), int(argv[1]), float(argv[2])
    anim = [float(x) for x in argv[3:6]] if len(argv) >= 6 else None
    dx = 1.0/res
    growth = {10: 1.1, 20: 1.08, 40: 1.06}.get(res, 1.06)
    full = res >= 20                       # the Betzy runs get the observer planes and rings
    d = f"d{res}"
    os.makedirs(d, exist_ok=True)

    with open(f"{d}/control.txt", "w") as f:
        f.write(CONTROL.format(dx=dx, ranks=ranks, growth=growth,
                               dxmax=0.5 if res == 10 else 0.25,
                               xfocus=1.0 if res == 10 else 2.0,
                               xwidth=4.0 if res == 10 else 6.0,
                               ywidth=2.0 if res == 10 else 2.5))

    # observers: ring r = 5 and distances along y (1/r^2 -> 1/r) and z
    obs  = [(5*math.cos(math.radians(a)), 5*math.sin(math.radians(a)), 0.01) for a in range(0, 360, 15 if full else 30)]
    obs += [(0.0, r, 0.01) for r in (3.0, 10.0, 30.0, 100.0, 1.0e3, 5.0e3, 2.0e4, 1.0e5)]
    obs += [(0.0, 0.01, 10.0), (0.0, 0.01, 2.0e4)]
    if full:
        # acoustic directivity ring at 20 km, and observer planes z = 0 for the animation:
        # near field +-40 m (5 m), acoustic field +-40 km (5 km)
        obs += [(2.0e4*math.cos(math.radians(a)), 2.0e4*math.sin(math.radians(a)), 0.01) for a in range(0, 360, 15)]
        for L in (5.0, 5.0e3):
            obs += [(i*L, j*L, 0.01) for i in range(-8, 9) for j in range(-8, 9) if max(abs(i), abs(j)) > 0]
    probes = [(0.01, 3.0, 0.01), (0.01, 0.01, 3.0), (-3.0, 0.01, 0.01), (3.5, 3.5, 0.01)]

    ctrl = CTRL_HEAD.format(tend=tend, ranks=ranks, print3d=1 if anim else 0)
    if anim:
        ctrl += f"P 35 {anim[0]} {anim[1]} {anim[2]}\n"
    ctrl += CTRL_BODY
    ctrl += "".join(f"P 64 {x:.4f} {y:.4f} {z:.4f}\n" for x, y, z in probes)
    ctrl += "\nU 10 1\nU 20 -1.5 2.5 -1.5 1.5 -1.5 1.5\nU 50 1.0 0.0 0.0\n"
    if full:
        ctrl += "U 31 0.05\n"
    ctrl += "".join(f"U 30 {x:.6g} {y:.6g} {z:.6g}\n" for x, y, z in obs)
    with open(f"{d}/ctrl.txt", "w") as f:
        f.write(ctrl)

    # lateral kick y = A (1 - cos(2 pi t/4))/2, A = 0.05 m, 0 <= t <= 4 s (X 240 1, X 11 0 2 0 0 0 0)
    with open(f"{d}/6DOF_motion.dat", "w") as f:
        for i in range(4001):
            t = i*0.001
            f.write(f"{t:.4f} 0.0 {0.05*0.5*(1-math.cos(2*math.pi*t/4.0)):.10e}\n")

    with open(f"{d}/job.slurm", "w") as f:
        f.write(JOB.format(res=res, ranks=ranks, nodes=max(1, ranks//128), tasks=min(ranks, 128),
                           walltime="06:00:00" if res >= 40 else "03:00:00"))
    print(f"{d}: {len(obs)} observers, {len(probes)} probes")


os.chdir(os.path.dirname(os.path.abspath(__file__)))

if len(sys.argv) > 1:
    make(sys.argv[1:])
else:
    for argv in (["10", "8", "150"], ["20", "128", "200"], ["40", "512", "200", "180", "200", "0.33"]):
        make(argv)
