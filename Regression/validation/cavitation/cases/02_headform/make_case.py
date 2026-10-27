# Validation case cavitation 02: hemispherical head form (Rouse & McNown 1948), sheet cavitation
# behind the shoulder at zero yaw, Re = U D/nu = 1.36e5, D = 1 m, U = 1 m/s, rho = 1000.
#
#   python3 make_case.py [<D/dx> <ranks> <end time> [sigma ...]] [--pref P] [--amr] [--reboud N]
# VTU frames (P 35) over the last 5 s, the last 10 s for end times of 20 s or more (unsteady
# cavities); Betzy wall time 6 h up to 15 s, 12 h beyond (sigma 0.3 at D/40: dt 0.0035 s, 2.8 s
# per step on 128 ranks). No restart: the state files (P 40) do not hold the vapour fraction.
# VTU frames (P 35) over the last 5 s, the last 10 s for end times of 20 s or more (unsteady cavities);
# Betzy wall time 6 h up to 15 s, 12 h beyond (sigma 0.3 at D/40: dt 0.0035 s, 2.8 s per step on 128 ranks).
# Without arguments (e.g. Run in Spyder): D/20 on 8 ranks to 20 s for sigma = none (no
# cavitation), 0.5, 0.4, 0.3. Writes d<D/dx>/<name>/ with control.txt, ctrl.txt, headform.stl,
# job.slurm (Betzy); the body file is also written as floating.stl for X 180.
#
# Body: hemisphere of radius 0.5 centred at the origin (nose tip at x = -0.5), cylinder to x = 8,
# flat back. Surface coordinate from the tip: s = 0.5 theta on the nose, s = pi/4 + x on the
# cylinder (shoulder at s/D = 0.785). Domain [-4, 12] x [-4, 4]^2, inflow U = 1 at x = -4, the
# x+ side a wall with an outlet pressure patch p = 0 (B 441, B 412) that fixes the pressure level;
# p_v = P - sigma rho U^2/2 (W 61), P the upstream pressure: the drag lowers it against the outlet,
# -15.5 Pa at D/40 with the outlet patch, -3.5 to -5 Pa with the outflow B 77 10 of --amr (the
# defaults; --pref P instead). From 5 s on W 69 (cavitation number control) corrects p_v slowly
# (tau 5 s) to p_ref - sigma rho U^2/2, p_ref at the upstream probe point (-3, 3, 0.01). A fast
# control (tau 1 s from the start) took the start-up pressure into p_v and the growing cavity kept
# it up: 14 m^3 of vapour at D/20 + AMR. The actual sigma follows from the upstream probe and the
# p_v column of the cavitation log (headform_analysis.py).
# Grid (B 127-129, cell based): uniform dx in x in [-0.8, 2.5], |y|, |z| < 0.75, growing to 0.25.
# k-epsilon (T 10 1) with wall functions at the direct-forcing body: laminar runs separate at the
# shoulder (D/20: reverse flow to s/D = 2, Cp_min -0.34 instead of -0.74).
# --reboud N: Reboud correction of the eddy viscosity in the mixture (W 70 N, N = 10): k-epsilon
# damps the re-entrant jet at the closure of the cavity, the cavities are too long and steady
# (D/20 + AMR: 0.10, 0.13, 0.34 D at sigma 0.5, 0.4, 0.3); written to d<D/dx>[_amr]_reboud/.
# --amr: CFD mesh refinement (G 1 1) in the box G 15 x in [-0.75, 2.45], |y|, |z| < 0.7 (inside the
# uniform part of the base grid): the nose, the shoulder and the cavities of sigma >= 0.3 at
# dx/2, the rest of the domain at dx (D/40 + AMR: D/80 at the body, ~3.2 M patch cells). Needs
# (cfd_amr::scope) N 40 3 instead of 4, no patch boundary: the x+ side the outflow C 14 2 with
# B 77 10 (pressure 0 there, Dirichlet, as the patch of the uniform runs), the patches placed on
# the ranks for the load (G 40 1), one time step for all levels (G 7 0, required by W 60). Output
# into the LAGOON store (P 18 1, the level-0 volume and the patches with alpha_v) instead of the
# VTU files, every 0.5 s (D/40 + AMR: ~0.2 GB per output); analysis: headform_analysis.py reads the
# finest patches from the store.
import math, os, sys

D, U, RHO, RE = 1.0, 1.0, 1000.0, 1.36e5
NU = U*D/RE

CONTROL = """C 11 1
C 12 3
C 13 3
C 14 {c14}
C 15 3
C 16 3

B 1 {dx}
B 10 -4.0 12.0 -4.0 4.0 -4.0 4.0

B 101 11
B 102 11
B 103 11
B 127 {dx} 0.25 0.85 3.3 1.1
B 128 {dx} 0.25 0.0 1.5 1.1
B 129 {dx} 0.25 0.0 1.5 1.1

M 10 {ranks}
"""

CTRL = """B 60 1
B 61 1
{outlet}D 10 4
D 20 2
D 30 1

F 30 0
F 40 0

N 40 {n40}
N 41 {tend}
N 47 0.2

M 10 {ranks}

P 10 1
{p18}P 24 1
P 35 {tvtu} {tend} {dtv}
P 64 -3.0 3.0 0.01
P 64 -3.0 0.01 3.0

T 10 1

W 1 {rho}
W 2 {nu:.6e}
W 10 {q:.1f}
W 22 0.0
{cav}
X 10 1
X 11 0 0 0 0 0 0
X 21 1000.0
X 180 1
{amr}"""

JOB = """#!/bin/bash
# Betzy: head form {name}, D/{res}, {ranks} ranks. In this directory: control.txt, ctrl.txt,
# floating.stl (written by make_case.py), REEF3D and DiveMESH built on Betzy. sbatch job.slurm
#SBATCH --job-name=hf{res}{tag}_{name}
#SBATCH --account=nnXXXXk
#SBATCH --partition=normal
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks}
#SBATCH --time={walltime}
#SBATCH --output=run_%j.out

set -e
for f in control.txt ctrl.txt floating.stl REEF3D DiveMESH; do
    [ -e "$f" ] || {{ echo "missing $f in $(pwd)"; exit 1; }}
done

module purge
module load foss/2023b

./DiveMESH > divemesh.log 2>&1
srun ./REEF3D > reef3d.log 2>&1
"""


def stl(fn, nth=72, nphi=24, xend=8.0):
    # closed triangulated body, normals outwards
    R = 0.5*D
    tris = []

    def ring(x, r):
        return [(x, r*math.cos(2*math.pi*j/nth), r*math.sin(2*math.pi*j/nth)) for j in range(nth)]

    # nose: rings from the tip (phi = 0) to the shoulder (phi = pi/2)
    rings = []
    for i in range(1, nphi+1):
        ph = 0.5*math.pi*i/nphi
        rings.append(ring(-R*math.cos(ph), R*math.sin(ph)))
    tip = (-R, 0.0, 0.0)
    for j in range(nth):
        tris.append((tip, rings[0][(j+1) % nth], rings[0][j]))
    # cylinder: axial rings at about the circumferential spacing
    ncyl = int(round(xend/(2*math.pi*R/nth)))
    for i in range(1, ncyl+1):
        rings.append(ring(xend*i/ncyl, R))
    for a, b in zip(rings[:-1], rings[1:]):
        for j in range(nth):
            j1 = (j+1) % nth
            tris.append((a[j], a[j1], b[j1]))
            tris.append((a[j], b[j1], b[j]))
    # back disk
    c = (xend, 0.0, 0.0)
    last = rings[-1]
    for j in range(nth):
        tris.append((c, last[j], last[(j+1) % nth]))

    def nrm(t):
        (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) = t
        ax, ay, az = x1-x0, y1-y0, z1-z0
        bx, by, bz = x2-x0, y2-y0, z2-z0
        nx, ny, nz = ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx
        l = math.sqrt(nx*nx+ny*ny+nz*nz) or 1.0
        return nx/l, ny/l, nz/l

    with open(fn, "w") as f:
        f.write("solid headform\n")
        for t in tris:
            f.write("facet normal %.6e %.6e %.6e\n outer loop\n" % nrm(t))
            for v in t:
                f.write("  vertex %.8e %.8e %.8e\n" % v)
            f.write(" endloop\nendfacet\n")
        f.write("endsolid headform\n")
    return len(tris)


def make(argv):
    pref = 0.0
    pref_given = "--pref" in argv
    if "--pref" in argv:
        k = argv.index("--pref")
        pref = float(argv[k+1])
        argv = argv[:k] + argv[k+2:]
    amr = "--amr" in argv
    argv = [a for a in argv if a != "--amr"]
    reboud = 0.0
    if "--reboud" in argv:
        k = argv.index("--reboud")
        reboud = float(argv[k+1])
        argv = argv[:k] + argv[k+2:]
    res = int(argv[0])
    ranks = int(argv[1])
    tend = float(argv[2])
    sigmas = argv[3:] if len(argv) > 3 else ["none", "0.5", "0.4", "0.3"]
    dx = D/res
    q = 0.5*RHO*U*U
    for sg in sigmas:
        name = "nocav" if sg == "none" else "sigma" + sg
        d = os.path.join(f"d{res}" + ("_amr" if amr else "") + ("_reboud" if reboud > 0 else ""), name)
        os.makedirs(d, exist_ok=True)
        # rho_l/rho_v = 1000 (W 62 1, nu_v 1e-5: mu_v as for water vapour): with the real ratio
        # 58000 the D/40 run diverged at t = 3.3 s (lateral velocity 567 m/s in the vapour);
        # at D/20 the two ratios give the same Cp within 0.006 and the same cavity
        # W 61 from the upstream pressure of the set-up, W 69 corrects it slowly from 5 s on
        p0 = pref if pref_given else (-4.5 if amr else -15.5)
        cav = "" if sg == "none" else (f"W 60 1\nW 61 {p0 - float(sg)*q:.1f}\nW 62 1.0\nW 63 1.0e-5\n"
                                       f"W 69 {float(sg)} {U} -3.0 3.0 0.01 5.0\n"
                                       + (f"W 70 {reboud}\n" if reboud > 0 else ""))
        with open(f"{d}/control.txt", "w") as f:
            f.write(CONTROL.format(dx=dx, ranks=ranks, c14=2 if amr else 21))
        with open(f"{d}/ctrl.txt", "w") as f:
            # W 10: discharge through the inlet face, U times its area 8 x 8
            outlet = "B 77 10\n" if amr else "B 441 1 4 11.99 12.01 -4.0 4.0 -4.0 4.0\nB 412 1 0.0\n"
            amrkeys = "G 1 1\nG 7 0\nG 15 -0.75 2.45 -0.7 0.7 -0.7 0.7\nG 40 1\n" if amr else ""
            f.write(CTRL.format(tend=tend, tvtu=max(0.0, tend-(10.0 if tend >= 20.0 else 5.0)), ranks=ranks, rho=RHO, nu=NU,
                                q=U*64.0, cav=cav, outlet=outlet, dtv=0.5 if amr else 0.25, n40=3 if amr else 4, p18="P 18 1\n" if amr else "",
                                amr=amrkeys))
        n = stl(f"{d}/floating.stl")
        with open(f"{d}/job.slurm", "w") as f:
            # AMR: ~3.2 M patch cells at D/40 and half the time step: about 8 times the uniform run
            wall = ("24:00:00" if tend <= 15.0 else "48:00:00") if amr else ("06:00:00" if tend <= 15.0 else "12:00:00")
            f.write(JOB.format(name=name, res=res, ranks=ranks, nodes=max(1, ranks//128),
                               tasks=min(ranks, 128), walltime=wall, tag="amr" if amr else ""))
        print(f"{d}: dx = {dx}, {n} triangles, p_v = {p0 - float(sg)*q if cav else 'no cavitation'}")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    make(sys.argv[1:] if len(sys.argv) > 1 else ["20", "8", "20"])
