# P5d T2: the lower Nidelva and the harbour of Trondheim on the two-block grid of nidelva_harbour_grid.py: the river
# discharge at the upstream end of the corridor (block 0, B 60 1, W 10 Q), the water level at the fjord edge of the
# harbour block (G 54 1 2 WL), SFLOW or NHFLOW, optionally with one refinement level over the river and the harbour
# basins (G 1 1, G 10 boxes).  The bed is schematic (see nidelva_harbour_grid.py).
# usage: python3 t2_harbour.py geo.dat rundir   (environment: SOLVER sflow|nhflow, NP, TEND (s), OUT (s, VTS output),
#        Q, WL, AMR (1: refinement, AMRBOX the G 10 box), NZ (NHFLOW layers), A520 (NHFLOW pressure), SED (1: bed load with the bend terms),
#        XCTRL (;), NOGRID (reuse the grid of rundir))
import sys, os, re, json, shutil, subprocess

T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(T, "../sflow_curvilinear"))
from p1lib import REEF3D, DIVEMESH

geo, d = sys.argv[1], os.path.abspath(sys.argv[2])
solver = os.environ.get("SOLVER", "sflow")
np_ = int(os.environ.get("NP", "2"))
tend = float(os.environ.get("TEND", "1800"))
out = float(os.environ.get("OUT", "300"))
Q = float(os.environ.get("Q", "100"))
WL = float(os.environ.get("WL", "10"))
nz = int(os.environ.get("NZ", "5"))
os.makedirs(d, exist_ok=True)

if not os.environ.get("NOGRID"):
    subprocess.run([sys.executable, os.path.join(T, "nidelva_harbour_grid.py"), geo, d, "--q", str(Q), "--wl", str(WL),
                    "--nz", str(nz)], check=True)
control = open(os.path.join(d, "control.txt")).read().replace("M 10 NP", "M 10 %d" % np_)
open(os.path.join(d, "control.txt"), "w").write(control)

if solver == "sflow":
    flow = "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nP 181 -1\nP 182 %g\n" % out
else:
    flow = "A 10 5\nA 520 %s\nP 10 1\nP 30 %g\n" % (os.environ.get("A520", "0"), out)
nhb = json.load(open(os.path.join(d, "grid.json"))).get("harbour_blocks", 1)
ctrl = (flow + "B 50 0.05\nB 60 1\nW 10 %g\nF 50 2\nF 60 %g\nF 62 %g\nG 50 1\nN 41 %g\nN 47 0.3\nM 10 %d\nW 22 -9.81\n"
        % (Q, WL, WL, tend, np_) + "".join("G 54 %d 2 %g\n" % (b, WL) for b in range(1, nhb+1)))
if os.environ.get("AMR"):
    # one level over the lower river, the river mouth and the east end of the canal (AMRBOX: xs xe ys ye)
    ctrl += "G 1 1\nG 2 0\nG 10 %s\n" % os.environ.get("AMRBOX", "1100 1950 -900 1500")
if os.environ.get("SED"):
    # bed load (d50 0.5 mm), the morphological factor MF, the bed slope and the spiral flow of the bends (P5b)
    ctrl += ("S 10 1\nS 11 2\nS 16 2\nS 20 0.0005\nS 21 100\nS 15 3\nS 17 1.0\nS 35 %s\nS 41 1\nS 43 4000\nS 44 1\nS 31 1\n"
             "S 100 0\nS 71 2.0\nS 86 1\nS 95 1\n" % os.environ.get("MF", "50"))
ctrl += os.environ.get("XCTRL", "").replace(";", "\n") + "\n"
open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)

with open(os.path.join(d, "dive.log"), "w") as log:
    subprocess.run([DIVEMESH], cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
mpi = os.environ.get("REEF3D_MPIRUN", "mpirun --allow-run-as-root --oversubscribe").split()
cmd = [REEF3D] if np_ == 1 else mpi + ["-np", str(np_), REEF3D]
with open(os.path.join(d, "reef.log"), "w") as log:
    r = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT)
print("rc", r.returncode, d)
