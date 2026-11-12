# P8c: one Weber junction run: weber_case.py writes the case, DIVEMesh, REEF3D with the AMR dump (REEF3D_AMR_DUMP) that
# weber_compare.py, weber_fig8.py and weber_fig910.py read.  Arguments after the case folder go to weber_case.py.
# usage: python3 weber_run.py case [weber_case.py options]   (environment: REEF3D_BIN, DIVEMESH_BIN, REEF3D_MPIRUN)
import os, sys, subprocess, shutil, time
T = os.path.dirname(os.path.abspath(__file__))
case = os.path.abspath(sys.argv[1])
args = sys.argv[2:]
np_ = int(args[args.index("--np")+1]) if "--np" in args else 2
if os.path.exists(case):
    shutil.rmtree(case)
subprocess.run([sys.executable, os.path.join(T, "weber_case.py"), case] + args, check=True)
DM = os.environ.get("DIVEMESH_BIN", shutil.which("DiveMESH") or "DiveMESH")
R3 = os.environ.get("REEF3D_BIN", os.path.join(T, "../../bin/REEF3D"))
with open(os.path.join(case, "dive.log"), "w") as log:
    subprocess.run([DM], cwd=case, stdout=log, stderr=subprocess.STDOUT, check=True)
env = dict(os.environ)
env["REEF3D_AMR_DUMP"] = os.path.join(case, "amr")
os.makedirs(env["REEF3D_AMR_DUMP"], exist_ok=True)
mpi = os.environ.get("REEF3D_MPIRUN", "mpirun --allow-run-as-root --oversubscribe").split()
cmd = [R3] if np_ == 1 else mpi + ["-np", str(np_), R3]
t0 = time.time()
with open(os.path.join(case, "reef.log"), "w") as log:
    r = subprocess.run(cmd, cwd=case, stdout=log, stderr=subprocess.STDOUT, env=env)
print("rc %d  %.0f s wall  %s" % (r.returncode, time.time()-t0, case))
