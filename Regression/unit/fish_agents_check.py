#!/usr/bin/env python3
# The agent fish (tools/fish_agents/fish_agents.py: the rule brain of 6DOF_brain.cpp in its Python form driving a
# cheap body), no REEF3D run needed.
#   * a uniform current of 0.2 m/s: 12 fish on stations 0.6 m apart, each starting 1 m downstream facing the
#     current, hold their stations (mean distance after 12 s below 0.15 m, every fish within 0.5 m)
#   * the same with the inflow and outflow sides seen as walls: the fish at 1 m from the inflow cannot hold
#     its station (the vision of the brain turns it away; the reason for "walls")
#   * two fish sent to the same station: their bodies never come closer than 2 r_fish (within 1 %)
#   * a flow from a REEF3D output (a VTU written here, a uniform 0.2 m/s on a 9 x 7 grid): the same tracks as
#     the uniform flow
#   python3 fish_agents_check.py        exit 1 on a failure
# Architect: Hans Bihs
import json, math, os, struct, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "tools", "fish_agents"))
import numpy as np                 # noqa: E402
import fish_agents as fa           # noqa: E402

nfail = npass = 0


def check(ok, what):
    global nfail, npass
    print(("  ok    " if ok else "  FAIL  ") + what)
    if ok:
        npass += 1
    else:
        nfail += 1


def stations(walls, flow):
    fish = [{"x": x + 1.0, "z": z, "h": 180.0, "target": [x, z]} for x in (1.0, 2.2, 3.4) for z in (1.2, 1.8, 2.4, 3.0)]
    return {"T": 20.0, "dt": 0.01, "out_dt": 0.1, "body": {}, "brain": {"U_des": 0.5, "b_max": 1.0, "D": 1.5},
            "flow": flow, "box": [0, 6, 0, 4.2], "walls": walls, "cylinders": [], "fish": fish}


def dists(fish, rows, t0):
    s = fa.summarise(fish, rows, t0)
    return [x["mean"] for x in s], [x["max"] for x in s]


# 1. a uniform current, the stations held
case = stations([0, 0, 1, 1], {"type": "uniform", "u": [0.2, 0.0]})
fish, rows = fa.run(case)
md, mx = dists(fish, rows, 12.0)
check(np.mean(md) < 0.15 and max(mx) < 0.5, "current 0.2, 12 fish: mean distance after 12 s %.3f m, largest %.3f m" % (np.mean(md), max(mx)))

# 2. the inflow seen as a wall: the column 1 m from it does not hold
fish2, rows2 = fa.run(stations([1, 1, 1, 1], {"type": "uniform", "u": [0.2, 0.0]}))
md2, _ = dists(fish2, rows2, 12.0)
check(np.mean(md2[:4]) > 0.5 and np.mean(md2[8:]) < 0.15,
      "inflow seen as a wall: the fish 1 m from it %.2f m off, those 3.4 m from it %.2f m" % (np.mean(md2[:4]), np.mean(md2[8:])))

# 3. two fish on one station: the bodies do not overlap
case3 = {"T": 15.0, "dt": 0.01, "out_dt": 0.01, "body": {}, "brain": {"U_des": 0.5, "b_max": 1.0, "D": 1.5},
         "flow": {"type": "uniform", "u": [0.2, 0.0]}, "box": [0, 6, 0, 4], "walls": [0, 0, 1, 1], "cylinders": [],
         "fish": [{"x": 3.0, "z": 1.7, "h": 180.0, "target": [2.0, 2.0]}, {"x": 3.0, "z": 2.3, "h": 180.0, "target": [2.0, 2.0]}]}
fish3, rows3 = fa.run(case3)
R = np.array(rows3)
a, b = R[R[:, 1] == 0], R[R[:, 1] == 1]
dmin = float(np.min(np.hypot(a[:, 2] - b[:, 2], a[:, 3] - b[:, 3])))
check(dmin > 0.99 * 2 * 0.15, "two fish on one station: closest %.3f m (2 r_fish = 0.30 m)" % dmin)


# 4. the flow from a REEF3D output: a uniform 0.2 m/s written as a VTU (binary appended, as REEF3D writes it)
def write_vtu(path, t, xs, zs, u):
    P = np.array([[x, 0.0, z] for z in zs for x in xs], "<f4")
    V = np.array([[u, 0.0, 0.0]] * len(P), "<f4")
    head = ('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n'
            '<UnstructuredGrid>\n<FieldData>\n<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1"> %g\n'
            '</DataArray>\n</FieldData>\n<Piece NumberOfPoints="%d" NumberOfCells="0">\n<Points>\n'
            '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="0"/>\n</Points>\n<PointData>\n'
            '<DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="appended" offset="%d"/>\n'
            '</PointData>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding="raw">\n_' % (t, len(P), 4 + P.nbytes))
    with open(path, "wb") as f:
        f.write(head.encode())
        for A in (P, V):
            f.write(struct.pack("<I", A.nbytes)); f.write(A.tobytes())
        f.write(b"\n</AppendedData>\n</VTKFile>\n")


with tempfile.TemporaryDirectory() as d:
    os.makedirs(os.path.join(d, "REEF3D_CFD_VTU"))
    for k in range(2):
        write_vtu(os.path.join(d, "REEF3D_CFD_VTU", "REEF3D-CFD-%08d-000001.vtu" % k), 10.0 + k, np.linspace(0, 6, 9), np.linspace(0, 4.2, 7), 0.2)
    case4 = stations([0, 0, 1, 1], {"type": "reef3d", "run": d, "t0": 10.0, "t1": 11.0})
    fish4, rows4 = fa.run(case4)
diff = float(np.max(np.abs(np.array(rows4, float)[:, 2:4] - np.array(rows, float)[:, 2:4])))
check(diff < 1e-5, "REEF3D flow (uniform 0.2 m/s in a VTU) against the uniform flow: largest difference of the tracks %.2e m" % diff)

print("%d passed, %d failed" % (npass, nfail))
sys.exit(1 if nfail else 0)
