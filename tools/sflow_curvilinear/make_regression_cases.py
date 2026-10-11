# writes the curvilinear SFLOW regression cases (G 50 2, curvgrid.dat)
import os, json, numpy as np
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../Regression/cases")

def grid(path, X, Y, Z):
    ni, nj = X.shape[0]-1, X.shape[1]-1
    with open(path, "w") as f:
        f.write("%d %d\n" % (ni, nj))
        for j in range(nj+1):
            for i in range(ni+1):
                f.write("%.17g %.17g %.17g\n" % (X[i, j], Y[i, j], Z[i, j]))

def case(name, control, ctrl, G, meta):
    d = os.path.join(C, name)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "control.txt"), "w").write(control)
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)
    grid(os.path.join(d, "curvgrid.dat"), *G)
    meta = dict(meta); meta["solver"] = "sflow"; meta["files"] = ["curvgrid.dat"]
    json.dump(meta, open(os.path.join(d, "case.json"), "w"), indent=1)

base_ctrl = """A 10 2
A 210 3
A 211 4
A 220 0
A 240 1
"""

# 1 dam break on a grid rotated by 90 deg and shifted
n, dx = 60, 0.5
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
case("sflow_2d_curv_dambreak_rot90",
     "B 1 0.5\nB 10 0.0 30.0 0.0 30.0 0.0 0.5\n\nM 10 1\nM 20 2\n",
     base_ctrl + "\nF 60 1.0\nF 72 997.0 1003.0 1997.0 2003.0 2.0\n\nG 50 2\n\nN 47 0.3\n\nW 22 -9.81\n",
     (1000.0-(Y0-15.0), 2000.0+(X0-15.0), np.zeros_like(X0)),
     {"np": 1, "steps": 60, "tags": ["2d", "sflow", "curvilinear", "quick"],
      "description": "SFLOW on a curvilinear grid (G 50 2): dam break in a closed basin whose grid is the Cartesian one rotated by 90 deg and shifted; the same eta per cell index as the Cartesian run (F 72 in physical coordinates).",
      "covers": ["sflow_curv normal velocities, fluxes, walls, divergence, bed-slope source, CFL", "lexer::curv_nodes (G 50 2)"]})

# 2 lake at rest with an island on a wavy, skewed grid
L = n*dx
X = X0 + 0.25*Y0 + 0.8*np.sin(2*np.pi*X0/L)*np.sin(np.pi*Y0/L)
Y = Y0 + 0.1*X0 + 0.6*np.sin(np.pi*X0/L)*np.sin(4*np.pi*Y0/L)
rng = np.random.default_rng(1)
Z = 1.3*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)
case("sflow_2d_curv_lake_island",
     "B 1 0.5\nB 10 0.0 30.0 0.0 30.0 0.0 0.5\n\nM 10 1\nM 20 2\n",
     base_ctrl + "\nF 60 1.0\n\nG 50 2\n\nN 47 0.3\n\nW 22 -9.81\n",
     (X, Y, Z),
     {"np": 1, "steps": 60, "tags": ["2d", "sflow", "curvilinear", "wetdry", "quick"],
      "description": "SFLOW on a wavy, skewed curvilinear grid (G 50 2): lake at rest with an emerged island (wetting and drying), bed levels from the grid nodes.",
      "covers": ["sflow_curv well balancing (closure of the face vectors)", "lexer::geometry_bed"]})

# 3 river channel with an island, rotated by 30 deg, UTM-like coordinates
ni, nj = 100, 20
I, J = np.meshgrid(np.arange(ni+1), np.arange(nj+1), indexing="ij")
X0, Y0 = I*1.0, J*1.0
Zb = 2.3*np.exp(-((X0-40.0)**2 + (Y0-14.0)**2)/30.0)
a = np.deg2rad(30.0)
river_control = "C 11 1\nC 12 3\nC 13 3\nC 14 2\n\nB 1 1.0\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\nM 20 2\n"
case("sflow_2d_curv_channel_rot30",
     river_control % (ni, nj),
     base_ctrl + "\nB 60 1\n\nF 50 2\nF 60 2.0\nF 62 2.0\n\nG 50 2\n\nN 47 0.3\n\nW 10 20.0\nW 22 -9.81\n",
     (5.0e5 + X0*np.cos(a) - Y0*np.sin(a), 7.0e6 + X0*np.sin(a) + Y0*np.cos(a), Zb),
     {"np": 1, "steps": 100, "tags": ["2d", "sflow", "curvilinear", "wetdry", "river"],
      "description": "SFLOW on a curvilinear grid (G 50 2): straight river reach rotated by 30 deg at UTM coordinates, discharge inflow (B 60 1) along the inflow face normal, stage outflow, an emerging island.",
      "covers": ["sflow_curv in- and outflow fluxes", "inflow2D/Qin2D/Qout2D on face normals", "magnitude Froude limiter"]})

# 4 90 deg bend with Manning friction on two ranks
R, B, Lin, Lout = 50.0, 20.0, 40.0, 60.0
Lb = 0.5*np.pi*R
ni = int(round(Lin+Lb+Lout)); nj = 20
s = np.linspace(-Lin, Lb+Lout, ni+1); nrm = np.linspace(-0.5*B, 0.5*B, nj+1)
X = np.zeros((ni+1, nj+1)); Y = np.zeros_like(X)
for i, sv in enumerate(s):
    if sv <= 0.0: c, nl = (sv, 0.0), (0.0, 1.0)
    elif sv <= Lb: t = sv/R; c, nl = (R*np.sin(t), R-R*np.cos(t)), (-np.sin(t), np.cos(t))
    else: c, nl = (R, R+(sv-Lb)), (-1.0, 0.0)
    X[i, :] = c[0] + nrm*nl[0]; Y[i, :] = c[1] + nrm*nl[1]
case("sflow_2d_curv_bend_mpi2",
     river_control % (ni, nj),
     base_ctrl + "A 218 1\n\nB 50 0.05\nB 60 1\n\nF 50 2\nF 60 2.0\nF 62 2.0\n\nG 50 2\n\nN 47 0.3\n\nW 10 20.0\nW 22 -9.81\n",
     (X, Y, np.zeros_like(X)),
     {"np": 2, "steps": 100, "tags": ["2d", "sflow", "curvilinear", "river", "mpi"],
      "description": "SFLOW on a curvilinear grid (G 50 2): 90 deg river bend (radius 50 m, width 20 m) with straight reaches, Manning friction, discharge inflow and stage outflow, two ranks.",
      "covers": ["sflow_curv on several ranks (face fluxes exchanged)", "superelevation in a bend (tests/t3)"]})
print("ok")
