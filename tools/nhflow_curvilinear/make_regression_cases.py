# writes the curvilinear NHFLOW regression cases (G 50 2, curvgrid.dat): Regression/cases/nhflow_3d_curv_*
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
    meta = dict(meta); meta["solver"] = "nhflow"; meta["files"] = ["curvgrid.dat"]
    json.dump(meta, open(os.path.join(d, "case.json"), "w"), indent=1)

walls = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\n"

# 1 dam break in a closed basin, the grid turned by 90 deg and shifted (A 520 2)
n, nz, dx = 40, 8, 0.05
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z0 = np.zeros_like(X0)
basin = walls + "B 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 %%d\n" % (n, n, nz, n*dx, n*dx)
case("nhflow_3d_curv_dambreak_rot90", basin % 1,
     "A 10 5\nA 520 2\nF 60 0.5\nF 72 999.89 1000.09 1999.91 2000.11 0.54\nG 50 2\nN 47 0.5\nM 10 1\nP 10 0\nW 22 -9.81\n",
     (1000.0-(Y0-1.0), 2000.0+(X0-1.0), Z0),
     {"np": 1, "steps": 30, "tags": ["3d", "nhflow", "curvilinear", "quick"],
      "description": "NHFLOW on a curvilinear grid (G 50 2): dam break in a closed basin whose grid is the Cartesian one turned by 90 deg and shifted, non-hydrostatic (A 520 2): per cell index the Cartesian run to round-off.",
      "covers": ["nhflow_curv: normal velocities, HLL fluxes, walls (pair_bc, bc_fluxes), sigma metrics, Poisson rows, pressure gradients"]})

# 2 the same, turned by 30 deg at UTM coordinates, two ranks
th = np.radians(30.0)
case("nhflow_3d_curv_dambreak_rot30_mpi2", basin % 2,
     "A 10 5\nA 520 2\nF 60 0.5\nF 72 %.4f %.4f %.4f %.4f 0.54\nG 50 2\nN 47 0.5\nM 10 2\nP 10 0\nW 22 -9.81\n" % (2.7e5+0.8, 2.7e5+1.0, 7.04e6+1.2, 7.04e6+1.4),
     (2.7e5 + X0*np.cos(th) - Y0*np.sin(th), 7.04e6 + X0*np.sin(th) + Y0*np.cos(th), Z0),
     {"np": 2, "steps": 30, "tags": ["3d", "nhflow", "curvilinear", "mpi"],
      "description": "NHFLOW on a curvilinear grid turned by 30 deg at UTM coordinates, two ranks: a water box released in a closed basin (F 72 in physical coordinates).",
      "covers": ["nhflow_curv on two ranks (geometry halo, chain-rule gradients across the partition)"]})

# 3 lake at rest with an island on a wavy, skewed grid (wetting and drying)
n2, nz2, dx2 = 60, 6, 0.5
L2 = n2*dx2
I, J = np.meshgrid(np.arange(n2+1), np.arange(n2+1), indexing="ij")
x, y = I*dx2, J*dx2
X = x + 0.25*y + 0.8*np.sin(2*np.pi*x/L2)*np.sin(np.pi*y/L2)
Y = y + 0.1*x + 0.6*np.sin(np.pi*x/L2)*np.sin(4*np.pi*y/L2)
Zw = 0.6*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0)
case("nhflow_3d_curv_wavy_lake",
     walls + "B 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\n" % (n2, n2, nz2, L2, L2),
     "A 10 5\nA 520 2\nF 60 1.0\nG 50 2\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\n",
     (X, Y, Zw),
     {"np": 1, "steps": 40, "tags": ["3d", "nhflow", "curvilinear", "quick"],
      "description": "NHFLOW on a wavy, skewed curvilinear grid: lake at rest over a submerged bump (A 520 2), velocities and pressure stay at round-off.",
      "covers": ["nhflow_curv well balancing (bed-slope source, non-orthogonal Poisson rows)"]})

# 4 wave channel turned by 30 deg at UTM coordinates: Dirichlet generation, active absorption
nx, ny, nzc, dxc, h = 200, 4, 10, 0.1, 0.5
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
s, nn = I*dxc, J*dxc
X0c, Y0c = 2.7e5, 7.04e6
chan = "C 11 6\nC 12 3\nC 13 3\nC 14 7\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\nM 20 2\n" % (nx, ny, nzc, nx*dxc, ny*dxc, h)
case("nhflow_3d_curv_waves_rot30", chan,
     "A 10 5\nA 520 2\nB 90 1\nB 92 1\nB 93 0.02 1.0\nB 98 3\nB 99 3\nB 105 30.0 %.3f %.3f\nF 60 %g\nG 50 2\nN 47 0.5\nM 10 1\nP 10 0\nW 22 -9.81\n" % (X0c, Y0c, h),
     (X0c + s*np.cos(th) - nn*np.sin(th), Y0c + s*np.sin(th) + nn*np.cos(th), np.zeros_like(s)),
     {"np": 1, "steps": 150, "tags": ["3d", "nhflow", "curvilinear", "waves", "iowave"],
      "description": "NHFLOW wave channel on a curvilinear grid turned by 30 deg at UTM coordinates: linear waves (kh 2) by Dirichlet generation at x- (B 98 3) along the face normal, active absorption at x+ (B 99 3), wave direction B 105.",
      "covers": ["iowave_nhflow_dirichlet on face normals", "nhflow_active_beach on face normals", "position::pos_x/pos_y on curvilinear grids", "bc_fluxes at in- and outflow"]})

# 5 the same channel stretched along x (orthogonal mapped grid), relaxation zones by B 107 / B 108
xi = I/nx
sm = nx*dxc*(xi + 0.05*np.sin(2*np.pi*xi))
def P(sv, nv):
    return X0c + sv*np.cos(th) - nv*np.sin(th), Y0c + sv*np.sin(th) + nv*np.cos(th)
a, b = P(0.0, -1.0), P(0.0, ny*dxc+1.0)
c_, e = P(nx*dxc, -1.0), P(nx*dxc, ny*dxc+1.0)
case("nhflow_3d_curv_waves_relax_map30", chan,
     ("A 10 5\nA 520 2\nB 90 1\nB 92 1\nB 93 0.02 1.0\nB 98 2\nB 99 1\nB 105 30.0 %.3f %.3f\n" % (X0c, Y0c))
     + "B 108 %.4f %.4f %.4f %.4f 3.0\nB 107 %.4f %.4f %.4f %.4f 5.0\n" % (a[0], b[0], a[1], b[1], c_[0], e[0], c_[1], e[1])
     + "F 60 %g\nG 50 2\nN 47 0.5\nM 10 1\nP 10 0\nW 22 -9.81\n" % h,
     (*P(sm, nn), np.zeros_like(s)),
     {"np": 1, "steps": 150, "tags": ["3d", "nhflow", "curvilinear", "waves", "iowave"],
      "description": "NHFLOW wave channel on an orthogonal mapped grid (stretched along the channel, turned by 30 deg): relaxation zones given by B 107 / B 108 lines in physical coordinates.",
      "covers": ["iowave relaxation (B 98 2, B 99 1) on curvilinear grids", "the gate test of phase P3 (dispersion on a mapped grid)"]})

# 6 river reach turned by 30 deg: discharge inflow, stage outflow, k-epsilon, implicit diffusion
nx, ny, nzr, dxr, hr = 60, 12, 8, 0.5, 1.5
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
s, nn = I*dxr, J*dxr
case("nhflow_3d_curv_river_kepsilon_rot30",
     "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\nM 20 2\n" % (nx, ny, nzr, nx*dxr, ny*dxr, hr),
     "A 10 5\nA 512 2\nA 519 1\nA 520 2\nA 560 1\nB 50 0.01\nB 60 1\nF 60 1.5\nF 62 1.5\nG 50 2\nN 47 0.5\nM 10 1\nP 10 0\nW 2 0.001\nW 10 9.0\nW 22 -9.81\n",
     (X0c + s*np.cos(th) - nn*np.sin(th), Y0c + s*np.sin(th) + nn*np.cos(th), np.zeros_like(s)),
     {"np": 1, "steps": 100, "tags": ["3d", "nhflow", "curvilinear", "river", "turbulence"],
      "description": "NHFLOW river reach on a curvilinear grid turned by 30 deg at UTM coordinates: discharge inflow along the inflow face normal (B 60 1), stage outflow (F 62), bed friction, k-epsilon (A 560 1) with the implicit diffusion (A 512 2).",
      "covers": ["ioflow nhflow inflow / Qin / Qout on face normals", "pair_bc river outflow", "nhflow_idiff curvilinear rows", "nhflow_scalar_ifou face fluxes", "nhflow_strain chain rule"]})

# 7 oblique tide through a basin turned by 30 deg (Riemann and Flather edges)
nt, nzt, dxt = 40, 3, 75.0
I, J = np.meshgrid(np.arange(nt+1), np.arange(nt+1), indexing="ij")
s, nn = I*dxt, J*dxt
edges = "B 520 1 3 0\nB 521 1 1 0.0 0.0 0.0\nB 523 1 1\nB 520 3 3 0\nB 521 3 3 0.0 0.0 0.0\nB 523 3 1\nB 520 2 4 0\nB 521 2 2 0.0 0.0 0.0\nB 523 2 1\nB 520 4 4 0\nB 521 4 4 0.0 0.0 0.0\nB 523 4 1\n"
case("nhflow_3d_curv_tide_rot30",
     "C 11 2\nC 12 2\nC 13 2\nC 14 2\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\nM 20 2\n" % (nt, nt, nzt, nt*dxt, nt*dxt),
     "A 10 5\nB 90 1\nB 92 0\nB 510 1 1 60.0 600.0\nB 511 1 0.25 600.0 90.0\nB 515 1 %.3f %.3f 0.0\n" % (X0c, Y0c) + edges
     + "F 60 10.0\nG 50 2\nN 47 0.5\nM 10 1\nP 10 0\nW 22 -9.81\n",
     (X0c + s*np.cos(th) - nn*np.sin(th), Y0c + s*np.sin(th) + nn*np.cos(th), np.zeros_like(s)),
     {"np": 1, "steps": 150, "tags": ["3d", "nhflow", "curvilinear", "tide", "iowave"],
      "description": "NHFLOW on a curvilinear grid turned by 30 deg at UTM coordinates: an oblique tide (direction 60 deg in the physical frame) through a basin with Riemann edges at x- and y- and Flather edges at x+ and y+.",
      "covers": ["iowave open edges on face normals (nhflow_open_edges)", "bc_fluxes at open edges", "background at the physical face centres"]})

# 8 mesh refinement: dam break with patches following the surface, turned by 90 deg, subcycling
case("nhflow_3d_curv_amr_adaptive_sub", basin % 1,
     "A 10 5\nA 520 2\nF 60 0.5\nF 72 999.89 1000.09 1999.91 2000.11 0.54\nG 50 2\nN 47 0.5\nM 10 1\nP 10 0\nW 22 -9.81\n"
     + "G 1 1\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\nG 7 1\n",
     (1000.0-(Y0-1.0), 2000.0+(X0-1.0), Z0),
     {"np": 1, "steps": 30, "tags": ["3d", "nhflow", "curvilinear", "amr"],
      "description": "NHFLOW mesh refinement on a curvilinear grid (turned by 90 deg): patches following the surface of a dam break, regrid every step, subcycling (G 7 1): per cell index the Cartesian run to round-off.",
      "covers": ["nhflow_amr on curvilinear grids: patch geometry, fills, refluxing per face length (sub_reflux)"]})

# 9 mesh refinement: lake at rest over a bump on the wavy grid, a static box across the bump
n3, nz3, dx3 = 48, 6, 0.5
L3 = n3*dx3
I, J = np.meshgrid(np.arange(n3+1), np.arange(n3+1), indexing="ij")
x, y = I*dx3, J*dx3
X = x + 0.25*y + 0.8*np.sin(2*np.pi*x/L3)*np.sin(np.pi*y/L3)
Y = y + 0.1*x + 0.6*np.sin(np.pi*x/L3)*np.sin(4*np.pi*y/L3)
case("nhflow_3d_curv_amr_wavy_lake",
     walls + "B 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\n" % (n3, n3, nz3, L3, L3),
     "A 10 5\nA 520 2\nF 60 1.0\nG 50 2\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\nG 1 1\nG 10 9.0 19.0 8.0 18.0\nG 4 4\n",
     (X, Y, 0.6*np.exp(-((X-14.0)**2+(Y-13.0)**2)/10.0)),
     {"np": 1, "steps": 40, "tags": ["3d", "nhflow", "curvilinear", "amr", "quick"],
      "description": "NHFLOW mesh refinement on a wavy, skewed curvilinear grid: lake at rest over a bump under a static patch (G 10 box in physical coordinates, six patches): velocities and pressure at round-off on both levels.",
      "covers": ["nhflow_amr: restriction per area (cw4), bed prolongation with the area mean of the parent"]})

# 10 mesh refinement with wetting and drying in the patches (G 30 1), turned by 90 deg
n4, nz4, dx4 = 64, 6, 0.0625
I, J = np.meshgrid(np.arange(n4+1), np.arange(n4+1), indexing="ij")
X0, Y0 = I*dx4, J*dx4
case("nhflow_3d_curv_amr_shore_rot90",
     walls + "B 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 0.35\n\nM 10 1\n" % (n4, n4, nz4, n4*dx4, n4*dx4),
     "A 10 5\nA 520 2\nF 60 0.2\nF 72 999.5 1000.5 1999.0 1999.5 0.3\nG 50 2\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\nG 1 1\nG 4 8\nG 30 1\nG 10 999.0 1001.0 2000.25 2001.5\n",
     (1000.0-(Y0-2.0), 2000.0+(X0-2.0), np.clip((X0-2.0)/2.0*0.3, 0.0, None)),
     {"np": 1, "steps": 40, "tags": ["3d", "nhflow", "curvilinear", "amr", "wetdry"],
      "description": "NHFLOW mesh refinement with wetting and drying in the patches (G 30 1) on a curvilinear grid turned by 90 deg: water running up a plane beach under a static patch.",
      "covers": ["nhflow_amr G 30 on curvilinear grids (F 72 and G 10 boxes in physical coordinates)"]})
print("written")
