# writes the regression cases of the curvilinear sediment transport (P5a, S 10 1 with G 50): one block (G 50 2,
# curvgrid.dat), mesh refinement, multi-block grids (G 50 1, blocks.dat): Regression/cases/*_sed_*
import os, sys, json, numpy as np
S5 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, S5)
sys.path.insert(0, os.path.join(S5, "../sflow_multiblock"))
from p5lib import control, nodes, box_phys
from p2lib import write_blocks
C = os.path.join(S5, "../../Regression/cases")

def grid(path, X, Y, Z):
    ni, nj = X.shape[0]-1, X.shape[1]-1
    with open(path, "w") as f:
        f.write("%d %d\n" % (ni, nj))
        for j in range(nj+1):
            for i in range(ni+1):
                f.write("%.17g %.17g %.17g\n" % (X[i, j], Y[i, j], Z[i, j]))

def case(name, solver, control_txt, ctrl_txt, meta, G=None, blocks=None):
    d = os.path.join(C, name)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "control.txt"), "w").write(control_txt)
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl_txt)
    meta = dict(meta)
    meta["solver"] = solver
    if G is not None:
        grid(os.path.join(d, "curvgrid.dat"), *G)
        meta["files"] = ["curvgrid.dat"]
    if blocks is not None:
        write_blocks(os.path.join(d, "blocks.dat"), blocks)
        meta["files"] = ["blocks.dat"]
    json.dump(meta, open(os.path.join(d, "case.json"), "w"), indent=1)

flow = {"nhflow": "A 10 5\nA 520 0\n", "sflow": "A 10 2\nA 209 0\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"}
bedload = "S 10 1\nS 11 2\nS 20 0.0005\nS 15 3\nS 17 1.0\nS 41 1\nS 43 1\nS 44 1\nS 31 1\n"
tags = ["curvilinear", "sediment"]

# 1 the bed-load dam break of t1 on the grid turned by 90 deg (NHFLOW, hydrostatic): per cell index the Cartesian run
nx, ny, nz, dx = 32, 24, 6, 0.125
X, Y = nodes(nx, ny, dx, "rot90", 1000.0, 2000.0)
b = box_phys((1.0, 1.75, 1.0, 2.0), "rot90", 1000.0, 2000.0)
case("nhflow_3d_curv_sed_dambreak_rot90", "nhflow", control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=1),
     flow["nhflow"] + "N 47 0.3\nN 48 0\nN 49 0.01\nF 60 0.2\nM 10 1\nP 10 0\nW 22 -9.81\n" + bedload + "S 35 20.0\n"
     + "F 72 %g %g %g %g 0.35\nG 50 2\n" % b,
     {"np": 1, "steps": 40, "tags": ["3d", "nhflow", "quick"] + tags,
      "description": "Sediment transport (S 10 1, bed load, S 31 1) with NHFLOW on a curvilinear grid (G 50 2): bed-load dam break over a flat erodible bed, the grid the Cartesian one turned by 90 deg and shifted; per cell index the Cartesian run to round-off (sediment_curvilinear/t1).",
      "covers": ["sediment_curv: face normals P, Q, contravariant cell velocity, face sign, chain-rule bed gradient", "sediment_exner_curv (topovel_curv), the symmetric filter of the curvilinear grid (sediment_xfilter)", "bedslope, bedload direction, bed shear on the curvilinear grid; the curvilinear time step of NHFLOW"]},
     G=(X, Y, np.zeros_like(X)))

# 2, 3 the river of t2 with the steep mound and the sand slide: SFLOW on the wavy grid, NHFLOW turned by 30 deg at UTM
nx, ny, nz, dx = 40, 12, 6, 0.1
L, W, h = nx*dx, ny*dx, 0.25
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
Z = np.maximum(0.08 - np.hypot(I*dx-1.5, J*dx-0.5*W)*0.15, 0.0)
bc = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n"
river = ("B 50 0.01\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nN 47 0.3\nM 10 NP\nP 10 0\nW 22 -9.81\n" % (0.4*W*h, h, h)
         + bedload + "S 35 50.0\nS 90 1\nS 81 8\nS 82 0\nG 50 2\n")
X, Y = nodes(nx, ny, dx, "id", wavy=0.25)
case("sflow_2d_curv_sed_river_wavy", "sflow", control(nx, ny, nz, L, W, 0.5, bc=bc, np_=1),
     flow["sflow"] + river.replace("NP", "1"),
     {"np": 1, "steps": 60, "tags": ["2d", "sflow", "river"] + tags,
      "description": "Sediment transport with SFLOW on a skewed (wavy) curvilinear grid: river reach with discharge inflow, stage outflow and Manning friction over a sand mound, bed load and the sand slide (S 90 1, angle of repose 8 deg); the bed volume changes by the sediment flux through the in- and outflow only (sediment_curvilinear/t2).",
      "covers": ["fill_PQ_sflow on the curvilinear grid, bedshear taubed_curv", "sandslide_f::slide_curv (distance and area ratio of the cells)"]},
     G=(X, Y, Z))
X, Y = nodes(nx, ny, dx, "rot", 2.7e5, 7.04e6, 30.0)
case("nhflow_3d_curv_sed_river_rot30_mpi2", "nhflow", control(nx, ny, nz, L, W, 0.5, bc=bc, np_=2),
     flow["nhflow"] + "A 519 1\nN 48 0\nN 49 0.01\n" + river.replace("NP", "2"),
     {"np": 2, "steps": 60, "tags": ["3d", "nhflow", "river", "mpi"] + tags,
      "description": "As sflow_2d_curv_sed_river_wavy with NHFLOW (hydrostatic, fixed time step) on the grid turned by 30 deg at UTM coordinates, two ranks: the same bed as one rank (sediment_curvilinear/t3).",
      "covers": ["sediment_curv across the partition (P, Q, the bed gradient in the halo)", "sand slide on several ranks with the curvilinear distances"]},
     G=(X, Y, Z))

# 4, 5 mesh refinement (t4): the dam break of t1 turned by 90 deg, NHFLOW with patches following the surface and the slide,
# SFLOW with two levels in a static box on two ranks
nx, ny, nz, dx = 32, 24, 6, 0.125
X, Y = nodes(nx, ny, dx, "rot90", 1000.0, 2000.0)
b = box_phys((1.0, 1.75, 1.0, 2.0), "rot90", 1000.0, 2000.0)
rb = box_phys((0.5, 2.5, 0.75, 2.25), "rot90", 1000.0, 2000.0)
amrbase = "N 47 0.3\nN 48 0\nN 49 0.01\nF 60 0.2\nM 10 NP\nP 10 0\nW 22 -9.81\n" + bedload + "S 35 20.0\nG 3 2\nG 4 8\nG 50 2\n" + "F 72 %g %g %g %g 0.35\n" % b
case("nhflow_3d_curv_sed_amr_adaptive", "nhflow", control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=1),
     flow["nhflow"] + amrbase.replace("NP", "1") + "G 1 1\nG 2 4\nG 20 0.01\nS 90 1\nS 81 8\nS 82 0\n",
     {"np": 1, "steps": 40, "tags": ["3d", "nhflow", "amr"] + tags,
      "description": "Sediment transport on the mesh refinement of NHFLOW (G 1 1, patches following the surface, regrid every 4 steps) on a curvilinear grid turned by 90 deg: the Exner fluxes of the coarse faces next to a patch from the patch, one dtsed, the filter and the sand slide over all grids, the bed of fresh patches from the old patches and their parents; the bed volume over the leaf cells is conserved (sediment_curvilinear/t4).",
      "covers": ["nhflow_amr_sed: the sediment step of all grids, flux matching, qb restriction and fill, sed_bed_ghosts", "reefamr_sed.h: the filter (sediment_xfilter) and the sand slide (sandslide amr_sweep) across the grids", "the water level of the cells around the patches after a bed update (fill_stage with S 10)"]},
     G=(X, Y, np.zeros_like(X)))
case("sflow_2d_curv_sed_amr_lev2_mpi2", "sflow", control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=2),
     flow["sflow"] + amrbase.replace("NP", "2") + "G 1 2\nG 2 0\nG 10 %g %g %g %g\n" % rb,
     {"np": 2, "steps": 40, "tags": ["2d", "sflow", "amr", "mpi"] + tags,
      "description": "Sediment transport on the mesh refinement of SFLOW: two levels in a static box on a curvilinear grid turned by 90 deg, two ranks (sediment_curvilinear/t4).",
      "covers": ["sflow_amr_sed on several ranks: remote fills of qb and of the bed, remote flux matching (face_run_at), the cross-grid exchanges of the filter"]},
     G=(X, Y, np.zeros_like(X)))

# 6, 7 multi-block grids (t5): the four turned quadrants of P4 around a regular vertex, a steep bump, the sand slide
n, nz, dx = 40, 6, 0.05
L = n*dx
hq = n//2
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.15*y + 0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.05*x + 0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.2*np.exp(-((X-1.3)**2+(Y-0.7)**2)/0.04)
def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+hq+1, j0:j0+hq+1], k=-k).copy() for A in (X, Y, Z))
quads = [quad(0, 0, 0), quad(hq, 0, 1), quad(0, hq, 2), quad(hq, hq, 3)]
mbcontrol = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, L, nz*dx)
mb = ("F 60 0.2\nF 72 0.65 1.45 0.65 1.45 0.38\nG 50 1\nG 52 1\nA 509 0\nN 47 0.3\nN 48 0\nN 49 0.005\nM 10 NP\nP 10 0\nW 22 -9.81\n"
      + bedload + "S 35 20.0\n")
case("nhflow_3d_mb_sed_amr", "nhflow", mbcontrol.replace("NP", "1"),
     flow["nhflow"] + mb.replace("NP", "1") + "S 90 1\nG 1 1\nG 2 0\nG 3 2\nG 4 8\nG 10 0.7 1.4 0.7 1.4\n",
     {"np": 1, "steps": 20, "tags": ["3d", "nhflow", "multiblock", "amr"] + tags,
      "description": "Sediment transport on a multi-block grid (G 50 1, DIVEMesh R 20 1) with NHFLOW: four curvilinear blocks turned by 0, 90, 180, 270 deg around a regular vertex, a steep bump, bed load, the filter and the sand slide, one refinement level in a static box over the vertex: the same bed as one block to round-off (sediment_curvilinear/t5).",
      "covers": ["the sediment of the blocks and of the ring (fill_run0), the interfaces of the blocks and of the patches across them open for the filter (the intermediate fields swapped across turned blocks) and the sand slide (xslice_add)", "the water level of the cells around the patches after a bed update"]},
     blocks=quads)
case("sflow_2d_mb_sed_mpi2", "sflow", mbcontrol.replace("NP", "2"),
     flow["sflow"] + mb.replace("NP", "2") + "S 90 4\n",
     {"np": 2, "steps": 20, "tags": ["2d", "sflow", "multiblock", "mpi"] + tags,
      "description": "As nhflow_3d_mb_sed_amr with SFLOW, without refinement, the weighted multi-direction sand slide (S 90 4), two ranks (blocks on other ranks than the native cells next to them): the same bed as one block to round-off (sediment_curvilinear/t5).",
      "covers": ["remote fills and the reverse exchange (xrun_rev) of the slid material", "sandslide_weighted_multidir in steps (amr_sweep)"]},
     blocks=quads)
print("written")
