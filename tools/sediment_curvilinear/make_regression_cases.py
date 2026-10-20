# writes the regression cases of the curvilinear sediment transport (P5a, S 10 1 with G 50): one block (G 50 2,
# curvgrid.dat), mesh refinement, multi-block grids (G 50 1, blocks.dat): Regression/cases/*_sed_*; the bend terms
# of P5b (S 86 transverse bed slope, S 95 spiral flow): *_sed_*bend*; the lags of P7 (S 95 2 the spiral flow with its
# adaptation length, S 64 the suspended part with the adaptation length of the suspension): *_sed_*lags*
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

# P5b: 8 the bed-load dam break of t1 turned by 90 deg with SFLOW, the transverse slope and the spiral flow, S 31 2
nx, ny, nz, dx = 32, 24, 6, 0.125
X, Y = nodes(nx, ny, dx, "rot90", 1000.0, 2000.0)
b = box_phys((1.0, 1.75, 1.0, 2.0), "rot90", 1000.0, 2000.0)
bend = "S 86 1\nS 95 1\n"
case("sflow_2d_curv_sed_bend_dambreak_rot90", "sflow", control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=1),
     flow["sflow"] + "N 47 0.3\nF 60 0.2\nM 10 1\nP 10 0\nW 22 -9.81\n" + bedload.replace("S 31 1", "S 31 2") + "S 35 20.0\n" + bend
     + "F 72 %g %g %g %g 0.35\nG 50 2\n" % b,
     {"np": 1, "steps": 40, "tags": ["2d", "sflow", "quick", "bend"] + tags,
      "description": "The bend terms of the bed load (P5b) with SFLOW on a curvilinear grid turned by 90 deg: the transverse bed slope effect (S 86 1, Talmon et al. 1995) and the spiral-flow correction of the transport direction (S 95 1), the cell flux vector of S 31 2; the bed-load dam break of t1, per cell index the Cartesian run to round-off (sediment_curvilinear/t1).",
      "covers": ["sediment_exner_bend: bend_prep (qb/f(theta), streamline curvature, tan(delta)), bend_flux, bend_div", "the turned face directions (sediment_curv::face_sign with td) and cell vectors (topovel_curv form 2)"]},
     G=(X, Y, np.zeros_like(X)))

# 9 the LFM bend of t6 (coarse, short): SFLOW, 180 deg, the bed load fed from the fixed inflow reach
R, B, h0, Q, d50, Cz, sl, ds = 4.25, 1.7, 0.2, 0.17, 0.00078, 26.4, 0.0018, 0.2
n1, n2, nb, nyb = int(round(11.0/ds)), int(round(6.0/ds)), int(round(np.pi*R/ds)), 9
s_n = np.concatenate([np.linspace(0, 11.0, n1+1), 11.0 + np.linspace(0, np.pi*R, nb+1)[1:], 11.0 + np.pi*R + np.linspace(0, 6.0, n2+1)[1:]])
n_n = np.linspace(-0.5*B, 0.5*B, nyb+1)
S, N = np.meshgrid(s_n, n_n, indexing="ij")
X, Y = S.copy(), N.copy()
mb_ = (S > 11.0) & (S < 11.0 + np.pi*R)
m2 = S >= 11.0 + np.pi*R
phi = (S[mb_] - 11.0)/R - 0.5*np.pi
X[mb_], Y[mb_] = 11.0 + (R - N[mb_])*np.cos(phi), R + (R - N[mb_])*np.sin(phi)
X[m2], Y[m2] = 11.0 - (S[m2] - 11.0 - np.pi*R), 2*R - N[m2]
ks = 12.0*h0/np.exp(Cz/np.sqrt(9.81)/2.5)
nxb = len(s_n) - 1
X9, Y9, Z9 = X, Y, sl*(s_n[-1] - S)
ctrl9 = ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nB 50 %g\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\n" % (ks, Q, h0, h0)
     + "S 10 1\nS 11 2\nS 16 2\nS 20 %g\nS 21 %g\nS 15 3\nS 17 1.0\nS 35 100\nS 41 1\nS 43 200\nS 44 1\nS 31 1\nS 100 0\nS 71 2.0\n" % (d50, ks/d50)
     + bend + "G 50 2\n")
case("sflow_2d_curv_sed_bend_lfm", "sflow", control(nxb, nyb, 1, nxb*ds, nyb*ds, 0.5, bc=bc, np_=1), ctrl9,
     {"np": 1, "steps": 400, "tags": ["2d", "sflow", "river", "bend"] + tags,
      "description": "River bend with bed load (P5b): the 180 deg LFM flume of Delft Hydraulics (R 4.25 m, B 1.7 m, h 0.2 m, Q 0.17 m3/s, d50 0.78 mm) on a coarse curvilinear grid, SFLOW, Keulegan bed shear, the transverse slope (S 86 1) and the spiral flow (S 95 1), the first 2 m of the inflow reach fixed (S 71) feed the bed load; the start of the point bar and the outer-bank scour (sediment_curvilinear/t6).",
      "covers": ["S 95 with a curved mean flow: streamline curvature from the cell velocities, alpha_I from the bed shear", "the slope flux at the walls and next to the fixed inflow reach"]},
     G=(X, Y, sl*(s_n[-1] - S)))

# P7: 12 the LFM bend of case 9 with the lags: the spiral flow builds up over L_a (S 95 2), 75 % of the load suspended
# with its adaptation length (S 64 0.75)
G9 = (X9, Y9, Z9)
case("sflow_2d_curv_sed_bend_lfm_lags", "sflow", control(nxb, nyb, 1, nxb*ds, nyb*ds, 0.5, bc=bc, np_=1),
     ctrl9.replace(bend, "S 86 1\nS 95 2\nS 64 0.75\nS 65 4\n"),
     {"np": 1, "steps": 400, "tags": ["2d", "sflow", "river", "bend", "lags"] + tags,
      "description": "The lags of the bend morphology (P7) in the LFM bend of sflow_2d_curv_sed_bend_lfm: the spiral-flow intensity relaxed along the flow over its adaptation length L_a = (m-1) h/(2 kappa^2) (S 95 2, Kalkwijk & Booij 1986), 75 % of the transport suspended and relaxed over L_s = h U/(S65 w_s) (S 64 0.75, S 65 4); the spiral flow and the slope turn the total load (sediment_curvilinear/t7, t8).",
      "covers": ["sediment_exner_lag: relax_along (upwind along the flow, contravariant on the curvilinear grid, Gauss-Seidel with alternating sweeps), susp_split, spiral_lag, spiral_apply", "qb/f of the slope flux with the lagged suspended load (bend_K)"]},
     G=G9)

# 10 the four turned quadrants of t5 with SFLOW on three ranks and the bend terms (blocks on other ranks than the cells
# next to them: the bend fields of the cells around the blocks from their sources)
n, nz, dx = 40, 6, 0.05
L = n*dx
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.15*y + 0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.05*x + 0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.06*np.exp(-((X-1.3)**2+(Y-0.7)**2)/0.04)
quads = [quad(0, 0, 0), quad(hq, 0, 1), quad(0, hq, 2), quad(hq, hq, 3)]
case("sflow_2d_mb_sed_bend_mpi3", "sflow", mbcontrol.replace("NP", "3"),
     flow["sflow"] + mb.replace("NP", "3").replace("S 31 1", "S 31 2") + bend,
     {"np": 3, "steps": 40, "tags": ["2d", "sflow", "multiblock", "mpi", "bend"] + tags,
      "description": "The bend terms (S 86 1, S 95 1, S 31 2) on the multi-block grid of t5 with SFLOW on three ranks: the same bed as one block to round-off, the bed volume conserved (sediment_curvilinear/t5).",
      "covers": ["the bend fields qb/f, tan(delta), the bed and the open flag of the slope flux filled around the blocks (sediment_exner::amr_field)"]},
     blocks=quads)

# 13 the four turned quadrants of t5 (SFLOW) with a static refinement box over the vertex and the lags: the relaxations
# hold the cells around every grid and the covered cells, passes until the interfaces no longer move them
lags = "S 95 2\nS 97 0.1\nS 64 0.5\nS 86 1\n"
case("sflow_2d_mb_sed_lags_amr", "sflow", mbcontrol.replace("NP", "1"),
     flow["sflow"] + mb.replace("NP", "1") + lags + "G 1 1\nG 2 0\nG 3 2\nG 4 8\nG 10 0.7 1.4 0.7 1.4\n",
     {"np": 1, "steps": 40, "tags": ["2d", "sflow", "multiblock", "amr", "bend", "lags"] + tags,
      "description": "The lags of P7 (S 95 2 with a short adaptation length S 97 0.1, S 64 0.5) with the transverse slope (S 86 1) on the multi-block grid of t5 with a refinement box over the vertex: the relaxations of every grid hold the cells around it and its covered cells, passes across the grids until the interfaces no longer move them; the same bed as one block to 1e-4 of the change (the order of the Gauss-Seidel sweeps differs, the dam break amplifies it at its wet front), the bed volume conserved (sediment_curvilinear/t5).",
      "covers": ["sflow_amr::sediment_step_sflow: lagfix, the passes of lag_refresh with the fills (sync), bend_K after the passes", "the lagged fields Isp and qsn filled around the grids (amr_field 5, 6), prolonged and restricted"]},
     blocks=quads)
case("sflow_2d_mb_sed_lags_mpi2", "sflow", mbcontrol.replace("NP", "2"),
     flow["sflow"] + mb.replace("NP", "2") + lags,
     {"np": 2, "steps": 40, "tags": ["2d", "sflow", "multiblock", "mpi", "bend", "lags"] + tags,
      "description": "As sflow_2d_mb_sed_lags_amr without refinement on two ranks: the lagged fields of the blocks on other ranks than the cells next to them (remote fills), the relaxation of the native grid across the ranks (sediment_curvilinear/t5).",
      "covers": ["remote fills of Isp and qsn", "relax_along with the exchange of the ranks, the global change of the passes"]},
     blocks=quads)

# 11 NHFLOW, two levels of patches following the surface, the transverse slope: the slope flux of two patches that meet
nx, ny, nz, dx = 32, 24, 6, 0.125
X, Y = nodes(nx, ny, dx, "rot90", 1000.0, 2000.0)
case("nhflow_3d_curv_sed_bend_amr_adaptive", "nhflow", control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=1),
     flow["nhflow"] + amrbase.replace("NP", "1") + "G 1 2\nG 2 4\nG 20 0.01\nS 86 1\n",
     {"np": 1, "steps": 40, "tags": ["3d", "nhflow", "amr", "bend"] + tags,
      "description": "The transverse bed slope effect (S 86 1) on the mesh refinement of NHFLOW, two levels of patches following the surface on a curvilinear grid turned by 90 deg: patches that meet compute the same slope flux (the bed of the cells around them from their source, not the flow bed of a fresh patch); the bed volume over the leaf cells is conserved (sediment_curvilinear/t4).",
      "covers": ["amr_qb with bend_prep, the fill of the bend fields (nhflow_amr_sed: block_up, fill_run)"]},
     G=(X, Y, np.zeros_like(X)))
print("written")
