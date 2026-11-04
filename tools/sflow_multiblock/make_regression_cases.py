# writes the multi-block regression cases (P2) into Regression/cases
import os, json, numpy as np
from p2lib import write_blocks
C = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../Regression/cases")

def case(name, meta, control=None, ctrl=None, blocks=None):
    d = os.path.join(C, name)
    os.makedirs(d, exist_ok=True)
    if control is not None:
        open(os.path.join(d, "control.txt"), "w").write(control)
    if ctrl is not None:
        open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)
    if blocks is not None:
        write_blocks(os.path.join(d, "blocks.dat"), blocks)
    json.dump(meta, open(os.path.join(d, "case.json"), "w"), indent=1)

# four blocks around a regular vertex, indices turned by 0, 90, 180, 270 deg, wavy skewed grid
n, dx = 64, 0.5
L = n*dx
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.15*y + 0.5*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.05*x + 0.4*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.3*np.exp(-((X-20.0)**2+(Y-12.0)**2)/8.0)
h = n//2
def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))
quads = [quad(0, 0, 0), quad(h, 0, 1), quad(0, h, 2), quad(h, h, 3)]

control = "B 1 0.5\nB 10 0.0 32.0 0.0 32.0 0.0 0.5\n\nM 10 1\nM 20 2\n\nR 20 1\n"
ctrl = ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 14.0 20.0 13.0 19.0 2.0\n\n"
        "G 50 1\n\nN 47 0.3\n\nM 10 1\n\n"
        "P 51 16.0 15.0\nP 51 19.5 18.5\nP 51 22.0 14.0\nP 51 15.0 20.0\nP 51 20.0 21.0\n\nW 22 -9.81\n")
files = ["blocks.dat"]
cov0 = ["reefamr_block: blocks 1.. as level-0 patches, block_resolve (two hops at the vertex), fill_run0, ring of the native grid",
        "DIVEMesh R 20 1: interfaces, transforms, BLKS section; geo_mesh BLKS reader"]
case("sflow_2d_mb_quadrants",
     {"description": "SFLOW multi-block level 0 (G 50 1, DIVEMesh R 20 1): four curvilinear blocks around a regular vertex, their indices turned by 0, 90, 180 and 270 deg, dam break over the vertex, bumpy bed. The same solution as one block to round-off (p2tests/t2).",
      "solver": "sflow", "np": 1, "steps": 40, "files": files, "tags": ["2d", "sflow", "curvilinear", "multiblock"], "covers": cov0},
     control, ctrl, quads)
case("sflow_2d_mb_quadrants_mpi2",
     {"base": "sflow_2d_mb_quadrants", "np": 2, "description": "As sflow_2d_mb_quadrants on 2 ranks: blocks on other ranks than the native cells next to them (remote fills of the blocks and of the ring).",
      "tags": ["2d", "sflow", "curvilinear", "multiblock", "mpi"], "covers": ["reefamr::blocks_plans remote requests, exchange_level0 of the ring"]})
amr = ["G 1 2", "G 2 0", "G 3 2", "G 4 8", "G 10 13.0 22.0 12.0 21.0"]
case("sflow_2d_mb_amr",
     {"base": "sflow_2d_mb_quadrants", "ctrl_add": amr,
      "description": "As sflow_2d_mb_quadrants with two refinement levels in a static box over the vertex: patches in all four blocks, cut at the block faces, filled and flux-matched across the turned interfaces. The same solution as one block with the same box to round-off (p2tests/t4).",
      "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr"],
      "covers": ["canonical fills across block faces, block_face_map (turned flux matching, sgn), tile maps over the block index space, blocks_propagate, bed_at/area_at of the blocks"]})
case("sflow_2d_mb_amr_mpi2",
     {"base": "sflow_2d_mb_amr", "np": 2, "description": "As sflow_2d_mb_amr on 2 ranks (patches of block 0 cut at the rank boxes, remote fills and flux matching across the interfaces); bitwise the run on 1 rank.",
      "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr", "mpi"], "covers": ["remote match entries with side2/sgn, gfl0 of multi-block grids"]})
case("sflow_2d_mb_amr_place",
     {"base": "sflow_2d_mb_amr", "np": 2, "ctrl_add": ["G 40 1"], "description": "As sflow_2d_mb_amr with placed patches (G 40 1): rectangles split per block and placed for the load of the ranks.",
      "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr", "mpi", "place"], "covers": ["place_patches with the clip boxes of the blocks"]})
case("sflow_2d_mb_amr_adaptive",
     {"base": "sflow_2d_mb_quadrants", "steps": 60, "ctrl_add": ["G 1 2", "G 2 4", "G 3 2", "G 4 8", "G 20 0.02"],
      "description": "As sflow_2d_mb_quadrants with two levels following the surface jumps (G 20), a regrid every 4 steps: the refined region moves across the interfaces. The same solution as one block to round-off (p2tests/t5).",
      "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr", "adaptive"],
      "covers": ["tagging of the blocks, regrid with the blocks kept first, level-0 refill before the regrid (regrid_prepare), ini_patch_old in the blocks"]})

# two blocks, the second turned by 90 deg, dam break running up a slope (wetting and drying), shoreline refinement
n2 = 64
I2, J2 = np.meshgrid(np.arange(n2+1), np.arange(n2+1), indexing="ij")
X2 = I2*0.5 + 0.1*J2*0.5; Y2 = J2*0.5 + 0.3*np.sin(np.pi*I2/n2)*np.sin(2*np.pi*J2/n2)
Z2 = 0.08*Y2 + 0.3*np.exp(-((X2-20.0)**2+(Y2-10.0)**2)/8.0)
h2 = n2//2
def rot_right(A):
    return np.array([[A[n2-b, a] for b in range(h2+1)] for a in range(n2+1)])
two = [(X2[:h2+1], Y2[:h2+1], Z2[:h2+1]), (rot_right(X2), rot_right(Y2), rot_right(Z2))]
case("sflow_2d_mb_wetdry_amr",
     {"description": "SFLOW multi-block wetting and drying with refinement: dam break running up a slope across the interface of two blocks (the second turned by 90 deg), shoreline (G 22) and surface-jump (G 20) refinement, regrid every 4 steps, 2 ranks. The same solution as one block to round-off (p2tests/t7).",
      "solver": "sflow", "np": 2, "steps": 60, "files": files, "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr", "wetdry", "mpi"],
      "covers": ["wet-dry fills across a turned interface, shoreline tags of the blocks"]},
     "B 1 0.5\nB 10 0.0 32.0 0.0 32.0 0.0 0.5\n\nM 10 1\nM 20 2\n\nR 20 1\n",
     ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 10.0 22.0 2.0 8.0 2.0\n\n"
      "G 50 1\nG 1 2\nG 2 4\nG 3 2\nG 4 8\nG 20 0.02\nG 22 1\n\nN 47 0.3\n\nM 10 1\n\n"
      "P 51 15.0 6.0\nP 51 17.0 11.0\nP 51 25.0 9.0\n\nW 22 -9.81\n"), two)

# three blocks around a singular vertex (a hexagon split into three quadrilaterals)
m, R = 32, 16.0
O = np.array([20.0, 20.0])
V = [O + R*np.array([np.cos(np.pi/3*k), np.sin(np.pi/3*k)]) for k in range(6)]
a = np.linspace(0, 1, m+1)
A, Cc = np.meshgrid(a, a, indexing="ij")
hexb = []
for b in range(3):
    P0, P1, P2, P3 = O, V[2*b], V[(2*b+1) % 6], V[(2*b+2) % 6]
    Xh = (1-A)*(1-Cc)*P0[0] + A*(1-Cc)*P1[0] + (1-A)*Cc*P3[0] + A*Cc*P2[0]
    Yh = (1-A)*(1-Cc)*P0[1] + A*(1-Cc)*P1[1] + (1-A)*Cc*P3[1] + A*Cc*P2[1]
    hexb.append((Xh, Yh, 0.2*np.exp(-((Xh-26.0)**2+(Yh-22.0)**2)/10.0)))
case("sflow_2d_mb_singular",
     {"description": "SFLOW multi-block grid with a singular vertex: a hexagon split into three blocks (three cells around the centre node), dam break over the vertex, two levels in a box and following the surface jumps: degenerate ring corner cells solid, no refinement around the vertex (p2tests/t6).",
      "solver": "sflow", "np": 1, "steps": 60, "files": files, "tags": ["2d", "sflow", "curvilinear", "multiblock", "amr"],
      "covers": ["singular vertices (blocks_refine_setup), degenerate ring cells (lexer::geometry_ini), corner cells through one face"]},
     "B 1 0.5\nB 10 0.0 16.0 0.0 16.0 0.0 0.5\n\nM 10 1\nM 20 2\n\nR 20 1\n",
     ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 17.0 23.0 17.0 23.0 1.6\n\n"
      "G 50 1\nG 1 2\nG 2 4\nG 3 2\nG 4 8\nG 10 10.0 30.0 10.0 30.0\nG 20 0.02\n\nN 47 0.3\n\nM 10 1\n\n"
      "P 51 20.5 20.5\nP 51 26.0 22.0\nP 51 14.0 18.0\n\nW 22 -9.81\n"), hexb)
print("cases written")
