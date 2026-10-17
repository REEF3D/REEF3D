# writes the regression cases of the open faces of the blocks 1.. of a multi-block grid (P5d, G 54): the river reach
# of t1_open_faces.py as three blocks, the discharge inflow through the -i face of block 1, the water level at the +i
# face of block 2, no in- or outflow on the native grid (Regression/cases/*_mb_open_faces*)
import os, sys, json, numpy as np
T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(T, "../sediment_curvilinear"))
sys.path.insert(0, os.path.join(T, "../sflow_multiblock"))
from p5lib import nodes
from p2lib import write_blocks
C = os.path.join(T, "../../Regression/cases")

nx, ny, nz, dx, NA = 48, 10, 4, 0.1, 16
L, W, h = nx*dx, ny*dx, 0.25
Q = 0.4*W*h
X, Y = nodes(nx, ny, dx, "rot", 1000.0, 2000.0, 30.0, wavy=0.25)
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
Z = np.maximum(0.08 - np.hypot(I*dx-2.4, J*dx-0.5*W)*0.15, 0.0)
blocks = [tuple(A[a:b+1].copy() for A in (X, Y, Z)) for a, b in ((NA, nx-NA), (0, NA), (nx-NA, nx))]

def control(np_):
    return ("C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 %d\n\nR 20 1\n"
            % (dx, L, W, nz*0.1, np_))

def ctrl(flow, np_, extra=""):
    return (flow + "B 50 0.01\nB 60 1\nW 10 %g\nF 50 2\nF 60 %g\nF 62 %g\nG 50 1\nG 54 1 1 %g\nG 54 2 2 %g\n"
            "N 47 0.3\nN 48 0\nN 49 0.005\nM 10 %d\nP 10 0\nW 22 -9.81\n" % (Q, h, h, Q, h, np_) + extra)

def case(name, solver, np_, steps, ctrl_txt, desc, covers, tags):
    d = os.path.join(C, name)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "control.txt"), "w").write(control(np_))
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl_txt)
    write_blocks(os.path.join(d, "blocks.dat"), blocks)
    json.dump({"np": np_, "steps": steps, "tags": tags, "description": desc, "covers": covers, "solver": solver,
               "files": ["blocks.dat"]}, open(os.path.join(d, "case.json"), "w"), indent=1)

SF = "A 10 2\nA 209 0\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"
case("sflow_2d_mb_open_faces", "sflow", 1, 200, ctrl(SF, 1),
     "SFLOW on a multi-block grid (G 50 1): a river reach on a skewed grid turned by 30 deg (a mound on the bed, friction) as "
     "three blocks along the river, the discharge inflow through the -i face of block 1 and the water level at the +i face "
     "of block 2 (G 54), none on the native grid: the same water level as one block with the native in- and outflow (B 60 1, "
     "F 62) to round-off (harbour_multiblock/t1).",
     ["G 54: the open faces of the blocks 1.. (reefamr::blocks_open, the bc of the 2D lists in build_bc2D)",
      "ioflow_amr_block: SFLOW discharge, inflow along the face normal, water level", "sflow_HLL::open_io"],
     ["2d", "sflow", "curvilinear", "multiblock", "river"])
case("sflow_2d_mb_open_faces_amr", "sflow", 1, 100, ctrl(SF, 1, "G 1 1\nG 2 0\nG 10 1001.0 1004.0 1997.0 2005.0\n"),
     "As sflow_2d_mb_open_faces with one refinement level over the whole reach (G 10): the patches cross the block "
     "interfaces and keep away from the open faces (the no-refinement band of the in- and outflow).",
     ["reefamr::blocks_refine_setup: no refinement next to the open faces of the blocks"],
     ["2d", "sflow", "curvilinear", "multiblock", "river", "amr"])
NH = "A 10 5\nA 520 2\n"
case("nhflow_3d_mb_open_faces_mpi2", "nhflow", 2, 100, ctrl(NH, 2),
     "As sflow_2d_mb_open_faces with non-hydrostatic NHFLOW (A 520 2) on two ranks: the open faces of the blocks as the "
     "native in- and outflow (the log profile of the inflow scaled to the discharge, the river outflow of start4V, the "
     "pressure rows of the IO cells), the same solution as one block to round-off (harbour_multiblock/t1).",
     ["ioflow_amr_block: NHFLOW discharge, log inflow profile, inflow and outflow levels",
      "nhflow_amr::open_lists (gcin, gcout, IO of the open faces); lexer::river_xm, river_xp in start*V and nhflow_curv"],
     ["3d", "nhflow", "curvilinear", "multiblock", "river", "mpi"])
print("written")
