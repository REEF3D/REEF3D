# writes the multi-block NHFLOW regression cases (P4, G 50 1 with DIVEMesh R 20 1): Regression/cases/nhflow_3d_mb_*
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sflow_multiblock"))
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

# four blocks around a regular vertex, indices turned by 0, 90, 180, 270 deg, wavy skewed grid (t2)
n, nz, dx = 40, 6, 0.05
L = n*dx
h = n//2
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.15*y + 0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.05*x + 0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.06*np.exp(-((X-1.3)**2+(Y-0.7)**2)/0.04)
def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))
quads = [quad(0, 0, 0), quad(h, 0, 1), quad(0, h, 2), quad(h, h, 3)]
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\n\nR 20 1\n" % (dx, L, L, nz*dx)
ctrl = "A 10 5\nA 520 2\nF 60 0.3\nF 72 0.65 1.45 0.65 1.45 0.34\nG 50 1\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\n"
tags = ["3d", "nhflow", "curvilinear", "multiblock"]
case("nhflow_3d_mb_quadrants",
     {"description": "NHFLOW multi-block grid (G 50 1, DIVEMesh R 20 1): four curvilinear blocks around a regular vertex, their indices turned by 0, 90, 180 and 270 deg, wavy skewed grid, bumpy bed, non-hydrostatic dam break over the vertex (A 520 2). The same solution as one block (G 52 1) to round-off (nhflow_multiblock/t2).",
      "solver": "nhflow", "np": 1, "steps": 20, "files": ["blocks.dat"], "tags": tags + ["quick"],
      "covers": ["nhflow_amr_block: blocks 1.. as level-0 patches, fill_level0 (stage input of the blocks and of the ring), fill_bed0, blocks_ini",
                 "composite pressure over the blocks: the ring no unknown, pr_fill0, block-local V-cycles, the coarse space across the blocks (G 53)",
                 "G 52 1: the i-walls as the j-walls (pair_bc, iwall_sym, the Neumann rows of the pressure correction); block_node_any through the interface first"]},
     control, ctrl, quads)
case("nhflow_3d_mb_quadrants_mpi2",
     {"base": "nhflow_3d_mb_quadrants", "np": 2,
      "description": "As nhflow_3d_mb_quadrants on 2 ranks: blocks on other ranks than the native cells next to them (remote fills of the blocks and of the ring, the ring in the partition halo).",
      "tags": tags + ["mpi"], "covers": ["ring_halo, remote fills of fill_run0 in the composite solve, the porosity below the bed in the partition halo (nhflow_f::ini)"]})
case("nhflow_3d_mb_amr",
     {"base": "nhflow_3d_mb_quadrants", "ctrl_add": ["G 1 1", "G 2 0", "G 3 2", "G 4 8", "G 10 0.55 1.75 0.6 1.6"],
      "description": "As nhflow_3d_mb_quadrants with one refinement level in a static box over the vertex: patches in all four blocks, cut at the block faces, filled and flux-matched across the turned interfaces (sign of the fluxes), the composite pressure with the block-local cycles below the patches (nhflow_multiblock/t5).",
      "tags": tags + ["amr"], "covers": ["flux_hook with the sign of turned coarse faces, prolongation from the blocks, restriction into the blocks, FAC with the blocks as level 0"]})
case("nhflow_3d_mb_idiff",
     {"base": "nhflow_3d_mb_quadrants", "ctrl_add": ["A 512 2", "W 2 0.001"],
      "description": "As nhflow_3d_mb_quadrants with the implicit horizontal diffusion (A 512 2): one composite diffusion solve over the blocks (the ring no unknown, the cells around the blocks filled from the increments).",
      "tags": tags, "covers": ["nhflow_amr_diff: leaf lists without the ring, df_sync level-0 fill"]})

# potential flow initialisation (I 11) through a bypass over 4 blocks (t4 with NB 3): channel blocked by a dry obstacle
nx, ny0, ny1 = 120, 12, 24
H = (ny0+ny1)*dx
L = nx*dx
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny0+ny1+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.03*np.sin(np.pi*y/H)*np.sin(2*np.pi*x/L)
Y = y + 0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/H)
Z = np.where(((x > 2.4) & (x < 3.6) & (y < ny0*dx+1e-9)) | ((y > ny0*dx+1e-9) & ((x < 1.2) | (x > 4.8))), 0.6, 0.0)
offs = [(0, 0, 0, nx, ny0)] + [(b*40, ny0, (b % 3)+1, 40, ny1) for b in range(3)]
blocks = [tuple(np.rot90(A[i0:i0+mi+1, j0:j0+mj+1], k=-k).copy() for A in (X, Y, Z)) for (i0, j0, k, mi, mj) in offs]
control = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n\nB 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\n\nR 20 1\n" % (dx, L, H, nz*dx)
ctrl = "A 10 5\nA 520 2\nB 60 1\nF 60 0.3\nF 62 0.3\nG 50 1\nI 11 1\nN 47 0.3\nM 10 1\nP 10 0\nW 10 0.05\nW 22 -9.81\n"
case("nhflow_3d_mb_potential",
     {"description": "NHFLOW multi-block grid with in- and outflow on block 0 (DIVEMesh C 11 1, C 14 2 on the faces of block 0 without an interface): a channel blocked by a dry obstacle, the water goes around it through three turned blocks; potential flow initialisation (I 11) as one composite solve over the blocks, then a few steps with the discharge inflow and the stage outflow (nhflow_multiblock/t4).",
      "solver": "nhflow", "np": 1, "steps": 5, "files": ["blocks.dat"], "tags": tags + ["potential", "river"],
      "covers": ["nhflow_amr::potential, nhflow_potential_f::amr_assemble/amr_velocities, the curvilinear potential rows (lapcoef, inflow flux along the face normal, chain rule velocities)",
                 "DIVEMesh R 20 1: several interfaces on one face of block 0, boundary types of the faces of block 0 without an interface"]},
     control, ctrl, blocks)
print("written")
