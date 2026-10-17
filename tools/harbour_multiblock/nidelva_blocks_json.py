# P5d: the Nidelva - harbour grid of nidelva_harbour_grid.py as the block topology of LAGOON (blocks.json): the
# corners, the edges (the banks of the river corridor through their points, the straight sides of the harbour
# blocks, the river mouth shared by the corridor and harbour block 1) and the three blocks with their open faces;
# a case folder with control.txt, ctrl.txt and the schematic bed as geo.dat (10 m), for LAGOON's Blocks tab
# (Make grid writes blocks.dat and the input files)
# usage: python3 nidelva_blocks_json.py griddir casedir [--points 40] [--q 100] [--wl 10]
import argparse, json, os, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("grid", help="folder of nidelva_harbour_grid.py (grid.npz, grid.json)")
ap.add_argument("out")
ap.add_argument("--points", type=int, default=40, help="points shaping each bank")
ap.add_argument("--q", type=float, default=100.0)
ap.add_argument("--wl", type=float, default=10.0)
ap.add_argument("--geo-step", type=int, default=2, help="every n-th raster point into geo.dat")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
g = np.load(os.path.join(a.grid, "grid.npz"))
meta = json.load(open(os.path.join(a.grid, "grid.json")))
CX, CY = g["CX"], g["CY"]                      # (ni + 1, nj + 1): i along the river, j to its left
BX, BY = g["BX"], g["BY"]                      # (nI + 1, nJ + 1): i north, j west
ni, nj = CX.shape[0] - 1, CX.shape[1] - 1
nI, nJ = BX.shape[0] - 1, BX.shape[1] - 1
split = meta.get("harbour_blocks", 2)
cut = int(round(nJ / split))                   # the two harbour blocks: j 0 .. cut (east), cut .. nJ (west)
m0, m1 = meta["mouth_j_of_box"]                # the mouth: j m0 .. m1 of the box's south face

def P(x, y):
    return [round(float(x), 4), round(float(y), 4)]

def bank(j):
    k = np.linspace(0, ni, a.points + 2).round().astype(int)[1:-1]
    return [P(CX[i, j], CY[i, j]) for i in k]

vertices = {
    "up_right": P(CX[0, 0], CY[0, 0]), "up_left": P(CX[0, nj], CY[0, nj]),
    "mouth_east": P(BX[0, m0], BY[0, m0]), "mouth_west": P(BX[0, m1], BY[0, m1]),
    "se": P(BX[0, 0], BY[0, 0]), "s_mid": P(BX[0, cut], BY[0, cut]), "sw": P(BX[0, nJ], BY[0, nJ]),
    "ne": P(BX[-1, 0], BY[-1, 0]), "n_mid": P(BX[-1, cut], BY[-1, cut]), "nw": P(BX[-1, nJ], BY[-1, nJ]),
}
assert np.allclose(vertices["mouth_east"], P(CX[-1, 0], CY[-1, 0]), atol=1e-3), "the corridor ends on the box"
edges = {
    "bank_right": {"start": "up_right", "end": "mouth_east", "cells": ni, "points": bank(0)},
    "bank_left": {"start": "up_left", "end": "mouth_west", "cells": ni, "points": bank(nj)},
    "upstream": {"start": "up_right", "end": "up_left", "cells": nj, "curve": "straight"},
    "mouth": {"start": "mouth_east", "end": "mouth_west", "cells": nj, "curve": "straight"},
    "south_east": {"start": "se", "end": "mouth_east", "cells": m0, "curve": "straight"},
    "south_mid": {"start": "mouth_west", "end": "s_mid", "cells": cut - m1, "curve": "straight"},
    "south_west": {"start": "s_mid", "end": "sw", "cells": nJ - cut, "curve": "straight"},
    "north_east": {"start": "ne", "end": "n_mid", "cells": cut, "curve": "straight"},
    "north_west": {"start": "n_mid", "end": "nw", "cells": nJ - cut, "curve": "straight"},
    "east": {"start": "se", "end": "ne", "cells": nI, "curve": "straight"},
    "middle": {"start": "s_mid", "end": "n_mid", "cells": nI, "curve": "straight"},
    "west": {"start": "sw", "end": "nw", "cells": nI, "curve": "straight"},
}
blocks = [
    {"name": "Nidelva", "faces": {"i0": ["upstream"], "i1": ["mouth"], "j0": ["bank_right"], "j1": ["bank_left"]},
     "bc": {"i0": {"type": "inflow", "value": a.q}}, "smoothing": 300},
    {"name": "harbour east", "faces": {"i0": ["south_east", "mouth", "south_mid"], "i1": ["north_east"],
                                       "j0": ["east"], "j1": ["middle"]},
     "bc": {"i1": {"type": "outflow", "value": a.wl}}, "smoothing": 0},
    {"name": "harbour west", "faces": {"i0": ["south_west"], "i1": ["north_west"], "j0": ["middle"], "j1": ["west"]},
     "bc": {"i1": {"type": "outflow", "value": a.wl}}, "smoothing": 0},
]
topology = {"version": 1, "vertices": vertices, "edges": edges, "blocks": blocks, "bed": "geo", "z": 0.0,
            "note": "Nidelva - Trondheim harbour (P5d); the bed in geo.dat is SCHEMATIC, datum + %g m" % a.wl}
json.dump(topology, open(os.path.join(a.out, "blocks.json"), "w"), indent=1)

# the schematic bed as geo.dat (x y z), every n-th raster point
xs, ys, bed = g["xs"], g["ys"], g["bed"]
s = a.geo_step
X, Y = np.meshgrid(xs[::s], ys[::s], indexing="ij")
np.savetxt(os.path.join(a.out, "geo.dat"), np.column_stack([X.ravel(), Y.ravel(), bed[::s, ::s].ravel()]), fmt="%.2f")
open(os.path.join(a.out, "control.txt"), "w").write(
    "// Nidelva - Trondheim harbour: the multi-block grid of blocks.json (LAGOON: Set-up -> Blocks -> Make grid)\n"
    "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 1 3\nB 10 -400 2300 -2200 2000 0.0 15\n\nM 10 2\n\nR 20 1\n")
open(os.path.join(a.out, "ctrl.txt"), "w").write(
    "// Nidelva - Trondheim harbour, SFLOW (schematic bed)\nA 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"
    "B 50 0.05\nF 50 2\nF 60 %g\nG 50 1\nN 41 3600\nN 47 0.3\nM 10 2\nP 181 -1\nP 182 300\nW 22 -9.81\n" % a.wl)
print("blocks.json: %d corners, %d edges, %d blocks; geo.dat %d points" % (len(vertices), len(edges), len(blocks), X.size))
