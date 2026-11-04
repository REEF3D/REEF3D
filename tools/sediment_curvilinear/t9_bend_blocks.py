# P7 T9: the lags of the bend morphology across block interfaces: the LFM bend of T6/T8 (DS 0.1) as one block and as
# two blocks along the flow (block 1 the inflow reach with the discharge, G 54 1 1 Q; block 0 the bend and the outflow
# reach with the native outflow; the sediment does not leave through an open face G 54 of a block): bed, water level
# per cell after STEPS steps; LAYOUT=three: block 1 the inflow reach, block 0 the bend, block 2 the outflow reach with the
# water level (G 54 2 2 h; P8a: the open faces of the blocks pass the sediment as the native in- and outflow), for the
# bend terms
#   flow   S 95 1 with the bed held (S 35 0): the flow of the block split
#   s97_0, s97_01   S 95 2 with S 97 0 and 0.1
#   s951   S 95 1 (no lags: the reference of the block split)
#   s952   S 95 2, S 97 1 (the spiral flow builds up over L_a = 1.5 m across the interface at the bend entrance)
#   susp   S 95 2, S 64 0.5 (L_s of the suspended load up to 1.2 m)
# and the passes of the lags per step (the log of the block run).  The bed volume of both against its change.
# usage: python3 t9_bend_blocks.py [steps]; environment CASES (s951 s952 susp; flow s97_0 s97_01), NA (cells of the inflow block, 80), LAYOUT (two, three), NOGUARD (no fixed inflow reach S 71), NB (cells of the outflow block of three, 40), MF (100), NORUN
import sys, re, numpy as np
from p5lib import *
sys.path.insert(0, os.path.join(S5, "../nhflow_multiblock"))
from p4lib import run4

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
cases = os.environ.get("CASES", "s951 s952 susp").split()
R, B, h0, Q, d50, Cz, slope = 4.25, 1.7, 0.2, 0.17, 0.00078, 26.4, 0.0018
L1, L2, ds = 11.0, 6.0, 0.1
n1, n2, nb = int(round(L1/ds)), int(round(L2/ds)), int(round(np.pi*R/ds))
ny = int(round(B/ds)) | 1
nx = n1 + nb + n2
NA = int(os.environ.get("NA", "80"))
NB = int(os.environ.get("NB", "40"))
LAYOUT = os.environ.get("LAYOUT", "two")
g = 9.81
ks = 12.0*h0/np.exp(Cz/np.sqrt(g)/2.5)
ks_m = (20.0*h0**(1.0/6.0)/Cz)**6

s_n = np.concatenate([np.linspace(0, L1, n1+1), L1 + np.linspace(0, np.pi*R, nb+1)[1:], L1 + np.pi*R + np.linspace(0, L2, n2+1)[1:]])
n_n = np.linspace(-0.5*B, 0.5*B, ny+1)
S, N = np.meshgrid(s_n, n_n, indexing="ij")
X, Y = np.zeros_like(S), np.zeros_like(S)
m1, m2 = S <= L1, S >= L1 + np.pi*R
mb = ~(m1 | m2)
X[m1], Y[m1] = S[m1], N[m1]
phi = (S[mb] - L1)/R - 0.5*np.pi
r = R - N[mb]
X[mb], Y[mb] = L1 + r*np.cos(phi), R + r*np.sin(phi)
t = S[m2] - L1 - np.pi*R
X[m2], Y[m2] = L1 - t, 2*R - N[m2]
Z = slope*(s_n[-1] - S)
AR = 0.5*np.abs((X[1:, 1:]-X[:-1, :-1])*(Y[:-1, 1:]-Y[1:, :-1]) - (X[:-1, 1:]-X[1:, :-1])*(Y[1:, 1:]-Y[:-1, :-1]))

XC = {"eh": "S 95 1;S 11 4", "ehplain": "S 95 0;S 86 0;S 11 4", "ehpnoneq": "S 95 0;S 86 0;S 11 4;S 33 1;S 40 0.3", "ehpnoneqt": "S 95 0;S 86 0;S 11 4;S 33 1;S 40 0.3;S 49 5000", "ehsnoneq": "S 95 0;S 86 1;S 11 4;S 33 1;S 40 0.3", "ehsusp": "S 95 0;S 86 0;S 11 4;S 64 0.5;S 65 4", "ehslope": "S 95 0;S 86 1;S 11 4", "ehnoneq": "S 95 1;S 11 4;S 33 1;S 40 0.3", "ehnoneq2": "S 95 1;S 11 4;S 33 1;S 40 0.3;S 49 3000", "noneq": "S 95 1;S 33 1;S 40 0.3", "flow": "S 95 1;S 35 0", "plain": "S 95 0;S 86 0", "noslope": "S 95 1;S 86 0", "s97_0": "S 95 2;S 97 0", "s97_01": "S 95 2;S 97 0.1", "s951": "S 95 1", "s952": "S 95 2;S 97 1", "susp": "S 95 2;S 64 0.5"}
LAY = {"one": ([(0, nx)], 1, 2, ""),
       "two": ([(NA, nx), (0, NA)], 21, 2, "G 54 1 1 %g\n" % Q),
       "three": ([(NA, nx-NB), (0, NA), (nx-NB, nx)], 21, 21, "G 54 1 1 %g\nG 54 2 2 %g\n" % (Q, h0))}
LAY["onep"] = LAY["one"]
LAY = {k: v for k, v in LAY.items() if k in ("one", LAYOUT)}
# PERTURB=i,j,eps: the bed node (i, j) of the run LAYOUT=onep raised by eps (how fast the flow and the bed amplify a
# round-off difference there)
if os.environ.get("PERTURB") and LAYOUT == "onep":
    pi_, pj_, pe_ = os.environ["PERTURB"].split(",")
    ZP = Z.copy(); ZP[int(pi_), int(pj_)] += float(pe_)
else:
    ZP = Z

def control(c11, c14):
    return ("C 11 %d\nC 12 21\nC 13 21\nC 14 %d\nC 15 21\nC 16 3\n\nB 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 1\n\nR 20 1\n"
            % (c11, c14, ds, nx*ds, ny*ds))

def ctrl(xc):
    return ("A 10 2\nA 209 0\nA 210 3\nA 211 4\nA 218 1\nA 220 0\nA 240 1\nB 50 %g\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nG 50 1\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\n"
            % (ks_m, Q, h0, h0)
            + "S 10 1\nS 11 2\nS 16 2\nS 20 %g\nS 21 %g\nS 15 3\nS 17 1.0\nS 35 %s\nS 41 1\nS 43 1000\nS 44 1\nS 31 1\nS 100 0\n"
            % (d50, ks/d50, os.environ.get("MF", "100"))
            + ("" if os.environ.get("NOGUARD") else "S 71 2.0\n")
            + "S 86 1\nS 87 9\nS 88 0.5\nS 89 0.3\nS 96 0.4\n" + xc.replace(";", "\n") + "\n")

def cells(d, bl):
    """bed and water level per cell (i along s, j across) of a run with the AMR dump; one block: the regression dump"""
    if not os.path.exists(os.path.join(d, "amr", "amr_r0.bin")):
        st = read_state(sorted(glob.glob(os.path.join(d, "dump", "state_*_r0.bin")))[-1])
        return np.asarray(st["bed"]).reshape(nx, ny), np.asarray(st["WL"]).reshape(nx, ny)
    ring = [0, 0, 0, 0]
    for l in log_lines(d, "ring"):
        m = re.search(r"ring (\d+) (\d+) (\d+) (\d+)", l)
        if m:
            ring = [int(v) for v in m.groups()]
    bed, wl = np.full((nx, ny), np.nan), np.full((nx, ny), np.nan)
    for G in amr_grids2(d, 1):
        if G["lev"] != 0:
            continue
        for a in range(G["nx"]):
            for b in range(G["ny"]):
                if G["blk"] == 0:
                    io, jo = G["I0"]+a-ring[0], G["J0"]+b-ring[2]
                    if io < 0 or io >= bl[0][1]-bl[0][0] or jo < 0 or jo >= ny:
                        continue
                    k = (bl[0][0]+io, jo)
                else:
                    k = (bl[G["blk"]][0]+a, b)
                bed[k], wl[k] = G["bed"][a, b], G["WL"][a, b]
    return bed, wl

import glob
z0 = 0.25*(Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])
for c in cases:
    res = {}
    for name, (bl, c11, c14, g54) in LAY.items():
        d = os.path.join(S5, "runs/t9_%s_%s" % (c, name))
        rc = 0
        if not os.environ.get("NORUN"):
            rc = run4(d, control(c11, c14), ctrl(XC[c]) + g54 + ("G 1 1\n" if name in ("one", "onep") else ""),
                      [tuple(A[a:b+1].copy() for A in (X, Y, (ZP if name == "onep" else Z))) for a, b in bl], np_=1, steps=steps)
        res[name] = cells(d, bl)
        ps = [int(m.group(1)) for l in log_lines(d, "lags across the grids") for m in [re.search(r"passes (\d+)", l)] if m]
        print("%-5s %-5s rc=%d %s" % (c, name, rc, ("passes per step: mean %.1f max %d" % (np.mean(ps), max(ps))) if ps else ""),
              log_lines(d, "!!!")[:2])
    (b1, w1), (b3, w3) = res["one"], res[LAYOUT]
    db = b1 - z0
    print("%-5s bed change max %.3e m; %s vs one: bed %.1e  water level %.1e  (bed change rms %.1e)"
          % (c, np.nanmax(np.abs(db)), LAYOUT, np.nanmax(np.abs(b3-b1)), np.nanmax(np.abs(w3-w1)), np.sqrt(np.nanmean((b3-b1)**2))))
    print("      bed volume change: one %.3e m3, blocks %.3e m3" % (np.nansum(AR*(b1-z0)), np.nansum(AR*(b3-z0))))
