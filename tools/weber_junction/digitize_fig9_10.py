# Digitize the measured points (dots) of Huang, Weber & Lai (2002), Figs. 9 and 10 (RUN5 velocity profiles), from the
# 300 dpi images embedded in the PDF: pdfimages -f 8 -l 9 -png paper.pdf im (im-001.png Fig. 9, im-002.png Fig. 10).
# Per panel: the base line (z = 0), the z ticks (0 .. 0.3), the zero lines of the three stations; per station the drawn
# profile (the tallest connected curve) is tracked row by row (its edges median-filtered) and removed; what is left are
# the dots: isolated blobs (centroid) and the parts of dots that touch the curve (centre 2 px inside the outer edge).
# usage: python3 digitize_fig9_10.py im-001.png im-002.png data/fig9_10_profiles.json
import numpy as np, json, sys
from PIL import Image
from scipy import ndimage

def runs1d(row):
    xs = np.nonzero(row)[0]
    if len(xs) == 0: return []
    out = []; s = xs[0]; p = xs[0]
    for x in xs[1:]:
        if x != p + 1: out.append((s, p)); s = x
        p = x
    out.append((s, p)); return out

def panels(a):
    H, Wd = a.shape
    out = []
    for (y0, y1) in ((150, 260), (440, 560), (735, 850)):
        for (x0, x1) in ((0, 495), (495, Wd)):
            sub = a[y0:y1, x0:x1]
            base = y0 + int(np.argmax(sub.sum(1)))
            vl = []
            for x in range(x0, x1):
                n = 0; y = base
                while y > 0 and (a[y, x] or a[y-1, x]): n += 1; y -= 1
                if n > 100: vl.append((x, n))
            grp = []
            for x, n in vl:
                if grp and x - grp[-1][-1][0] <= 2: grp[-1].append((x, n))
                else: grp.append([(x, n)])
            xs = [(float(np.mean([g[0] for g in gg])), max(g[1] for g in gg)) for gg in grp]
            ax = xs[0][0]
            c = int(round(ax))
            c = max(range(c-2, c+3), key=lambda x: a[base-140:base, x].sum())
            ticks = [y for y in range(base-160, base+4) if a[y, c+2:c+7].sum() >= 4]
            tg = []
            for y in ticks:
                if tg and y - tg[-1][-1] <= 2: tg[-1].append(y)
                else: tg.append([y])
            tc = [float(np.mean(g)) for g in tg]
            pz = (tc[-1] - tc[0])/0.3
            full = [x for x, n in xs[1:] if n >= 144]
            sp = 104.5
            zl = []
            for k in (1, 2, 3):
                cand = [x for x in full if abs(x - (ax + k*sp)) < 11]
                zl.append(min(cand, key=lambda x: abs(x - (ax + k*sp))) if cand else ax + k*sp)
            out.append(dict(base=base, axis=ax, z0=tc[-1], pz=pz, zero=zl))
    return out

def station(a, p, zx, pu, win):
    top = int(p["z0"] - 0.352*p["pz"]); bot = int(p["z0"] - 0.006*p["pz"])
    zs = sorted([p["axis"]] + list(p["zero"]) + [p["axis"] + 4*104.5])
    k = zs.index(zx)
    xa = int(round(max(zx - win[0], zs[k-1] + 4))); xb = int(round(min(zx + win[1], zs[k+1] - 3)))
    w = a[top:bot, xa:xb].copy()
    jz = int(round(zx)) - xa
    w[:, max(jz-2, 0):jz+3] = False
    lab, n = ndimage.label(w, structure=np.ones((3, 3)))
    if n == 0: return []
    sl = ndimage.find_objects(lab)
    hts = [s[0].stop - s[0].start for s in sl]
    ic = int(np.argmax(hts))
    line = lab == ic + 1
    resid = np.zeros_like(w)
    if hts[ic] > 25:
        L = np.full(w.shape[0], np.nan); R = np.full(w.shape[0], np.nan)
        for y in range(w.shape[0]):
            xs = np.nonzero(line[y])[0]
            if len(xs): L[y] = xs.min(); R[y] = xs.max()
        Lm = ndimage.median_filter(np.nan_to_num(L, nan=-1), 9, mode="nearest")
        Rm = ndimage.median_filter(np.nan_to_num(R, nan=-1), 9, mode="nearest")
        for y in range(w.shape[0]):
            xs = np.nonzero(line[y])[0]
            for x in xs:
                if np.isnan(L[y]) or Lm[y] < 0: continue
                if x < Lm[y] - 0.5 or x > Rm[y] + 0.5: resid[y, x] = True
        other = w & ~line
    else:
        other = w; Lm = Rm = None
    pts = []
    # isolated blobs
    lo, no = ndimage.label(other, structure=np.ones((3, 3)))
    for i, s in enumerate(ndimage.find_objects(lo)):
        h = s[0].stop - s[0].start; ww = s[1].stop - s[1].start
        ar = (lo[s] == i+1).sum()
        if h > 7 or ww > 7 or ar < 6: continue
        cy, cx = ndimage.center_of_mass(lo[s] == i+1)
        pts.append((cy + s[0].start, cx + s[1].start, "iso"))
    # dots touching the curve
    lr, nr = ndimage.label(ndimage.binary_dilation(resid, np.ones((1, 2))) & resid | resid, structure=np.ones((3, 3)))
    for i, s in enumerate(ndimage.find_objects(lr)):
        m = lr[s] == i+1
        ar = m.sum(); h = s[0].stop - s[0].start
        if ar < 3 or h < 2 or h > 7: continue
        cy = ndimage.center_of_mass(m)[0] + s[0].start
        ys = int(round(cy))
        mid = 0.5*(Lm[ys] + Rm[ys])
        cxr = ndimage.center_of_mass(m)[1] + s[1].start
        if cxr > mid: cx = s[1].stop - 1 - 2.0
        else: cx = s[1].start + 2.0
        pts.append((cy, cx, "touch"))
    # merge points closer than 3 px (both sides of a dot on the curve)
    pts.sort()
    merged = []
    for q in pts:
        if merged and abs(q[0] - merged[-1][0]) < 3 and abs(q[1] - merged[-1][1]) < 4:
            mq = merged[-1]; merged[-1] = (0.5*(q[0] + mq[0]), 0.5*(q[1] + mq[1]), "touch")
        else: merged.append(q)
    out = []
    for cy, cx, kind in merged:
        z = (p["z0"] - (cy + top))/p["pz"]
        u = (cx + xa - zx)/pu
        out.append((round(z, 4), round(u, 3), kind))
    return sorted(out)

def digitize(f, scale_val, win):
    a = np.array(Image.open(f)) == 0
    best = 0
    for y in range(25, 45):
        xs = np.nonzero(a[y, :200])[0]
        if len(xs) > best: best = len(xs); pu = (xs.max() - xs.min() + 1)/scale_val
    res = []
    for p in panels(a):
        res.append([station(a, p, zx, pu, win) for zx in p["zero"]])
    return res, pu

XS = (1.0, 0.0, -1.0, -2.0, -6.0)
if __name__ == "__main__":
    # usage: dig910.py fig9.png fig10.png out.json [raw.json]
    out = {"source": "Huang, Weber & Lai (2002), J. Hydraul. Eng. 128(3), Figs. 9 and 10, RUN5 (q* 0.75), experiment "
                     "(Shumate 1998): the dots, digitized from the 300 dpi images of the PDF; reading +-2 px "
                     "(Fig. 9 +-0.05 U0, Fig. 10 +-0.015 U0), points on top of the drawn model curve partly missed",
           "coords": "Huang: x/W positive upstream (the branch at -1 .. 0), y/W from the branch-side wall of the main "
                     "channel (negative in the branch), z/W from the bed; U0 = 0.628 m/s",
           "stations": []}
    raw = {}
    for fig, f, sv, qs, win in (("fig9", sys.argv[1], 2.0, ("-u", "v"), (15, 80)), ("fig10", sys.argv[2], 0.6, ("v", "-u"), (40, 70))):
        r, pu = digitize(f, sv, win)
        raw[fig] = dict(px_per_U0=pu, panels=r)
        for ip, pd in enumerate(r):
            for k, pts in enumerate(pd):
                if ip < 5: x, y, q = XS[ip], (0.25, 0.5, 0.75)[k], qs[0]
                else: x, y, q = (-0.75, -0.5, -0.25)[k], -1.0, qs[1]
                if not pts: continue
                out["stations"].append(dict(fig=fig, panel="abcdef"[ip], x=x, y=y, q=q,
                                            z=[p[0] for p in pts], val=[p[1] for p in pts]))
        print(fig, "px/U0", pu, [[len(s) for s in p] for p in r])
    out["note_fig10"] = ("the scale bar of Fig. 10 reads v/U0 = 6.0; read as 0.6 (with 6.0 the transverse velocity at "
                         "x/W = -1 would be 3.4 U0, ten times the branch inflow; with 0.6 it is 0.34 U0); x/W = 1 of "
                         "Fig. 10: the dots lie on the zero line (v = 0), not digitized")
    json.dump(out, open(sys.argv[3], "w"), indent=0)
    if len(sys.argv) > 4: json.dump(raw, open(sys.argv[4], "w"))
