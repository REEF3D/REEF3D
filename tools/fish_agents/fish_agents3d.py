#!/usr/bin/env python3
# Agent fish in 3D: many fish, each driven by the brain of REEF3D (the rule brain of 6DOF_brain.cpp in its Python
# form, tools/fish_brain/brain_server.py, in 3D: it steers in the horizontal plane), with the cheap body of
# fish_agents.py, in a 3D flow computed by REEF3D (a sea cage, a channel).
#
#   python3 fish_agents3d.py case.json [--out tracks.csv] [--summary summary.json] [--after 60]
#
# The body: a point with a heading in the horizontal plane (fish_agents.body_step: the gait follows the brain,
# the speed s U_g follows f a L with tau_up / tau_down, coasts with tau_coast, brakes with C_brake; the heading
# turns with k b U/L). The depth: the fish swims towards the depth of its target with a vertical speed
# (z_target - z)/tau_z, at most w_max (a salmon changes depth slowly, by its swim bladder and by pitching).
# The fish moves with U along its heading plus the vertical speed plus the velocity of the flow at its position.
# The flow: uniform, or the velocity of a 3D REEF3D run (REEF3D_CFD_VTU), trilinear on its grid, one output or a
# series looped in time.
# The senses of the brain: the track, the target, 7 vision rays in the horizontal plane from the head to the
# walls of the box that are seen (an inflow or outflow is not), to vertical solid cylinders and to the nets of
# cages (vertical cylinders seen from inside and outside, within their depth); no lateral line. The other fish
# are not seen; two fish closer than 2 r_fish are pushed apart. A net holds the fish in its cage (the bottom of
# the net at its lower edge).
# Targets: a fixed point, or a point that circles the cage ("orbit"): the fish swim round the cage as farmed
# salmon do, each on its own radius and depth, at the speed U over the ground.
#
# case.json:
#   {"T": 120, "dt": 0.02, "out_dt": 0.5, "seed": 1,
#    "body": {"L": 1.0, "f0": 1.0, "s": 0.23, "k": 1.1, "tau_up": 4.1, "tau_down": 4.5, ...},   (fish_agents.BODY)
#    "brain": {"U_des": 0.2, "b_max": 1.0, "D": 1.5, "s0": 0.23, "tau_u": 0.8, "tau_d": 4.5},
#    "depth": {"tau": 5.0, "w_max": 0.05},
#    "flow": {"type": "uniform", "u": [0.1, 0.0, 0.0]}  or  {"type": "reef3d", "run": "<dir>", "t0": 300, "t1": 400},
#    "box": [x0, x1, y0, y1, z0, z1], "walls": [0, 0, 1, 1],
#    "cylinders": [[cx, cy, r, z_bottom, z_top]], "nets": [[cx, cy, r, z_bottom, z_top]],
#    "fish": [{"x": .., "y": .., "z": .., "h": deg, "target": [x, y, z]}, ...]
#         or {"n": 60, "orbit": {"c": [cx, cy], "r": [r0, r1], "z": [z0, z1], "U": 0.15, "dir": 1}}}
# Output: the tracks (t id x y z heading_deg U f a b mode) every out_dt; the summary per fish (distance to its
# target, mean frequency).
#
# Architect: Hans Bihs
import argparse, glob, json, math, os, random, re, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fish_agents as fa                                           # noqa: E402
bs = fa.bs
ANG = fa.ANG


# ---------------------------------------------------------------------------------------------------------------
# the flow
class UniformFlow3:
    def __init__(self, u):
        self.u = np.array(u, float)

    def __call__(self, t, x, y, z):
        return (np.full_like(x, self.u[0]), np.full_like(x, self.u[1]), np.full_like(x, self.u[2]))


class Reef3dFlow3:
    """the velocity of a 3D REEF3D run on its grid, outputs between t0 and t1, looped in time"""

    def __init__(self, run, t0=None, t1=None):
        byk = {}
        for f in sorted(glob.glob(os.path.join(run, "REEF3D_CFD_VTU", "REEF3D-CFD-*.vtu"))):
            byk.setdefault(int(re.search(r"REEF3D-CFD-(\d+)-\d+\.vtu$", f).group(1)), []).append(f)
        self.t, self.V = [], []
        for k in sorted(byk):
            P, V, tk = [], [], 0.0
            for f in byk[k]:
                tk, p, v = fa.read_vtu(f)
                P.append(p); V.append(v)
            if (t0 is not None and tk < t0 - 1e-9) or (t1 is not None and tk > t1 + 1e-9):
                continue
            P = np.vstack(P); V = np.vstack(V)
            ax = [np.unique(np.round(P[:, d], 6)) for d in range(3)]
            idx = [np.searchsorted(ax[d], np.round(P[:, d], 6)) for d in range(3)]
            G = np.zeros((len(ax[2]), len(ax[1]), len(ax[0]), 3), np.float32)
            G[idx[2], idx[1], idx[0]] = V
            self.xs, self.ys, self.zs = ax
            self.t.append(tk); self.V.append(G)
        if not self.t:
            raise SystemExit("no REEF3D output between t0 and t1 in " + run)
        self.t = np.array(self.t)
        self.period = (self.t[-1] - self.t[0]) + (self.t[1] - self.t[0] if len(self.t) > 1 else 1.0)

    def _interp(self, G, x, y, z):
        out = []
        for a, c in ((self.xs, x), (self.ys, y), (self.zs, z)):
            i = np.clip(np.searchsorted(a, c) - 1, 0, len(a) - 2)
            out.append((i, np.clip((c - a[i]) / (a[i + 1] - a[i]), 0, 1)))
        (i, tx), (j, ty), (k, tz) = out
        r = 0.0
        for dk, wk in ((0, 1 - tz), (1, tz)):
            for dj, wj in ((0, 1 - ty), (1, ty)):
                for di, wi in ((0, 1 - tx), (1, tx)):
                    r = r + G[k + dk, j + dj, i + di] * (wk * wj * wi)[:, None]
        return r[:, 0], r[:, 1], r[:, 2]

    def __call__(self, t, x, y, z):
        if len(self.t) == 1:
            return self._interp(self.V[0], x, y, z)
        tl = self.t[0] + (t % self.period)
        k = int(np.clip(np.searchsorted(self.t, tl) - 1, 0, len(self.t) - 1))
        k1 = (k + 1) % len(self.t)
        dtk = (self.t[k1] - self.t[k]) if k1 > k else (self.t[1] - self.t[0])
        w = float(np.clip((tl - self.t[k]) / dtk, 0, 1))
        a = self._interp(self.V[k], x, y, z); b = self._interp(self.V[k1], x, y, z)
        return tuple((1 - w) * a[c] + w * b[c] for c in range(3))

    def mean(self):
        return np.mean(self.V, axis=0)


# ---------------------------------------------------------------------------------------------------------------
# vision in the horizontal plane
def ray_lengths3(px, py, pz, dx, dy, box, walls, circles):
    """the distance along (dx, dy) from (px, py) to the box walls that are seen and to the first vertical
    circle (cx, cy, r, z_bottom, z_top) hit, from outside or inside, at the depth pz"""
    best = math.inf
    for lo, hi, p, d, seen in ((box[0], box[1], px, dx, walls[:2]), (box[2], box[3], py, dy, walls[2:])):
        if abs(d) > 1e-12:
            for wall, sn in zip((lo, hi), seen):
                if sn:
                    tt = (wall - p) / d
                    if tt > 0:
                        best = min(best, tt)
    for c in circles:
        if not (c[3] <= pz <= c[4]):
            continue
        ox, oy = px - c[0], py - c[1]
        bq = ox * dx + oy * dy
        disc = bq * bq - (ox * ox + oy * oy - c[2] * c[2])
        if disc >= 0.0:
            for tt in (-bq - math.sqrt(disc), -bq + math.sqrt(disc)):
                if tt > 0.0:
                    best = min(best, tt)
                    break
    return best


# ---------------------------------------------------------------------------------------------------------------
class Fish3:
    def __init__(self, i, x, y, z, h, target, body, brain, t_start):
        self.i, self.x, self.y, self.z, self.h = i, x, y, z, h
        self.target = target                   # a function of t: (x, y, z)
        self.g = {"f": [body["f0"], 0.0, body["f0"], 0.25], "a": [1.0, 0.0, 1.0, 0.25], "b": [0.0, 0.0, 0.0, 0.25]}
        self.Ug = 0.0
        self.second = None
        self.next = t_start
        self.st = {"side": 0.0, "refr": -1.0}
        self.pr = {"hello": 1, "body": i, "L": body["L"], "f0": body["f0"], "rho": 1000.0, "U_des": brain["U_des"],
                   "b_max": brain["b_max"], "D": brain["D"], "th_esc": 0.0, "t_start": t_start, "twoD": 0,
                   "s0": brain.get("s0", 0.0), "k_turn": brain.get("k_turn", 0.0),
                   "tau_u": brain.get("tau_u", 0.0), "tau_d": brain.get("tau_d", 0.0),
                   "n_u": brain.get("n_u", 0.0), "n_d": brain.get("n_d", 0.0)}
        self.acc_U = 0.0; self.acc_n = 0
        self.mode = -1


def orbit_target(c, r, z, w, ph):
    return lambda t: (c[0] + r * math.cos(ph + w * t), c[1] + r * math.sin(ph + w * t), z)


def run(case):
    rnd = random.Random(case.get("seed", 1))
    body = dict(fa.BODY, **case.get("body", {}))
    brain = dict({"U_des": 0.2, "b_max": 1.0, "D": 1.5}, **case.get("brain", {}))
    depth = dict({"tau": 5.0, "w_max": 0.05}, **case.get("depth", {}))
    L, f0 = body["L"], body["f0"]
    T, dt, out_dt = case.get("T", 60.0), case.get("dt", 0.02), case.get("out_dt", 0.5)
    box = case["box"]
    walls = case.get("walls", [1, 1, 1, 1])
    cyl = [list(c) for c in case.get("cylinders", [])]
    nets = [list(c) for c in case.get("nets", [])]
    fl = case.get("flow", {"type": "uniform", "u": [0.0, 0.0, 0.0]})
    flow = UniformFlow3(fl["u"]) if fl["type"] == "uniform" else Reef3dFlow3(fl["run"], fl.get("t0"), fl.get("t1"))
    t_start = case.get("t_start", max(1.0, 0.5 / f0))              # the brain starts at the end of the ramp

    fish = []
    spec = case["fish"]
    if isinstance(spec, dict):
        o = spec["orbit"]
        w_dir = o.get("dir", 1)
        for i in range(spec["n"]):
            r = rnd.uniform(*o["r"]); z = rnd.uniform(*o["z"]); ph = rnd.uniform(0.0, 2.0 * math.pi)
            w = w_dir * o["U"] / r
            tgt = orbit_target(o["c"], r, z, w, ph)
            x0, y0, _ = tgt(0.0)
            h = ph + w_dir * 0.5 * math.pi                     # along the orbit
            fish.append(Fish3(i, x0, y0, z, h, tgt, body, brain, t_start + rnd.uniform(0, 0.5 / f0)))
    else:
        for i, s in enumerate(spec):
            T_ = tuple(s["target"])
            fish.append(Fish3(i, s["x"], s["y"], s["z"], math.radians(s.get("h", 0.0)), (lambda t, T_=T_: T_), body, brain, t_start))

    rows = []
    t, t_out = 0.0, 0.0
    n = len(fish)
    while t < T - 1e-9:
        X = np.array([f.x for f in fish]); Y = np.array([f.y for f in fish]); Z = np.array([f.z for f in fish])
        uf, vf, wf = flow(t, X, Y, Z)
        for k, f in enumerate(fish):
            for key in "fab":
                v, r, tg, ts = f.g[key]
                lam = 4.744 / ts
                r += dt * (lam * lam * (tg - v) - 2 * lam * r); v += dt * r
                f.g[key] = [v, r, tg, ts]
            if f.second and t >= f.second[0]:
                for key, v in zip("fab", f.second[1:4]):
                    f.g[key][2] = v; f.g[key][3] = f.second[4]
                f.second = None
            fg, ag, bg = f.g["f"][0], f.g["a"][0], f.g["b"][0]
            f.Ug, kt = fa.body_step(f.Ug, fg, ag, bg, t, dt, body)
            Uw = body["s"] * f.Ug
            f.h += dt * kt * bg * Uw / L
            ex, ey = math.cos(f.h), math.sin(f.h)
            wz = max(-depth["w_max"], min(depth["w_max"], (f.target(t)[2] - f.z) / depth["tau"]))
            nx, ny, nz = f.x + dt * (Uw * ex + uf[k]), f.y + dt * (Uw * ey + vf[k]), f.z + dt * (wz + wf[k])
            # held in the box, outside the solids, inside the net of its cage
            nx = min(max(nx, box[0] + 0.05), box[1] - 0.05); ny = min(max(ny, box[2] + 0.05), box[3] - 0.05)
            nz = min(max(nz, box[4] + 0.05), box[5] - 0.05)
            for c in cyl:
                dd = math.hypot(nx - c[0], ny - c[1])
                if c[3] <= nz <= c[4] and dd < c[2] + 0.02:
                    s_ = (c[2] + 0.02) / max(dd, 1e-9); nx, ny = c[0] + (nx - c[0]) * s_, c[1] + (ny - c[1]) * s_
            for c in nets:
                if math.hypot(f.x - c[0], f.y - c[1]) < c[2]:          # in this cage
                    dd = math.hypot(nx - c[0], ny - c[1])
                    if dd > c[2] - 0.1:
                        s_ = (c[2] - 0.1) / dd; nx, ny = c[0] + (nx - c[0]) * s_, c[1] + (ny - c[1]) * s_
                    nz = min(max(nz, c[3] + 0.1), c[4] - 0.1)
            f.acc_U += ((nx - f.x) * ex + (ny - f.y) * ey) / dt; f.acc_n += 1
            f.x, f.y, f.z = nx, ny, nz
        if n > 1 and body["r_fish"] > 0.0:
            P = np.array([[f.x, f.y, f.z] for f in fish])
            D = P[:, None, :] - P[None, :, :]
            dd = np.linalg.norm(D, axis=2); np.fill_diagonal(dd, np.inf)
            ov = np.clip(2 * body["r_fish"] - dd, 0.0, None)
            if ov.any():
                w = 0.5 * ov / np.maximum(dd, 1e-9)
                P += (w[:, :, None] * D).sum(axis=1)
                for k, f in enumerate(fish):
                    f.x, f.y, f.z = P[k]
        t += dt
        for k, f in enumerate(fish):
            if t < f.next:
                continue
            ex, ey = math.cos(f.h), math.sin(f.h)
            hx, hy = f.x + 0.5 * L * ex, f.y + 0.5 * L * ey
            rays = [ray_lengths3(hx, hy, f.z, math.cos(f.h + an), math.sin(f.h + an), box, walls, cyl + nets) for an in ANG]
            T_ = f.target(t)
            sv = {"t": t, "f": f.g["f"][0], "a": f.g["a"][0], "b": f.g["b"][0], "pos": [f.x, f.y, f.z],
                  "axis": [ex, ey, 0.0], "lat": [-ey, ex, 0.0], "U": f.acc_U / max(1, f.acc_n),
                  "stride": f.st.get("kx", [0.5])[0], "target": list(T_), "rays": rays, "ray_angles": ANG}
            d = bs.rules(sv, f.st, f.pr)
            ts = d.get("ts", 0.0) or 0.5 / max(f.g["f"][0], 0.1)
            for key in "fab":
                f.g[key][2] = d[key]; f.g[key][3] = ts
            f.second = (d["t2"], d["f2"], d["a2"], d["b2"], d["ts2"]) if "t2" in d else None
            f.next = t + d.get("next", min(0.5 / d["f"], 0.5))
            f.mode = d.get("mode", 0)
            f.acc_U = 0.0; f.acc_n = 0
        if t >= t_out - 1e-9:
            for f in fish:
                rows.append((t, f.i, f.x, f.y, f.z, math.degrees(f.h), body["s"] * f.Ug, f.g["f"][0], f.g["a"][0], f.g["b"][0], f.mode))
            t_out += out_dt
    return fish, rows


def summarise(fish, rows, t0):
    R = np.array(rows)
    out = []
    for f in fish:
        m = (R[:, 1] == f.i) & (R[:, 0] > t0)
        tg = np.array([f.target(t) for t in R[m, 0]])
        d = np.sqrt(((R[m][:, 2:5] - tg) ** 2).sum(axis=1)) if len(tg) else np.array([])
        out.append({"id": f.i, "mean": float(d.mean()) if len(d) else None, "within_1": float((d < 1.0).mean()) if len(d) else None,
                    "f_mean": float(R[m, 7].mean()) if m.any() else None})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--out", default="tracks3d.csv")
    ap.add_argument("--summary", default="")
    ap.add_argument("--after", type=float, default=60.0, help="the summary over t > after")
    a = ap.parse_args()
    case = json.load(open(a.case))
    fish, rows = run(case)
    with open(a.out, "w") as fo:
        fo.write("t id x y z heading U f a b mode\n")
        for r in rows:
            fo.write("%.3f %d %.4f %.4f %.4f %.3f %.5f %.4f %.4f %.4f %d\n" % r)
    s = summarise(fish, rows, a.after)
    if a.summary:
        json.dump(s, open(a.summary, "w"), indent=1)
    print("%d fish, %.1f s: mean distance to the target after %.0f s %.2f m, within 1 m %.0f %%, mean frequency %.2f Hz" %
          (len(fish), case.get("T", 60.0), a.after, np.mean([x["mean"] for x in s]), 100 * np.mean([x["within_1"] for x in s]),
           np.mean([x["f_mean"] for x in s])))


if __name__ == "__main__":
    main()
