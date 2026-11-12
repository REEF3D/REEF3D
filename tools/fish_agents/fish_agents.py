#!/usr/bin/env python3
# Agent fish: many fish, each driven by the brain of REEF3D (the rule brain of 6DOF_brain.cpp in its Python form,
# tools/fish_brain/brain_server.py), with a cheap body instead of a resolved one, in a flow computed by REEF3D.
#
#   python3 fish_agents.py case.json [--out tracks.csv] [--summary summary.json]
#
# The body (2D, in the x-z plane of a 2D REEF3D run): a point with a heading. Its gait (frequency f, amplitude
# factor a, curvature bias b) follows the targets of the brain with the critically damped dynamics of the gait
# control (X 199). The speed through the water U = s U_g along the heading, U_g follows f a L: rising with
# tau_up (thrust), falling with tau_down when the gait slows, with tau_coast when the tail is still (a < 0.5),
# with tau_brake when it is still and the body bent (|b| > 0.5), then not below U_brake (the bent body still
# drifts and turns). The heading turns with k b U/L (k_brake while it brakes). The fish moves
# with U along its heading plus the velocity of the flow at its position. These constants come from the
# resolved runs: the 2D swimmer of REEF3D (L = 1 m, 1 Hz): k 1.3 (validation 08), tau_coast 11 s (validation 13);
# s 0.49, tau_up = tau_down 0.58 s (fitted apart: 0.61 and 0.49 s), tau_brake 0.86 s, U_brake 0.13 m/s, k_brake
# 1.77 fitted to its station run with the brain of patch 0020 (validation 17, calibrate.py). The brain's own
# model of the speed (X 201 tau_u tau_d) is set apart from the body in "brain".
#
# The flow: uniform, or the velocity of a 2D REEF3D run (REEF3D_CFD_VTU), one output or a series looped in time,
# bilinear on the grid of the run.
# The senses: the track (the brain estimates the current from it), the target, 7 vision rays to the walls of
# the box (the sides flagged in walls; an inflow or outflow is not seen), the solid cylinders and the other
# fish (circles of radius r_fish, seen within D_fish of the head, 0: not seen); no lateral line. Two fish
# closer than 2 r_fish are pushed apart (the bodies do not overlap).
# A fish that would enter a solid or leave the box is held at its edge.
#
# case.json:
#   {"T": 30, "dt": 0.01, "out_dt": 0.1, "seed": 1,
#    "body": {"L": 1.0, "f0": 1.0, "s": 0.49, "k": 1.3, "k_brake": 1.77, "tau_up": 0.58, "tau_down": 0.58,
#             "tau_coast": 11.0, "tau_brake": 0.86, "U_brake": 0.13, "r_fish": 0.15, "D_fish": 0.0},
#    "brain": {"U_des": 0.5, "b_max": 1.0, "D": 1.5, "s0": 0.5, "k_turn": 1.3, "tau_u": 0.8, "tau_d": 0.8},
#    "flow": {"type": "uniform", "u": [0.2, 0.0]}  or  {"type": "reef3d", "run": "<dir>", "t0": 20, "t1": 40},
#    "box": [0, 8, 0, 6], "walls": [0, 0, 1, 1], "cylinders": [[4.0, 3.0, 0.4]],
#    "fish": [{"x": 1.6, "z": 3.0, "h": 0.0, "target": [3.5, 3.0]}, ...]
#         or {"n": 50, "region": [x0, x1, z0, z1], "h": [h0, h1], "target": [x, z] | "start", "spread": 0.5}}
# Output: the tracks (t id x z heading_deg U f a b mode) every out_dt; the summary per fish.
#
# Architect: Hans Bihs
import argparse, json, math, os, random, re, sys, glob
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fish_brain"))
import brain_server as bs                                          # noqa: E402

ANG = [0.0, 0.4363, -0.4363, 0.8727, -0.8727, 1.309, -1.309]       # the rays of the brain (sixdof_brain::rays)


# ---------------------------------------------------------------------------------------------------------------
# the flow
class UniformFlow:
    def __init__(self, u):
        self.u = np.array(u, float)

    def __call__(self, t, x, z):
        return np.full_like(x, self.u[0]), np.full_like(x, self.u[1])


def _appended(raw):
    head, data = raw.split(b'<AppendedData encoding="raw">', 1)
    return head.decode("latin-1"), data[data.index(b"_") + 1:]


def _block(data, off, dtype):
    n = int(np.frombuffer(data[off:off + 4], dtype="<u4")[0])
    return np.frombuffer(data[off + 4:off + 4 + n], dtype=dtype)


def read_vtu(path):
    """time, points, velocity of one piece of a REEF3D CFD output (binary appended VTU)"""
    head, data = _appended(open(path, "rb").read())
    t = float(re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', head).group(1))
    off_v = int(re.search(r'Name="velocity" NumberOfComponents="3" format="appended" offset="(\d+)"', head).group(1))
    off_p = int(re.search(r'<Points>\s*<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="(\d+)"', head).group(1))
    return t, _block(data, off_p, "<f4").reshape(-1, 3), _block(data, off_v, "<f4").reshape(-1, 3)


class Reef3dFlow:
    """the x-z velocity of a 2D REEF3D run on its grid, outputs between t0 and t1, looped in time"""

    def __init__(self, run, t0=None, t1=None):
        files = sorted(glob.glob(os.path.join(run, "REEF3D_CFD_VTU", "REEF3D-CFD-*.vtu")))
        byk = {}
        for f in files:
            k = int(re.search(r"REEF3D-CFD-(\d+)-\d+\.vtu$", f).group(1))
            byk.setdefault(k, []).append(f)
        self.t, self.U, self.W = [], [], []
        for k in sorted(byk):
            P, V, tk = [], [], 0.0
            for f in byk[k]:
                tk, p, v = read_vtu(f)
                P.append(p); V.append(v)
            if (t0 is not None and tk < t0 - 1e-9) or (t1 is not None and tk > t1 + 1e-9):
                continue
            P = np.vstack(P); V = np.vstack(V)
            m = np.abs(P[:, 1] - P[:, 1].min()) < 1e-9
            P, V = P[m], V[m]
            xs = np.unique(np.round(P[:, 0], 7)); zs = np.unique(np.round(P[:, 2], 7))
            ix = np.searchsorted(xs, np.round(P[:, 0], 7)); iz = np.searchsorted(zs, np.round(P[:, 2], 7))
            U = np.zeros((len(zs), len(xs))); W = np.zeros_like(U)
            U[iz, ix] = V[:, 0]; W[iz, ix] = V[:, 2]
            self.xs, self.zs = xs, zs
            self.t.append(tk); self.U.append(U); self.W.append(W)
        if not self.t:
            raise SystemExit("no REEF3D output between t0 and t1 in " + run)
        self.t = np.array(self.t)
        self.period = (self.t[-1] - self.t[0]) + (self.t[1] - self.t[0] if len(self.t) > 1 else 1.0)

    def _interp(self, F, x, z):
        xs, zs = self.xs, self.zs
        i = np.clip(np.searchsorted(xs, x) - 1, 0, len(xs) - 2)
        j = np.clip(np.searchsorted(zs, z) - 1, 0, len(zs) - 2)
        tx = np.clip((x - xs[i]) / (xs[i + 1] - xs[i]), 0, 1); tz = np.clip((z - zs[j]) / (zs[j + 1] - zs[j]), 0, 1)
        return (F[j, i] * (1 - tx) * (1 - tz) + F[j, i + 1] * tx * (1 - tz) + F[j + 1, i] * (1 - tx) * tz + F[j + 1, i + 1] * tx * tz)

    def __call__(self, t, x, z):
        if len(self.t) == 1:
            return self._interp(self.U[0], x, z), self._interp(self.W[0], x, z)
        tl = self.t[0] + (t % self.period)
        k = int(np.clip(np.searchsorted(self.t, tl) - 1, 0, len(self.t) - 1))
        k1 = (k + 1) % len(self.t)
        dtk = (self.t[k1] - self.t[k]) if k1 > k else (self.t[1] - self.t[0])
        w = np.clip((tl - self.t[k]) / dtk, 0, 1)
        u = (1 - w) * self._interp(self.U[k], x, z) + w * self._interp(self.U[k1], x, z)
        v = (1 - w) * self._interp(self.W[k], x, z) + w * self._interp(self.W[k1], x, z)
        return u, v


# ---------------------------------------------------------------------------------------------------------------
# geometry: rays to the box, to the cylinders and to the other fish
def ray_lengths(px, pz, dx, dz, box, circles, walls=(1, 1, 1, 1)):
    """the distance along (dx, dz) from (px, pz) to the box walls that are seen (walls: x0, x1, z0, z1; an
    inflow or outflow boundary is not seen) and to the first circle (cx, cz, r) hit"""
    best = math.inf
    for lo, hi, p, d, seen in ((box[0], box[1], px, dx, walls[:2]), (box[2], box[3], pz, dz, walls[2:])):
        if abs(d) > 1e-12:
            for wall, sn in zip((lo, hi), seen):
                if not sn:
                    continue
                tt = (wall - p) / d
                if tt > 0:
                    x, z = px + tt * dx, pz + tt * dz
                    if box[0] - 1e-9 <= x <= box[1] + 1e-9 and box[2] - 1e-9 <= z <= box[3] + 1e-9:
                        best = min(best, tt)
    if len(circles):
        ox, oz = px - circles[:, 0], pz - circles[:, 1]
        bq = ox * dx + oz * dz
        cq = ox * ox + oz * oz - circles[:, 2] ** 2
        disc = bq * bq - cq
        ok = disc >= 0
        if ok.any():
            tt = -bq[ok] - np.sqrt(disc[ok])
            tt = tt[tt > 0]
            if len(tt):
                best = min(best, float(tt.min()))
    return best


# ---------------------------------------------------------------------------------------------------------------
class Fish:
    def __init__(self, i, x, z, h, target, body, brain, t_start):
        self.i, self.x, self.z, self.h = i, x, z, h
        self.target = target
        self.g = {"f": [body["f0"], 0.0, body["f0"], 0.25], "a": [1.0, 0.0, 1.0, 0.25], "b": [0.0, 0.0, 0.0, 0.25]}
        self.Ug = 0.0                         # the speed of the gait (starts from rest)
        self.second = None
        self.next = t_start
        self.st = {"side": 0.0, "refr": -1.0}   # the state of the brain (as brain_server.py starts it)
        self.pr = {"hello": 1, "body": i, "L": body["L"], "f0": body["f0"], "rho": 1000.0, "U_des": brain["U_des"],
                   "b_max": brain["b_max"], "D": brain["D"], "th_esc": 0.0, "t_start": t_start, "twoD": 1,
                   "s0": brain.get("s0", 0.0), "k_turn": brain.get("k_turn", 0.0),
                   "tau_u": brain.get("tau_u", 0.0), "tau_d": brain.get("tau_d", 0.0)}
        self.acc_U = 0.0; self.acc_n = 0
        self.mode = -1


def run(case):
    rnd = random.Random(case.get("seed", 1))
    body = dict({"L": 1.0, "f0": 1.0, "s": 0.49, "k": 1.3, "k_brake": 1.77, "tau_up": 0.58, "tau_down": 0.58,
                 "tau_coast": 11.0, "tau_brake": 0.86, "U_brake": 0.13, "r_fish": 0.15, "D_fish": 0.0},
                **case.get("body", {}))
    brain = dict({"U_des": 0.5, "b_max": 1.0, "D": 1.5}, **case.get("brain", {}))
    L, f0 = body["L"], body["f0"]
    T, dt, out_dt = case.get("T", 20.0), case.get("dt", 0.01), case.get("out_dt", 0.1)
    box = case.get("box", [0, 8, 0, 6])
    walls = case.get("walls", [1, 1, 1, 1])                           # seen: x0, x1, z0, z1 (0 for inflow, outflow)
    cyl = np.array(case.get("cylinders", []), float).reshape(-1, 3)
    fl = case.get("flow", {"type": "uniform", "u": [0.0, 0.0]})
    flow = UniformFlow(fl["u"]) if fl["type"] == "uniform" else Reef3dFlow(fl["run"], fl.get("t0"), fl.get("t1"))
    t_start = case.get("t_start", 1.0 + 1.0 / f0)                    # the brain starts after the ramp and a beat

    fish = []
    spec = case["fish"]
    if isinstance(spec, dict):
        x0, x1, z0, z1 = spec["region"]
        h0, h1 = spec.get("h", [0.0, 0.0])
        for i in range(spec["n"]):
            for _ in range(1000):
                x, z = rnd.uniform(x0, x1), rnd.uniform(z0, z1)
                if all(math.hypot(x - c[0], z - c[1]) > c[2] + body["r_fish"] for c in cyl) and \
                        all(math.hypot(x - f.x, z - f.z) > 2 * body["r_fish"] for f in fish):
                    break
            tg = spec.get("target", "start")
            if tg == "start":
                tgt = (x, z)
            else:
                sp = spec.get("spread", 0.0)
                tgt = (tg[0] + rnd.uniform(-sp, sp), tg[1] + rnd.uniform(-sp, sp))
            fish.append(Fish(i, x, z, math.radians(rnd.uniform(h0, h1)), tgt, body, brain, t_start + rnd.uniform(0, 0.5 / f0)))
    else:
        for i, s in enumerate(spec):
            fish.append(Fish(i, s["x"], s["z"], math.radians(s.get("h", 0.0)), tuple(s["target"]), body, brain, t_start))

    rows = []
    t, t_out = 0.0, 0.0
    n = len(fish)
    while t < T - 1e-9:
        X = np.array([f.x for f in fish]); Z = np.array([f.z for f in fish])
        uf, wf = flow(t, X, Z)
        for k, f in enumerate(fish):
            # the gait: critically damped towards the targets
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
            ramp = min(1.0, t / 1.0)
            target_U = fg * ag * L * ramp
            if ag < 0.5:
                if abs(bg) > 0.5:
                    # braking: the bent body slows down fast, but not to rest (2D CFD: 0.1-0.2 m/s)
                    tau = body["tau_brake"]
                    target_U = max(target_U, body["U_brake"] / body["s"])
                else:
                    tau = body["tau_coast"]
            else:
                tau = body["tau_up"] if target_U >= f.Ug else body["tau_down"]
            f.Ug += dt / tau * (target_U - f.Ug)
            Uw = body["s"] * f.Ug
            f.h += dt * (body["k_brake"] if ag < 0.5 else body["k"]) * bg * Uw / L
            ex, ez = math.cos(f.h), math.sin(f.h)
            nx, nz = f.x + dt * (Uw * ex + uf[k]), f.z + dt * (Uw * ez + wf[k])
            # held at the box and outside the solids
            nx = min(max(nx, box[0] + 0.05), box[1] - 0.05); nz = min(max(nz, box[2] + 0.05), box[3] - 0.05)
            for c in cyl:
                dd = math.hypot(nx - c[0], nz - c[1])
                if dd < c[2] + 0.02:
                    s_ = (c[2] + 0.02) / max(dd, 1e-9)
                    nx, nz = c[0] + (nx - c[0]) * s_, c[1] + (nz - c[1]) * s_
            f.acc_U += ((nx - f.x) * ex + (nz - f.z) * ez) / dt; f.acc_n += 1
            f.x, f.z = nx, nz
        # the bodies do not overlap: pairs closer than 2 r_fish are pushed apart (no brain involved)
        if n > 1 and body["r_fish"] > 0.0:
            P = np.array([[f.x, f.z] for f in fish])
            dx = P[:, None, 0] - P[None, :, 0]; dz = P[:, None, 1] - P[None, :, 1]
            dd = np.hypot(dx, dz); np.fill_diagonal(dd, np.inf)
            ov = np.clip(2 * body["r_fish"] - dd, 0.0, None)
            if ov.any():
                w = 0.5 * ov / np.maximum(dd, 1e-9)
                P[:, 0] += (w * dx).sum(axis=1); P[:, 1] += (w * dz).sum(axis=1)
                for k, f in enumerate(fish):
                    f.x, f.z = P[k, 0], P[k, 1]
        t += dt
        # decisions
        for k, f in enumerate(fish):
            if t < f.next:
                continue
            ex, ez = math.cos(f.h), math.sin(f.h)
            hx, hz = f.x + 0.5 * L * ex, f.z + 0.5 * L * ez
            # the other fish are seen only within D_fish of the head (a school keeps its spacing, it does not
            # flee from it)
            others = np.array([[g.x, g.z, body["r_fish"]] for g in fish if g is not f
                               and math.hypot(g.x - hx, g.z - hz) < body["D_fish"] + body["r_fish"]]).reshape(-1, 3)
            circles = np.vstack([cyl, others]) if len(others) else cyl
            rays = [ray_lengths(hx, hz, math.cos(f.h + an), math.sin(f.h + an), box, circles, walls) for an in ANG]
            sv = {"t": t, "f": f.g["f"][0], "a": f.g["a"][0], "b": f.g["b"][0], "pos": [f.x, 0.0, f.z],
                  "axis": [ex, 0.0, ez], "lat": [-ez, 0.0, ex], "U": f.acc_U / max(1, f.acc_n),
                  "stride": f.st.get("kx", [0.5])[0], "target": [f.target[0], 0.0, f.target[1]],
                  "rays": rays, "ray_angles": ANG}
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
                rows.append((t, f.i, f.x, f.z, math.degrees(f.h), body["s"] * f.Ug, f.g["f"][0], f.g["a"][0], f.g["b"][0], f.mode))
            t_out += out_dt
    return fish, rows


def summarise(fish, rows, t0):
    out = []
    R = np.array(rows)
    for f in fish:
        m = (R[:, 1] == f.i) & (R[:, 0] > t0)
        d = np.hypot(R[m, 2] - f.target[0], R[m, 3] - f.target[1])
        out.append({"id": f.i, "mean": float(d.mean()) if len(d) else None, "max": float(d.max()) if len(d) else None,
                    "within_0.5": float((d < 0.5).mean()) if len(d) else None, "target": list(f.target)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--out", default="tracks.csv")
    ap.add_argument("--summary", default="")
    ap.add_argument("--after", type=float, default=12.0, help="the summary over t > after")
    a = ap.parse_args()
    case = json.load(open(a.case))
    fish, rows = run(case)
    with open(a.out, "w") as fo:
        fo.write("t id x z heading U f a b mode\n")
        for r in rows:
            fo.write("%.3f %d %.5f %.5f %.3f %.5f %.4f %.4f %.4f %d\n" % r)
    s = summarise(fish, rows, a.after)
    if a.summary:
        json.dump(s, open(a.summary, "w"), indent=1)
    md = [x["mean"] for x in s if x["mean"] is not None]
    wi = [x["within_0.5"] for x in s if x["within_0.5"] is not None]
    print("%d fish, %.1f s: mean distance to the target after %.0f s %.2f m, within 0.5 m %.0f %%" %
          (len(fish), case.get("T", 20.0), a.after, float(np.mean(md)), 100 * float(np.mean(wi))))


if __name__ == "__main__":
    main()
