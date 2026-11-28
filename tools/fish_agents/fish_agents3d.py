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
# A school ("school"): no fish has its own orbit. Each fish follows its intention, a point that moves at the
# speed U over the ground in the direction the fish wants: round the cage (dir) within a band of radius and depth,
# away from neighbours closer than d_sep (separation), along the swimming direction of its neighbours within d_ali
# (alignment), towards neighbours within d_coh when it has few (cohesion); the fish see each other all round,
# within the distances. The point is held within leash of the fish, lead ahead of it at the start. Validation 27.
# A body "at_rest": 1 starts at rest over the ground (the flow it does not follow yet decays with tau_coast), as a
# resolved body does in REEF3D and as the start model of the brain (X 202, n_u) assumes.
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
#         or {"n": 60, "orbit": {"c": [cx, cy], "r": [r0, r1], "z": [z0, z1], "U": 0.15, "dir": 1}}
#         or {"n": 60, "school": {"c": [cx, cy], "r": [r0, r1], "z": [z0, z1], "U": 0.15, "dir": 1,
#                                 "d_sep": 0.7, "d_ali": 2.0, "d_coh": 3.0, "w_path": 1.0, "w_band": 1.0,
#                                 "w_sep": 1.5, "w_ali": 1.0, "w_coh": 0.5, "lead": 1.5, "leash": 2.5}}}
#    "flow": {..., "ramp": [[t, s], ...], "inside": 0.6}   (optional: the flow scaled by s(t), linear in t, and by
#         "inside" within the nets: the current damped by the net and the fish, 20 % by a bare net, up to 50 % in a
#         stocked cage, Hvas et al. 2021; 40 cm/s inside against 70 cm/s outside, Johansson et al. 2014)
#    "physiology": {"W": 2500, "T": 12} or {"ucrit": 1.0}, and the constants of tools/fish_physiology/physiology.py
#         (optional): the oxygen debt and the fatigue of each fish from its speed through the water; an exhausted fish
#         swims no faster than its aerobic limit until it has recovered
#    "stages": {"U_pref": [0.30, 0.08], "h": 0.1, "lead": 1.0, "burst_glide": 1, "u_bg": 0, "delta": 0.05,
#               "kappa": 0.05, "net_escape": 1, "t_pin": 2.0}   (optional, with "school" and "physiology"): the
#         behaviour of salmon in a cage as the current rises (Hvas et al. 2021, Fig. 1; Johansson et al. 2014):
#         (a) 0 circling: the school at its preferred speed;
#         (b, c) 1 standing on the current: a fish whose current exceeds its own preferred swimming speed
#              (U_pref, mean and SD over the fish [m/s], with the hysteresis h) leaves the circle and holds its
#              station facing the current, lead L upstream of where it stopped (at least lead L + 1 m inside the
#              net); when it falls below, it circles again;
#         (-) 2 burst and glide: standing in a current above its gait transition (U_a, or "u_bg" U_crit), the fish
#              no longer swims steadily ("burst_glide": 1): it beats at its fastest gait (2 f0) until its speed
#              through the water is (1 + delta) of the current, glides (the tail still) until it is (1 - delta) of
#              it, and holds its station on average (the thresholds rise by kappa per L the fish is behind it);
#              with "burst_glide": 0 it swims steadily and the stage only marks it as anaerobic;
#         (d) 3 exhausted: it cannot hold the current and falls back on the downstream net.
#         "net_escape": 1 (default): a fish that is not exhausted and is pressed against the net for t_pin s
#         turns towards its station and swims off at its fastest gait (the net is no free support); the contacts
#         are counted (welfare_metrics: the time on the net, the contacts per fish and hour).
#    "body": {"C_glide": 0.13}   (optional): the glide with the tail still by a drag law, dU/dt = -C_glide U^2/L,
#         instead of tau_coast (a real salmon at Re 5e5: C_f 0.0066 on a wetted area 0.47 L^2, about 0.13)
#    with "physiology" the cost of swimming is that of the gait: the demand at the speed the gait would hold,
#         s f a L (no faster than the fish can swim), and SMR with the tail still (a glide, a brake)
# Output: the tracks (t id x y z heading_deg U f a b mode) every out_dt, with "physiology" also (stage F MO2 U_crit
# u_cur): the stage, the fatigue (debt over capacity), MO2 [mg O2/kg/h], the fish's U_crit and the current at it;
# the summary per fish (distance to its target, mean frequency); for a school also school_metrics() (spacing,
# polarisation, milling, the net, contacts) and with "physiology" welfare_metrics() (the time above 60 and 80 %
# of U_crit, the shares in each stage, the fish exhausted).
#
# Architect: Hans Bihs
import argparse, glob, json, math, os, random, re, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fish_physiology"))
import fish_agents as fa                                           # noqa: E402
import physiology as phy                                           # noqa: E402
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
                   "n_u": brain.get("n_u", 0.0), "n_d": brain.get("n_d", 0.0), "s_sd": brain.get("s_sd", 0.0),
                   "tau_c": brain.get("tau_c", 0.0)}
        self.acc_U = 0.0; self.acc_n = 0
        self.mode = -1
        self.vx = self.vy = self.vz = 0.0     # the velocity over the ground of the last step
        self.contacts = 0
        self.stage = 0                        # the cage stage (stages): 0 circling, 1 standing, 2 anaerobic, 3 exhausted
        self.station = None
        self.bg = None                        # burst and glide: "burst" | "glide" | None (steady)
        self.pin = 0.0                        # the time pressed against the net
        self.escape = None                    # the end of an escape from the net
        self.on_net = False; self.net_t = 0.0; self.net_n = 0
        self.Uc = 0.0


SCHOOL = {"U": 0.15, "dir": 1, "d_sep": 0.7, "d_ali": 2.0, "d_coh": 3.0, "w_path": 1.0, "w_band": 1.0,
          "w_sep": 1.5, "w_ali": 1.0, "w_coh": 0.5, "lead": 1.5, "leash": 2.5}


def school_step(fish, sc, dt):
    """move the intention points of a school: each at the speed U over the ground in the wanted direction
    (round the cage within the band, separation, alignment, cohesion with the neighbours seen), held within leash
    of its fish"""
    P = np.array([[f.x, f.y, f.z] for f in fish]); V = np.array([[f.vx, f.vy, f.vz] for f in fish])
    D = P[None, :, :] - P[:, None, :]                         # D[i, j] = p_j - p_i
    dd = np.linalg.norm(D, axis=2); np.fill_diagonal(dd, np.inf)
    cx, cy = sc["c"]
    for i, f in enumerate(fish):
        if f.stage > 0:                       # standing on the current, or exhausted: its own station (stage_step)
            continue
        q = np.array(f.goal)
        # the path: round the cage at the radius and depth band (springs only outside the band)
        rx, ry = q[0] - cx, q[1] - cy
        rr = max(math.hypot(rx, ry), 1e-6)
        e_t = np.array([-ry, rx, 0.0]) / rr * sc["dir"]
        e_r = np.array([rx, ry, 0.0]) / rr
        band = np.zeros(3)
        if rr < sc["r"][0]:
            band += e_r * (sc["r"][0] - rr)
        elif rr > sc["r"][1]:
            band -= e_r * (rr - sc["r"][1])
        if q[2] < sc["z"][0]:
            band[2] += sc["z"][0] - q[2]
        elif q[2] > sc["z"][1]:
            band[2] -= q[2] - sc["z"][1]
        w = sc["w_path"] * e_t + sc["w_band"] * band
        # the neighbours
        d_i = dd[i]
        ns = d_i < sc["d_sep"]
        if ns.any():
            sep = -(D[i, ns] / (d_i[ns, None] ** 2)).sum(axis=0) * sc["d_sep"]
            w += sc["w_sep"] * sep
        na = d_i < sc["d_ali"]
        if na.any() and sc["w_ali"] > 0.0:
            vm = V[na].mean(axis=0); nv = np.linalg.norm(vm)
            if nv > 1e-6:
                w += sc["w_ali"] * vm / nv
        nc = d_i < sc["d_coh"]
        if nc.sum() < 3 and sc["w_coh"] > 0.0:                    # few neighbours: towards those seen
            nc2 = d_i < 2.0 * sc["d_coh"]
            if nc2.any():
                dc = D[i, nc2].mean(axis=0); nd = np.linalg.norm(dc)
                if nd > 1e-6:
                    w += sc["w_coh"] * dc / nd
        nw = np.linalg.norm(w)
        e = w / nw if nw > 1e-9 else e_t
        q = q + dt * sc["U"] * e
        # held within leash of the fish
        p = P[i]; dq = q - p; ndq = np.linalg.norm(dq)
        if ndq > sc["leash"]:
            q = p + dq * sc["leash"] / ndq
        f.goal = list(q)


def school_metrics(rows, net=None, t0=60.0, r_fish=0.15):
    """the school after t0: nearest-neighbour distance, polarisation |mean e|, milling (mean tangential unit
    velocity round the cage centre), radius and depth, the share within 0.5 m of the net, contacts (pairs closer
    than 2 r_fish per fish and minute)"""
    R = np.array(rows)
    ts = np.unique(R[:, 0]); ts = ts[ts > t0]
    nnd, pol, mil, rad, dep, near, cont = [], [], [], [], [], [], 0
    for t in ts:
        m = R[:, 0] == t
        P = R[m][:, 2:5]; h = np.radians(R[m][:, 5])
        D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2); np.fill_diagonal(D, np.inf)
        nnd.extend(D.min(axis=1)); cont += int((D < 2 * r_fish).sum() // 2)
        e = np.stack([np.cos(h), np.sin(h)], axis=1)
        pol.append(np.linalg.norm(e.mean(axis=0)))
        if net is not None:
            rx, ry = P[:, 0] - net[0], P[:, 1] - net[1]; rr = np.hypot(rx, ry)
            mil.append(np.mean((-ry * e[:, 0] + rx * e[:, 1]) / np.maximum(rr, 1e-9)))
            rad.extend(rr); near.extend(rr > net[2] - 0.5)
        dep.extend(P[:, 2])
    nf = len(np.unique(R[:, 1])); dur = (ts[-1] - ts[0]) if len(ts) > 1 else 1.0
    return {"nnd_mean": float(np.mean(nnd)), "nnd_p10": float(np.percentile(nnd, 10)), "nnd_p90": float(np.percentile(nnd, 90)),
            "polarisation": float(np.mean(pol)), "milling": float(np.mean(mil)) if mil else None,
            "radius_mean": float(np.mean(rad)) if rad else None, "radius_p10": float(np.percentile(rad, 10)) if rad else None,
            "radius_p90": float(np.percentile(rad, 90)) if rad else None, "depth_mean": float(np.mean(dep)),
            "near_net": float(np.mean(near)) if near else None, "contacts_per_fish_min": float(cont / nf / (dur / 60.0))}


def body_step3(Ug, fg, ag, bg, t, dt, body):
    """fish_agents.body_step, with the glide (tail still, body straight) by a drag law when body["C_glide"] > 0"""
    if body.get("C_glide", 0.0) > 0.0 and ag < 0.5 and abs(bg) <= 0.5:
        return Ug - dt * body["C_glide"] * body["s"] * Ug * Ug / body["L"], body["k_brake"]
    return fa.body_step(Ug, fg, ag, bg, t, dt, body)


def stage_step(fish, st, ph, cur, U_pref, L):
    """the cage stages: circling (0), standing on the current (1), on anaerobic power (2), exhausted (3); cur the
    horizontal current at each fish (2 x n), ph the physiology of the fish"""
    h = st["h"]
    for k, f in enumerate(fish):
        c = math.hypot(cur[0][k], cur[1][k])
        if ph.exhausted[k]:
            new = 3
        elif f.stage == 0:
            new = 1 if c > U_pref[k] * (1.0 + h) else 0
        else:
            new = 0 if c < U_pref[k] * (1.0 - h) else 1
        if new == 0 and f.stage > 0:                         # back to the circle: the intention ahead of the fish
            f.goal = [f.x + st["lead"] * L * math.cos(f.h), f.y + st["lead"] * L * math.sin(f.h), f.z]
            f.station = None
        if new > 0:
            if f.station is None:
                # where it stopped, but at least lead L + 1 m inside the net of its cage (not leaning on the net)
                f.station = [f.x, f.y, f.z]
                if st.get("net") is not None:
                    nx_, ny_, nr_ = st["net"][:3]
                    dr = math.hypot(f.x - nx_, f.y - ny_); rmax = max(0.0, nr_ - st["lead"] * L - 1.0)
                    if dr > rmax:
                        f.station = [nx_ + (f.x - nx_) * rmax / dr, ny_ + (f.y - ny_) * rmax / dr, f.z]
            if c > 1e-6:                                     # the station, lead L upstream
                f.goal = [f.station[0] - st["lead"] * L * cur[0][k] / c, f.station[1] - st["lead"] * L * cur[1][k] / c, f.station[2]]
            if st.get("burst_glide", 1):
                if new == 1 and c > st["u_bg_abs"][k] * (1.0 + 0.5 * h):
                    new = 2
                elif new == 1 and f.stage == 2 and c > st["u_bg_abs"][k] * (1.0 - 0.5 * h):
                    new = 2
            elif new == 1 and f.Uw > ph.u_aerobic()[k]:
                new = 2
        if new != 2 or not st.get("burst_glide", 1):
            f.bg = None
        elif f.bg is None:
            f.bg = "burst"
        f.stage = new


def gait_override(f, k, t, st, cur, body, L):
    """the gait targets that replace the brain's: the escape from the net, the burst and glide"""
    f0 = body["f0"]
    # the escape: pressed against the net for t_pin s, not exhausted: turn towards the station and swim off
    if st.get("net_escape", 1) and f.stage in (0, 1, 2) and st.get("net") is not None:
        if f.escape is None and f.pin >= st["t_pin"]:
            f.escape = t + 10.0
        if f.escape is not None:
            nx_, ny_ = st["net"][0], st["net"][1]
            gx, gy = (f.goal[0], f.goal[1]) if f.station is not None else (nx_, ny_)
            want = math.atan2(gy - f.y, gx - f.x)
            err = math.atan2(math.sin(want - f.h), math.cos(want - f.h))
            off = st["net"][2] - math.hypot(f.x - nx_, f.y - ny_)
            if t > f.escape or (abs(err) < 0.5 and off > 0.5):
                f.escape = None
            else:
                f.g["f"][2], f.g["f"][3] = 2.0 * f0, 0.15
                f.g["a"][2], f.g["a"][3] = 1.0, 0.15
                f.g["b"][2], f.g["b"][3] = max(-1.0, min(1.0, 2.0 * err)), 0.15
                f.second = None
                return
    # the burst and glide
    if f.bg is None:
        return
    c = math.hypot(cur[0][k], cur[1][k])
    e = 0.0
    if f.station is not None and c > 1e-6:                    # how far behind its station the fish is [L]
        e = ((f.x - f.station[0]) * cur[0][k] + (f.y - f.station[1]) * cur[1][k]) / c / L
    hi = c * (1.0 + st["delta"] + st["kappa"] * e)
    lo = c * (1.0 - st["delta"] + st["kappa"] * e)
    if f.bg == "burst" and f.Uw >= hi:
        f.bg = "glide"
    elif f.bg == "glide" and f.Uw <= lo:
        f.bg = "burst"
    b = max(-0.45, min(0.45, f.g["b"][2]))                    # the brain's steering, mild enough to glide
    if f.bg == "burst":
        f.g["f"][2], f.g["a"][2] = 2.0 * f0, 1.0
    else:
        f.g["a"][2] = 0.0
    f.g["f"][3] = f.g["a"][3] = 0.15
    f.g["b"][2] = b
    f.second = None


def welfare_metrics(rows, t0=0.0, net=None, fish=None):
    """with "physiology": per fish the share of the time above 60 and 80 % of its U_crit (the speed through the
    water), the shares of the fish in each stage over time, the fish exhausted at least once; with the net, the share
    of the time a fish that is not exhausted leans on the net (within 0.2 m of it, slower than 0.6 of its current:
    the net holds it, which costs it nothing in the model)"""
    R = np.array(rows)
    R = R[R[:, 0] > t0]
    lean = None
    if net is not None:
        rr = np.hypot(R[:, 2] - net[0], R[:, 3] - net[1])
        lean = float(((rr > net[2] - 0.3) & (R[:, 6] < 0.6 * R[:, 15]) & (R[:, 15] > 0.1) & (R[:, 11] < 3)).mean())
    ids = np.unique(R[:, 1])
    u = R[:, 6] / R[:, 14]
    a60 = [float((u[R[:, 1] == i] > 0.6).mean()) for i in ids]
    a80 = [float((u[R[:, 1] == i] > 0.8).mean()) for i in ids]
    ex = [bool((R[R[:, 1] == i, 11] == 3).any()) for i in ids]
    share = {k: float((R[:, 11] == k).mean()) for k in range(4)}
    out = {"above60_mean": float(np.mean(a60)), "above80_mean": float(np.mean(a80)), "exhausted_fish": float(np.mean(ex)),
           "stage_share": share, "F_max": float(R[:, 12].max()), "MO2_mean": float(R[:, 13].mean()), "leaning_on_net": lean}
    if fish is not None:
        dur = max(R[:, 0].max() - R[:, 0].min(), 1e-9)
        out["net_time_s_mean"] = float(np.mean([f.net_t for f in fish]))
        out["net_contacts_per_fish_h"] = float(np.mean([f.net_n for f in fish]) / (dur / 3600.0))
    return out


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
    ramp = np.array(fl["ramp"], float) if "ramp" in fl else None
    inside = fl.get("inside", 1.0)
    t_start = case.get("t_start", max(1.0, 0.5 / f0))              # the brain starts at the end of the ramp

    fish = []
    spec = case["fish"]
    if isinstance(spec, dict) and "school" in spec:
        sc = dict(SCHOOL, **spec["school"])
        for i in range(spec["n"]):
            r = rnd.uniform(*sc["r"]); z = rnd.uniform(*sc["z"]); ph = rnd.uniform(0.0, 2.0 * math.pi)
            x0, y0 = sc["c"][0] + r * math.cos(ph), sc["c"][1] + r * math.sin(ph)
            h = ph + sc["dir"] * 0.5 * math.pi + rnd.uniform(-0.5, 0.5)
            f = Fish3(i, x0, y0, z, h, None, body, brain, t_start + rnd.uniform(0, 0.5 / f0))
            f.goal = [x0 + sc["lead"] * math.cos(h), y0 + sc["lead"] * math.sin(h), z]
            f.target = (lambda t, f=f: tuple(f.goal))
            fish.append(f)
    elif isinstance(spec, dict):
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
    school = dict(SCHOOL, **spec["school"]) if isinstance(spec, dict) and "school" in spec else None
    at_rest = body.get("at_rest", 0)
    # the physiology and the cage stages
    ph, stg, U_pref = None, None, None
    if "physiology" in case:
        pc = dict(case["physiology"])
        uc = pc.pop("ucrit", None)
        if uc is None:
            uc = float(phy.ucrit(pc.pop("W", 1000.0), pc.pop("T", 10.0)))
        else:
            pc.pop("W", None); pc.pop("T", None)
        ph = phy.Physiology(n, uc, seed=case.get("seed", 1), **pc)
        if "stages" in case and school is not None:
            stg = dict({"U_pref": [0.3, 0.08], "h": 0.1, "lead": 1.0, "burst_glide": 1, "u_bg": None, "delta": 0.05,
                        "kappa": 0.05, "net_escape": 1, "t_pin": 2.0}, **case["stages"])
            stg["net"] = nets[0] if nets else None
            stg["u_bg_abs"] = (ph.ucrit * stg["u_bg"]) if stg["u_bg"] else ph.u_aerobic()
            rs = np.random.default_rng(case.get("seed", 1) + 101)
            U_pref = np.clip(stg["U_pref"][0] + stg["U_pref"][1] * rs.standard_normal(n), 0.05, None)
    phst = {"F": np.zeros(n), "MO2": np.zeros(n)}
    for f in fish:
        f.Uw = 0.0
    cur = None
    while t < T - 1e-9:
        if school is not None:
            school_step(fish, school, dt)
        X = np.array([f.x for f in fish]); Y = np.array([f.y for f in fish]); Z = np.array([f.z for f in fish])
        uf, vf, wf = flow(t, X, Y, Z)
        if ramp is not None or inside != 1.0:
            sc_ = np.full(n, float(np.interp(t, ramp[:, 0], ramp[:, 1])) if ramp is not None else 1.0)
            for c in nets:
                sc_ = np.where(np.hypot(X - c[0], Y - c[1]) < c[2], sc_ * inside, sc_)
            uf, vf, wf = uf * sc_, vf * sc_, wf * sc_
        cur = (uf, vf)
        umax = ph.u_max() if ph is not None else None
        for k, f in enumerate(fish):
            if stg is not None:
                gait_override(f, k, t, stg, cur, body, L)
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
            f.Ug, kt = body_step3(f.Ug, fg, ag, bg, t, dt, body)
            if umax is not None and body["s"] * f.Ug > umax[k]:     # exhausted: no faster than its aerobic limit
                f.Ug = umax[k] / body["s"]
            Uw = body["s"] * f.Ug
            f.Uw = Uw
            # the cost of swimming: the speed the gait would hold (SMR with the tail still)
            f.Uc = 0.0 if ag < 0.5 else body["s"] * fg * ag * L
            if umax is not None:
                f.Uc = min(f.Uc, umax[k])
                if umax[k] < 1e9 and f.on_net:                     # exhausted on the net: it rests there
                    f.Uc = 0.0
            f.h += dt * kt * bg * Uw / L
            ex, ey = math.cos(f.h), math.sin(f.h)
            wz = max(-depth["w_max"], min(depth["w_max"], (f.target(t)[2] - f.z) / depth["tau"]))
            # (at rest over the ground at t = 0: the flow it does not follow yet decays with tau_coast)
            fr = 1.0 - math.exp(-t / body["tau_coast"]) if at_rest else 1.0
            nx, ny, nz = f.x + dt * (Uw * ex + fr * uf[k]), f.y + dt * (Uw * ey + fr * vf[k]), f.z + dt * (wz + fr * wf[k])
            # held in the box, outside the solids, inside the net of its cage
            nx = min(max(nx, box[0] + 0.05), box[1] - 0.05); ny = min(max(ny, box[2] + 0.05), box[3] - 0.05)
            nz = min(max(nz, box[4] + 0.05), box[5] - 0.05)
            for c in cyl:
                dd = math.hypot(nx - c[0], ny - c[1])
                if c[3] <= nz <= c[4] and dd < c[2] + 0.02:
                    s_ = (c[2] + 0.02) / max(dd, 1e-9); nx, ny = c[0] + (nx - c[0]) * s_, c[1] + (ny - c[1]) * s_
            held = False
            for c in nets:
                if math.hypot(f.x - c[0], f.y - c[1]) < c[2]:          # in this cage
                    dd = math.hypot(nx - c[0], ny - c[1])
                    if dd > c[2] - 0.1:
                        s_ = (c[2] - 0.1) / dd; nx, ny = c[0] + (nx - c[0]) * s_, c[1] + (ny - c[1]) * s_
                        held = True
                    nz = min(max(nz, c[3] + 0.1), c[4] - 0.1)
            if stg is not None:
                if held and not f.on_net:
                    f.net_n += 1
                f.on_net = held
                f.pin = f.pin + dt if held else 0.0
                if held:
                    f.net_t += dt
            f.acc_U += ((nx - f.x) * ex + (ny - f.y) * ey) / dt; f.acc_n += 1
            f.vx, f.vy, f.vz = (nx - f.x) / dt, (ny - f.y) / dt, (nz - f.z) / dt
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
        if ph is not None:
            phst = ph.step(np.array([f.Uc for f in fish]), dt)
            if stg is not None:
                stage_step(fish, stg, ph, cur, U_pref, L)
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
            f.next = t + d.get("next", min(0.5 / d["f"], 0.5 / f.pr["f0"]))
            f.mode = d.get("mode", 0)
            f.acc_U = 0.0; f.acc_n = 0
        if t >= t_out - 1e-9:
            for k, f in enumerate(fish):
                r = (t, f.i, f.x, f.y, f.z, math.degrees(f.h), body["s"] * f.Ug, f.g["f"][0], f.g["a"][0], f.g["b"][0], f.mode)
                if ph is not None:
                    r += (f.stage, float(phst["F"][k]), float(phst["MO2"][k]), float(ph.ucrit[k]),
                          float(math.hypot(cur[0][k], cur[1][k])) if cur is not None else 0.0)
                rows.append(r)
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
        ext = len(rows[0]) > 11
        fo.write("t id x y z heading U f a b mode" + (" stage F MO2 U_crit u_cur" if ext else "") + "\n")
        for r in rows:
            fo.write("%.3f %d %.4f %.4f %.4f %.3f %.5f %.4f %.4f %.4f %d" % r[:11] + (" %d %.4f %.1f %.4f %.4f" % r[11:] if ext else "") + "\n")
    spec = case["fish"]
    if isinstance(spec, dict) and "school" in spec:
        net = case["nets"][0] if case.get("nets") else None
        m = school_metrics(rows, net, a.after, dict(fa.BODY, **case.get("body", {}))["r_fish"])
        R = np.array(rows); fm = float(R[R[:, 0] > a.after, 7].mean())
        m["f_mean"] = fm
        if "physiology" in case:
            m["welfare"] = welfare_metrics(rows, a.after, net, fish)
            print("welfare after %.0f s: %s" % (a.after, json.dumps(m["welfare"])))
        if a.summary:
            json.dump(m, open(a.summary, "w"), indent=1)
        print("%d fish, %.1f s, after %.0f s: nearest neighbour %.2f m (10-90 %%: %.2f-%.2f), polarisation %.2f, milling %s, "
              "radius %s m, depth %.2f m, near the net %s, contacts %.2f per fish and minute, mean frequency %.2f Hz" %
              (len(fish), case.get("T", 60.0), a.after, m["nnd_mean"], m["nnd_p10"], m["nnd_p90"], m["polarisation"],
               "%.2f" % m["milling"] if m["milling"] is not None else "-", "%.2f" % m["radius_mean"] if m["radius_mean"] else "-",
               m["depth_mean"], "%.0f %%" % (100 * m["near_net"]) if m["near_net"] is not None else "-", m["contacts_per_fish_min"], fm))
        return
    s = summarise(fish, rows, a.after)
    if a.summary:
        json.dump(s, open(a.summary, "w"), indent=1)
    print("%d fish, %.1f s: mean distance to the target after %.0f s %.2f m, within 1 m %.0f %%, mean frequency %.2f Hz" %
          (len(fish), case.get("T", 60.0), a.after, np.mean([x["mean"] for x in s]), 100 * np.mean([x["within_1"] for x in s]),
           np.mean([x["f_mean"] for x in s])))


if __name__ == "__main__":
    main()
