#!/usr/bin/env python3
# External brain for a deforming 6DOF body of REEF3D::CFD (X 191 id port, with X 193 and X 192).
#
#   python3 brain_server.py --port 5555 [--policy rules|own|file.py] [--log brain_log.csv]
#
# REEF3D (rank 0) connects to 127.0.0.1:port and sends one line of JSON per decision (once per half
# tail beat); the server answers with one line of JSON:
#
#   hello  {"hello":1,"body":0,"L":..,"f0":..,"rho":..,"U_des":..,"b_max":..,"D":..,"th_esc":..,"t_start":..,"twoD":1,"s0":..,"k_turn":..,"tau_u":..,"tau_d":..,"n_u":..,"n_d":..,"s_sd":..}
#          -> {"ok":1,"policy":"rules"}
#   senses {"t":..,"f":..,"a":..,"b":..,            the time of the decision, the gait now
#           "pos":[x,y,z],"axis":[..],"lat":[..],   the CoG, the mean body axis and lateral axis in the plane
#                                                   of the undulation (unit vectors, since the last decision)
#           "U":..,"stride":..,                     the mean speed along the axis over the ground [m/s], the stride s of the own brain [L]
#           "target":[x,y,z],                       the target (X 192), if any
#           "rays":[..],"ray_angles":[..],          vision: distances [m] at the angles to the axis [rad]
#           "ll_p":[..],"ll_s":[..],"ll_side":[..], lateral line: mean pressure, s/L, side of each probe
#           "own":{"f","a","b","ts","mode","second","t2","f2","a2","b2","ts2"}}   the decision of the own brain
#          -> {"f":..,"a":..,"b":..,"ts":..}       the gait targets from now: frequency [Hz], amplitude
#                                                   factor, curvature bias, settling time [s] (0: half a beat)
#             optional "t2","f2","a2","b2","ts2":  a second target at time t2 (e.g. the end of a C-start)
#             optional "next":                     seconds to the next decision (default min(0.5/f, 0.5))
#   bye    {"bye":1}
#
# Policies: "own" answers the decision of the own brain (a test of the link: the run is the run of X 193);
# "rules" is the rule brain of REEF3D (6DOF_brain.cpp) in Python, to start from (the same decisions); a file defines
# decide(senses, state, params) -> dict (state: a dict kept between the decisions of the body).
#
# Architect: Hans Bihs
import argparse, json, math, socketserver, sys, threading, importlib.util

LOCK = threading.Lock()


# the constants of the rule brain (6DOF_brain.cpp)
TAU, KTURN, SLOW, ARRIVE, FRONT = 1.5, 1.3, 1.5, 0.25, 0.35
TAUU, S0, TFLOOR, PSIFLOOR, REAFF, BURST, TAUC0 = 0.8, 0.5, 2.5, 1.0, 0.5, 0.5, 2.0
CKNOWN, CDRIFT, BRAKE, NKNOWN = 0.2, 0.1, 2.5, 3
UESC = 0.3                     # an escape: the pressures in units of q at this speed at least [L/s]


def _kf_update(st, H, y, R, dt):
    """The Kalman filter on [s, c_p0, c_p1] (sixdof_brain::kf_update)."""
    x, P = st["kx"], st["kP"]
    ssd = st.get("ssd", 0.0)
    P[0][0] += 1.0e-4 * dt * ssd * ssd / 0.09 if ssd > 0.0 else 1.0e-4 * dt
    P[1][1] += 2.5e-3 * dt
    P[2][2] += 2.5e-3 * dt
    for m in range(2):
        Ph = [P[i][0] * H[m][0] + P[i][1] * H[m][1] + P[i][2] * H[m][2] for i in range(3)]
        S = H[m][0] * Ph[0] + H[m][1] * Ph[1] + H[m][2] * Ph[2] + R
        K = [Ph[i] / S for i in range(3)]
        e = y[m] - (H[m][0] * x[0] + H[m][1] * x[1] + H[m][2] * x[2])
        for i in range(3):
            x[i] += K[i] * e
        for i in range(3):
            for j in range(3):
                P[i][j] -= K[i] * Ph[j]
    x[0] = min(1.6 * st["s0"], max(0.5 * st["s0"], x[0]))


def rules(sv, st, pr):
    """The rule brain of REEF3D (sixdof_brain::decide after the senses), with its own state: the estimate of
    the stride and the current from the gait and the track, the velocity through the water for the target."""
    L, f0, rho = pr["L"], pr["f0"], pr["rho"]
    Udes, bmax, Dvis, the = pr["U_des"], (pr["b_max"] if pr["b_max"] > 0 else 1.0), pr["D"], pr["th_esc"]
    s0 = pr.get("s0", 0.0) or S0                 # X 200: the prior of the stride and the turn rate of the bias
    kturn = pr.get("k_turn", 0.0) or KTURN
    tauu = pr.get("tau_u", 0.0) or TAUU         # X 201: how fast the speed follows a faster and a slower gait
    taud = pr.get("tau_d", 0.0) or TAUU
    nu = max(0.0, pr.get("n_u", 0.0) or 0.0)     # X 202: the speed follows the gait in tail beats (replaces X 201)
    nd = (pr.get("n_d", 0.0) or 0.0) or nu
    ssd = max(0.0, pr.get("s_sd", 0.0) or 0.0)    # X 203: the spread of the stride about s0 (0: 0.3)
    t, f, a, b = sv["t"], sv["f"], sv["a"], sv["b"]
    ex, el, c = sv["axis"], sv["lat"], sv["pos"]
    nrm = 1 if pr["twoD"] else 2
    p0, p1 = 0, (2 if nrm == 1 else 1)
    fmx = 2.0 * f0

    # 1. the speed of the gait since the last decision; the estimate of the stride and the current
    if "kx" not in st:
        st["kx"], st["s0"] = [s0, 0.0, 0.0], s0
        st["kP"] = [[(ssd * ssd * s0 * s0 if ssd > 0.0 else 0.09 * s0 * s0), 0.0, 0.0], [0.0, 0.3 * 0.3, 0.0], [0.0, 0.0, 0.3 * 0.3]]
        st["ssd"] = ssd
        st["ug"], st["tdec"], st["pdec"], st["bdec"], st["Tprev"], st["vT"] = -1.0, t, list(c), b, None, [0.0, 0.0, 0.0]
        st["adec"] = a
        st["tauc"] = TAUC0
        st["refr"] = max(st.get("refr", -1.0), t + 1.0 / f0)   # no escape in the first tail beat (the start-up)
        st["E0"] = 0.0
        if nu > 0.0:
            # X 202: the body at rest over the ground at t = 0; the speed of the gait from rest over the ramp, the
            # water that passed the body at rest decaying with the coasting time (sixdof_brain::decide)
            n = max(1, int(math.ceil(t / 0.01))); h = t / n; u = 0.0
            for i in range(n):
                u += h * f0 / nu * (f0 * L * min(1.0, (i + 0.5) * h / max(t, 1.0e-9)) - u)
            st["ug"], st["E0"] = u, math.exp(-t / st["tauc"])
    dt = t - st["tdec"]
    if dt > 0.0:
        y = [(c[p0] - st["pdec"][p0]) / dt, (c[p1] - st["pdec"][p1]) / dt]
        if st["ug"] < 0.0:
            # the first interval: the speed of the gait from the track (no current assumed), no update
            st["ug"] = min(f * a * L, max(0.0, y[0] * ex[p0] + y[1] * ex[p1]) / st["kx"][0])
        elif min(a, st["adec"]) < 0.95:
            # gliding: no update of the estimate, the speed of the gait from the track
            E0 = st.get("E0", 0.0)
            Em = E0 * st["tauc"] / dt * (1.0 - math.exp(-dt / st["tauc"])) if E0 > 0.0 else 0.0
            st["E0"] = E0 * math.exp(-dt / st["tauc"])
            uw = (y[0] - (1.0 - Em) * st["kx"][1]) * ex[p0] + (y[1] - (1.0 - Em) * st["kx"][2]) * ex[p1]
            u0 = st["ug"] * st["kx"][0]
            # a coasting body does not move backwards through the water: drifting back faster than the current
            # read allows (straight, no brake, the tail still), the current along the axis is at least that drift
            if uw < 0.0 and max(a, st["adec"]) < 0.05 and max(abs(b), abs(st["bdec"])) < 0.05:
                st["kx"][1] += uw * ex[p0]
                st["kx"][2] += uw * ex[p1]
                uw = 0.0
            st["ug"] = max(0.0, uw) / st["kx"][0]
            # the coasting time from the slowing down of the glide (both decisions gliding)
            cm = math.sqrt(st["kx"][1] ** 2 + st["kx"][2] ** 2)
            if max(a, st["adec"]) < 0.05 and max(abs(b), abs(st["bdec"])) < 0.05 and uw > 3.0 * cm + 0.05 * L and uw < u0:
                tm = min(30.0, max(0.5, dt / math.log(u0 / uw)))
                st["tauc"] = 0.7 * st["tauc"] + 0.3 * tm
        else:
            fal = f * a * L
            if nu > 0.0:                                           # n_u, n_d tail beats of the gait (X 202)
                tu = (nu if fal >= st["ug"] else nd) / max(f, 0.1 * f0)
            else:
                tu = tauu if fal >= st["ug"] else taud             # up with tau_u, down with tau_d (X 201)
            e = math.exp(-dt / tu)
            ug1 = fal + (st["ug"] - fal) * e
            ugm = fal + (st["ug"] - fal) * (tu / dt) * (1.0 - e)
            st["ug"] = ug1
            bm = 0.5 * (b + st["bdec"])
            R = 0.07 * 0.07 * (1.0 + 4.0 * bm * bm)
            E0 = st.get("E0", 0.0)                                 # the water that passed the body at the start
            Em = E0 * st["tauc"] / dt * (1.0 - math.exp(-dt / st["tauc"])) if E0 > 0.0 else 0.0
            st["E0"] = E0 * math.exp(-dt / st["tauc"])
            H = [[ugm * ex[p0], 1.0 - Em * ex[p0] * ex[p0], -Em * ex[p1] * ex[p0]],
                 [ugm * ex[p1], -Em * ex[p0] * ex[p1], 1.0 - Em * ex[p1] * ex[p1]]]
            _kf_update(st, H, y, R, dt)
            st["nkf"] = st.get("nkf", 0) + 1
    s = st["kx"][0]
    cur = [0.0, 0.0, 0.0]
    cur[p0], cur[p1] = st["kx"][1], st["kx"][2]
    Uw = max(sv["U"] if st["ug"] < 0.0 else st["ug"] * s, 0.1 * L)
    q = 0.5 * rho * Uw * Uw

    # 2. the target: the wanted velocity over the ground, the velocity through the water that gives it
    has_t = "target" in sv
    T = sv["target"] if has_t else c
    r = [T[k] - c[k] for k in range(3)]
    r[nrm] = 0.0
    dist = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) if has_t else 0.0
    if has_t and st["Tprev"] is not None and dt > 0.0:
        al = min(1.0, dt / 1.0)
        st["vT"] = [(1.0 - al) * st["vT"][k] + al * (T[k] - st["Tprev"][k]) / dt for k in range(3)]
    st["Tprev"] = list(T) if has_t else None
    vw = [0.0, 0.0, 0.0]
    if has_t:
        sp = (Udes if Udes > 0.0 else f0 * s * L) * min(1.0, dist / (SLOW * L))
        if dist > SLOW * L:
            ca = (r[0] * ex[0] + r[1] * ex[1] + r[2] * ex[2]) / dist
            if ca > 0.0:
                sp = sp / max(ca, 0.5)
        vw = [st["vT"][k] + (sp * r[k] / dist if dist > 1.0e-9 else 0.0) for k in range(3)]
    w = [vw[k] - cur[k] for k in range(3)]
    w[nrm] = 0.0
    Ureq = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
    psi = math.atan2(w[0] * el[0] + w[1] * el[1] + w[2] * el[2], w[0] * ex[0] + w[1] * ex[1] + w[2] * ex[2]) if Ureq > 0.02 * L else 0.0

    # 3. the frequency from the speed through the water, at least what a turn-around needs; the bias
    ft, bt, mode = max(0.3 * f0, min(fmx, f)), 0.0, 0
    if has_t:
        Uturn = min(abs(psi) / TFLOOR * L / (abs(kturn) * bmax), fmx * s * L) if abs(psi) > PSIFLOOR else 0.0
        Ucmd = max(Ureq, Uturn)
        ft = max(0.3 * f0, min(fmx, Ucmd / (s * L)))
        bt = (psi / TAU) * L / (kturn * Uw)
        mode = 3 if dist < ARRIVE * L else 0
    bt = max(-bmax, min(bmax, bt))

    # the lateral line at the head (s/L < 0.35), on q through the water, the signal of the own turn removed
    dpl = phd = 0.0
    lp, ls, lsd = sv.get("ll_p", []), sv.get("ll_s", []), sv.get("ll_side", [])
    if lp:
        dsum = hsum = asum = 0.0
        nd = nh = 0
        for k in range(len(lp)):
            asum += lp[k]
            if ls[k] < FRONT:
                dsum += lsd[k] * lp[k]
                hsum += lp[k]
                nh += 1
                if lsd[k] > 0:
                    nd += 1
        if nd > 0 and nh > 0:
            dpl = (dsum / nd) / q + REAFF * 0.5 * (b + st["bdec"])
            phd = (hsum / nh - asum / len(lp)) / q

    # vision: the closest ray, the side with the longer rays (kept while the threat lasts)
    threat = side = 0.0
    rays, ang = sv.get("rays", []), sv.get("ray_angles", [])
    if rays and Dvis > 0.0 and len(rays) == len(ang):
        D = Dvis * L
        fp = fm = 0.0
        for k in range(len(rays)):
            dk = min(rays[k], D)
            threat = max(threat, (1.0 - dk / D) * (1.0 - dk / D))
            if ang[k] > 0.0:
                fp += dk
            if ang[k] < 0.0:
                fm += dk
        side = 1.0 if fp > fm else (-1.0 if fm > fp else (-1.0 if psi >= 0.0 else 1.0))
        if threat > 0.0 and st["side"] != 0.0:
            side = st["side"]
        st["side"] = side if threat > 0.0 else 0.0

    # glide (a = 0) in a slow current, nothing seen: early when the coast reaches the target ahead (the
    # current known), or at the target needing less than the slowest swimming speed; braking (mode 6)
    at = 1.0
    Umin = 0.3 * f0 * s * L
    # coasting: the distance it carries the fish (its speed through the water along the axis times the
    # coasting time); the target ahead on the axis (within 30 deg)
    dstop = max(0.0, st["ug"]) * s * st["tauc"]
    rx, ry = sum(r[k] * ex[k] for k in range(3)), sum(r[k] * el[k] for k in range(3))
    onaxis = has_t and rx > 0.0 and abs(math.atan2(ry, rx)) < 0.5236
    cm = math.sqrt(cur[p0] * cur[p0] + cur[p1] * cur[p1])
    calm = cm < Umin and math.sqrt(sum(v * v for v in st["vT"])) < Umin
    # an early glide needs the current known and its drift over the coasting time small
    # at the target: the current learnt over NKNOWN updates (a glide stops the estimate)
    cvar = st["kP"][1][1] + st["kP"][2][2] < (CKNOWN * L) ** 2
    known = cvar and cm * st["tauc"] < CDRIFT * L
    learnt = cvar and st.get("nkf", 0) >= NKNOWN
    if not has_t or threat > 0.0 or not calm:
        st["glide"] = False
    elif st.get("glide", False):
        # beyond BURST when the coast no longer reaches the target ahead; with the target behind it keeps
        # braking (the bent body turns it towards the target) up to twice BURST
        st["glide"] = not (dist > BURST * L and (dist > 2.0 * BURST * L if rx <= 0.0 else dstop < dist - BURST * L))
        if not st["glide"]:
            st["tburst"] = t + 1.0 / f0          # a burst lasts at least one tail beat
    elif t >= st.get("tburst", -1.0):
        st["glide"] = (known and onaxis and dstop >= dist - ARRIVE * L) or (learnt and dist < ARRIVE * L and Ureq < Umin)
    if st["glide"]:
        # the tail still: the slowest frequency, so that the bend is released slowly; braking keeps the
        # frequency of the turn, so that the bend builds quickly
        fturn = ft
        ft, at, bt, mode = 0.3 * f0, 0.0, 0.0, 5
        # braking: the coast would carry the fish past the target, away from it, or past it to the side (the
        # target ahead but off the axis): the body bent, to the side chosen at the start of the braking
        # (bending to and fro would swim); the fish turns towards the target while it coasts
        if BRAKE > 0.0 and ((dstop > dist + ARRIVE * L or (not onaxis and dist > ARRIVE * L)) if rx > 0.0
                            else dist > ARRIVE * L):
            if st.get("bside", 0.0) == 0.0:
                st["bside"] = 1.0 if psi * kturn >= 0.0 else -1.0
            ft, bt, mode = fturn, st["bside"] * BRAKE, 6
        else:
            st["bside"] = 0.0
    else:
        st["bside"] = 0.0

    st["tdec"], st["pdec"], st["bdec"], st["adec"] = t, list(c), bt, at
    st["est"] = (s, cur[0], cur[1], cur[2], Uw, st["tauc"])

    # 1. escape: a C-start away from the side of the higher pressure (not while coasting or braking)
    # (the pressures against q at UESC at least: at a low speed the own turns and the start-up are no stimulus)
    qesc = q / max(q, 0.5 * rho * (UESC * L) ** 2)
    if the > 0.0 and t >= st["refr"] and not (st["glide"] or a < 0.95) and (abs(dpl) * qesc > the or phd * qesc > the):
        sg = 1.0 if dpl >= 0.0 else -1.0
        st["refr"] = t + 2.0
        st["bdec"] = -2.5 * sg
        st["adec"] = 0.5
        st["glide"], st["bside"] = False, 0.0
        return {"f": f, "a": 0.5, "b": -2.5 * sg, "ts": 0.08,
                "t2": t + 0.15, "f2": min(1.5 * f, fmx), "a2": 1.0, "b2": 0.0, "ts2": 0.2, "next": 0.4, "mode": 2}

    # 2. avoidance: a turn to the free side, slower
    if threat > 0.0:
        wv = min(1.0, threat / 0.25)
        bt = (1.0 - wv) * bt + wv * side * bmax
        ft = max(0.3 * f0, ft * (1.0 - 0.3 * wv))
        st["bdec"] = bt
        if wv > 0.2:
            mode = 1

    return {"f": ft, "a": at, "b": bt, "ts": 0.0, "next": min(0.5 / ft, 0.5), "mode": mode}


def own(sv, st, pr):
    """The decision of the own brain of REEF3D, unchanged."""
    o = sv["own"]
    d = {"f": o["f"], "a": o["a"], "b": o["b"], "ts": o["ts"]}
    if o["second"]:
        d.update({"t2": o["t2"], "f2": o["f2"], "a2": o["a2"], "b2": o["b2"], "ts2": o["ts2"]})
    d["next"] = 0.4 if o["mode"] == 2 else min(0.5 / o["f"], 0.5)
    d["mode"] = o["mode"]
    return d


class Handler(socketserver.StreamRequestHandler):
    def handle(self):
        pr, st, log = None, {"side": 0.0, "refr": -1.0}, None
        policy = self.server.policy
        for line in self.rfile:
            m = json.loads(line)
            if "hello" in m:
                pr = m
                if self.server.logname:
                    log = open(self.server.logname.replace("{body}", str(m["body"])), "w")
                    log.write("t,x,z,dist,f,a,b,mode,own_f,own_b\n")
                out = {"ok": 1, "policy": self.server.policy_name}
            elif "bye" in m:
                break
            else:
                d = policy(m, st, pr)
                out = {k: v for k, v in d.items() if k != "mode"}
                if log:
                    T = m.get("target", m["pos"])
                    dist = math.hypot(T[0] - m["pos"][0], T[2 if pr["twoD"] else 1] - m["pos"][2 if pr["twoD"] else 1])
                    log.write("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%.6f,%.6f\n" % (m["t"], m["pos"][0], m["pos"][2], dist, d["f"], d["a"], d["b"],
                                                                                   d.get("mode", -1), m["own"]["f"], m["own"]["b"]))
                    log.flush()
            self.wfile.write((json.dumps(out, allow_nan=False, separators=(",", ":"), default=float) + "\n").encode())
            self.wfile.flush()
        if log:
            log.close()


class Server(socketserver.ThreadingTCPServer):
    allow_reuse_address = True
    daemon_threads = True


class OneServer(socketserver.TCPServer):
    allow_reuse_address = True


def main():
    ap = argparse.ArgumentParser(description=__doc__ if __doc__ else "")
    ap.add_argument("--port", type=int, default=5555)
    ap.add_argument("--policy", default="rules", help="rules, own, or a .py file with decide(senses, state, params)")
    ap.add_argument("--log", default="", help="CSV of the decisions, {body} replaced by the body number")
    ap.add_argument("--once", action="store_true", help="stop after the first connection closes")
    a = ap.parse_args()
    if a.policy == "rules":
        pol = rules
    elif a.policy == "own":
        pol = own
    else:
        spec = importlib.util.spec_from_file_location("policy", a.policy)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        pol = mod.decide
    srv = (OneServer if a.once else Server)(("127.0.0.1", a.port), Handler)
    srv.policy, srv.policy_name, srv.logname = pol, a.policy, a.log
    print("fish brain on 127.0.0.1:%d, policy %s" % (a.port, a.policy), flush=True)
    if a.once:
        srv.handle_request()
    else:
        srv.serve_forever()


if __name__ == "__main__":
    main()
