#!/usr/bin/env python3
# External brain for a deforming 6DOF body of REEF3D::CFD (X 191 id port, with X 193 and X 192).
#
#   python3 brain_server.py --port 5555 [--policy rules|own|file.py] [--log brain_log.csv]
#
# REEF3D (rank 0) connects to 127.0.0.1:port and sends one line of JSON per decision (once per half
# tail beat); the server answers with one line of JSON:
#
#   hello  {"hello":1,"body":0,"L":..,"f0":..,"rho":..,"U_des":..,"b_max":..,"D":..,"th_esc":..,"t_start":..,"twoD":1}
#          -> {"ok":1,"policy":"rules"}
#   senses {"t":..,"f":..,"a":..,"b":..,            the time of the decision, the gait now
#           "pos":[x,y,z],"axis":[..],"lat":[..],   the CoG, the mean body axis and lateral axis in the plane
#                                                   of the undulation (unit vectors, since the last decision)
#           "U":..,"stride":..,                     the mean speed along the axis [m/s], the stride U/(f L) of the own brain
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
# "rules" is the rule brain of REEF3D (6DOF_brain.cpp) in Python, to start from; a file defines
# decide(senses, state, params) -> dict (state: a dict kept between the decisions of the body).
#
# Architect: Hans Bihs
import argparse, json, math, socketserver, sys, threading, importlib.util

LOCK = threading.Lock()


def rules(sv, st, pr):
    """The rule brain of REEF3D (sixdof_brain::decide after the senses), with its own state."""
    L, f0, rho = pr["L"], pr["f0"], pr["rho"]
    Udes, bmax, Dvis, the = pr["U_des"], (pr["b_max"] if pr["b_max"] > 0 else 1.0), pr["D"], pr["th_esc"]
    t, f, a, b = sv["t"], sv["f"], sv["a"], sv["b"]
    ex, el, c = sv["axis"], sv["lat"], sv["pos"]
    U = sv["U"]
    Uf = max(U, 0.1 * L)
    q = 0.5 * rho * Uf * Uf
    fmax = 2.0 * f0
    nrm = 1 if pr["twoD"] else 2

    # the stride U/f while swimming straight at the full amplitude, after the start
    if abs(b) < 0.15 and abs(a - 1.0) < 0.05 and t > pr["t_start"] + 2.0 / f0 and U > 0.0 and f > 0.0:
        s = U / (f * L)
        st["stride"] = s if st["nstride"] == 0 else 0.8 * st["stride"] + 0.2 * s
        st["nstride"] += 1

    # the target: distance and bearing (> 0: towards the lateral axis)
    has_t = "target" in sv
    T = sv["target"] if has_t else c
    r = [T[k] - c[k] for k in range(3)]
    r[nrm] = 0.0
    dist = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) if has_t else 0.0
    psi = math.atan2(sum(r[k] * el[k] for k in range(3)), sum(r[k] * ex[k] for k in range(3))) if has_t else 0.0

    # the lateral line at the head (s/L < 0.35)
    dpl = phd = 0.0
    lp, ls, lsd = sv.get("ll_p", []), sv.get("ll_s", []), sv.get("ll_side", [])
    if lp:
        dsum = hsum = asum = 0.0
        nd = nh = 0
        for k in range(len(lp)):
            asum += lp[k]
            if ls[k] < 0.35:
                dsum += lsd[k] * lp[k]
                hsum += lp[k]
                nh += 1
                if lsd[k] > 0:
                    nd += 1
        if nd > 0 and nh > 0:
            dpl = (dsum / nd) / q
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

    # 1. escape: a C-start away from the side of the higher pressure
    if the > 0.0 and t >= st["refr"] and (abs(dpl) > the or phd > the):
        sg = 1.0 if dpl >= 0.0 else -1.0
        st["refr"] = t + 2.0
        return {"f": f, "a": 0.5, "b": -2.5 * sg, "ts": 0.08,
                "t2": t + 0.15, "f2": min(1.5 * f, fmax), "a2": 1.0, "b2": 0.0, "ts2": 0.2, "next": 0.4, "mode": 2}

    # 3. the target: the turn rate psi/1.5 s, the frequency from the closing speed
    ft, bt, mode = f, 0.0, 0
    if has_t:
        bt = (psi / 1.5) * L / (1.41 * Uf)
        if Udes > 0.0:
            Ud = Udes * min(1.0, dist / (1.5 * L))
            stl = max(st["stride"], 0.1) * L
            if st["dprev"] >= 0.0 and t > st["tprev"] and st["fctl"] > 0.0:
                ft = st["fctl"] + 0.3 * (Ud - (st["dprev"] - dist) / (t - st["tprev"])) / stl
            else:
                ft = Ud / stl
            if abs(psi) > 1.0472:
                ft = max(ft, 0.8 * f0)
            ft = max(0.3 * f0, min(fmax, ft))
            st["fctl"] = ft
        st["dprev"], st["tprev"] = dist, t
        mode = 3 if dist < 0.25 * L else 0

    bt = max(-bmax, min(bmax, bt))

    # 2. avoidance: a turn to the free side, slower
    if threat > 0.0:
        w = min(1.0, threat / 0.25)
        bt = (1.0 - w) * bt + w * side * bmax
        ft = max(0.3 * f0, ft * (1.0 - 0.3 * w))
        if w > 0.2:
            mode = 1

    return {"f": ft, "a": 1.0, "b": bt, "ts": 0.0, "next": min(0.5 / ft, 0.5), "mode": mode}


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
        pr, st, log = None, {"stride": 0.55, "nstride": 0, "dprev": -1.0, "tprev": 0.0, "side": 0.0, "fctl": -1.0, "refr": -1.0}, None
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
