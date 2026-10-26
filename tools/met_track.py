#!/usr/bin/env python3
# Architect: Hans Bihs
"""
REEF3D parametric tropical cyclone vortex (W 200 3): best track -> fitted vortex track met-vortex.trk.

The equivalent of ADCIRC's ASWIP for REEF3D: reads a best track, removes the storm translation from the observed
winds, fits the vortex parameters and writes them per record; REEF3D evaluates the vortex per cell and time step
(src/met_vortex.cpp, the same formulas as evaluate() here, which is the reference of the unit test).

  fit        best track -> met-vortex.trk
             ATCF b-deck (bal142018.dat: one line per time and isotach 34/50/64 kt, radii NE SE SW NW) or
             IBTrACS CSV (ISO_TIME, LAT, LON, USA_WIND, USA_PRES, USA_RMW, USA_POCI, USA_R34_NE ... USA_R64_NW)
  synthetic  an idealised straight track in ATCF format (tests, validation)
  evaluate   wind (east, north) and pressure of the vortex at points (reference tables for the C++ code)
  info       the records of a met-vortex.trk

Vortex models (W 215; all fitted here, the model is chosen in ctrl.txt):
  1 Holland (1980): symmetric gradient wind, B = rho_a e Vg^2 / dp, Rmax = RMW of the track
      Vg(r) = sqrt(B dp/rho_a (Rm/r)^B exp(-(Rm/r)^B) + (r f/2)^2) - r f/2,   p = pc + dp exp(-(Rm/r)^B)
  2 GAHM, generalized asymmetric Holland model (Gao 2018, ADCIRC NWS 20): per quadrant (NE SE SW NW) and isotach
      (34, 50, 64 kt) Rmax, Bg, phi such that the profile passes through the isotach,
      Vg(r) = sqrt(Vg^2 (1+1/Ro) exp(phi (1-(Rm/r)^Bg)) (Rm/r)^Bg + (r f/2)^2) - r f/2,  Ro = Vg/(f Rm)
      Bg = (Vg^2 + Vg Rm f) rho_a e^phi / (phi dp),  phi = 1 + Vg Rm f / (Bg (Vg^2 + Vg Rm f))   (iterated)
      p = pc + dp exp(-phi (Rm/r)^Bg); between isotachs linear in r, between quadrants inverse-square in azimuth
  3 Holland, Belanger and Fritz (2010): surface wind V = Vms [(Rm/r)^bs exp(1 - (Rm/r)^bs)]^x, bs = rho_a e Vms^2/dp,
      x = 0.5 inside Rm, rising linearly to xn at rn (the mean 34 kt radius) outside; xn from V(rn) = 34 kt
Surface wind (models 1, 2): Wrf Vg (W 216, 0.9), turned inwards by the inflow angle of ADCIRC (W 218 1: 10 r/Rm inside
Rm, 10 + 75 (r/Rm - 1) to 1.2 Rm, then 25 deg), plus the translation scaled by Vg/Vmax (W 219 1), times the 1-min to
10-min factor (W 217, 0.8928). The fit removes the translation the same way: Vmax_g = (Vmax - |Vt|)/Wrf, and at an
isotach the gradient wind Vg with |Wrf Vg d + (Vg/Vmax_g) Vt| = V_isotach (d: the surface direction of the quadrant).

File met-vortex.trk (ASCII, SI units):
  REEF3D-MET-TRACK 1 nrec rho_a
  '#' lines: comments
  per record 64 numbers:
    YYYYMMDD.HHMMSS lon lat ut vt pc pn vmaxg rmw B80 vms bs xn rn Bg0 phi0
    then for the quadrants NE SE SW NW and the isotachs 34 50 64 kt: R Rm Bg phi  (R = 0: no isotach)
  lon, lat [deg] of the centre, ut, vt translation (east, north) [m/s], pc, pn central and outer pressure [Pa],
  vmaxg gradient Vmax without translation [m/s], rmw [m], B80 Holland B, vms surface Vmax without translation,
  bs Holland 2010 b, xn, rn the outer exponent and its radius [m], Bg0 phi0 GAHM with Rm = rmw (quadrants without isotachs)

Examples
  python3 met_track.py fit bal142018.dat -o met-vortex.trk --start 20181009.000000 --end 20181011.000000
  python3 met_track.py synthetic -o synthetic.dat --lon0 -86 --lat0 26 --heading 0 --speed 5 --vmax 120 --pc 940
  python3 met_track.py evaluate met-vortex.trk --model 2 --time 20181010.120000 --points pts.txt
"""
import argparse
import datetime as dt
import math
import sys

import numpy as np

KT = 0.514444
NM = 1852.0
OMEGA = 7.2921159e-5
R_EARTH = 6371000.0
EPOCH = dt.datetime(1970, 1, 1)
QUAD = [45.0, 135.0, 225.0, 315.0]          # NE SE SW NW, bearing clockwise from north
ISO = [34.0, 50.0, 64.0]


# ------------------------------------------------------------------ geometry and time
def stamp(t):
    return t.strftime("%Y%m%d.%H%M%S")


def epoch(t):
    return (t - EPOCH).total_seconds()


def from_stamp(s):
    s = "%.6f" % float(s)
    d, f = s.split(".")
    return dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]), int(f[0:2]), int(f[2:4]), int(f[4:6]))


def dist_bearing(lon1, lat1, lon2, lat2):
    """great-circle distance [m] and bearing [deg clockwise from north] from point 1 to point 2"""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    a = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    d = 2 * R_EARTH * math.asin(math.sqrt(min(1.0, a)))
    b = math.degrees(math.atan2(math.sin(dl) * math.cos(p2), math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)))
    return d, b % 360.0


def coriolis(lat):
    return 2 * OMEGA * math.sin(math.radians(lat))


# ------------------------------------------------------------------ profiles
def fang(r, rm):
    """inflow angle [deg] of ADCIRC"""
    if r < rm:
        return 10.0 * r / rm
    if r < 1.2 * rm:
        return 10.0 + 75.0 * (r / rm - 1.0)
    return 25.0


def holland80(r, vmaxg, rm, B, dp, rho, f):
    x = (rm / r) ** B
    af = abs(f)
    return math.sqrt(B * dp / rho * x * math.exp(-x) + (r * af / 2) ** 2) - r * af / 2


def gahm(r, vmaxg, rm, Bg, phi, f):
    af = abs(f)
    ro = vmaxg / (af * rm) if af > 0 else 1e30
    x = (rm / r) ** Bg
    return math.sqrt(vmaxg ** 2 * (1 + 1 / ro) * math.exp(phi * (1 - x)) * x + (r * af / 2) ** 2) - r * af / 2


def gahm_bphi(vmaxg, rm, dp, rho, f):
    """Bg, phi of GAHM by iteration (as ADCIRC calcRmaxesFull)"""
    af = abs(f)
    B = rho * math.e * vmaxg ** 2 / dp
    phi = 1.0
    for _ in range(400):
        phi_new = 1 + vmaxg * rm * af / (B * (vmaxg ** 2 + vmaxg * rm * af))
        B_new = (vmaxg ** 2 + vmaxg * rm * af) * rho * math.exp(phi_new) / (phi_new * dp)
        done = abs(B_new - B) < 1e-10 and abs(phi_new - phi) < 1e-12
        B, phi = B_new, phi_new
        if done:
            break
    return B, phi


def holland10(r, vms, rm, bs, xn, rn):
    g = (rm / r) ** bs * math.exp(1 - (rm / r) ** bs)
    if r <= rm or rn <= rm:
        x = 0.5
    else:
        x = 0.5 + (r - rm) * (xn - 0.5) / (rn - rm)
    return vms * g ** x


# ------------------------------------------------------------------ evaluation (reference of src/met_vortex.cpp)
def quad_params(rec, q, r):
    """Rm, Bg, phi of quadrant q at radius r: isotachs linear in r, outside the outermost the outermost"""
    iso = [(rec["R"][q][k], rec["Rm"][q][k], rec["Bg"][q][k], rec["phi"][q][k]) for k in range(3) if rec["R"][q][k] > 0]
    if not iso:
        return rec["rmw"], rec["Bg0"], rec["phi0"]
    iso.sort(key=lambda a: -a[0])                         # outermost first
    if r >= iso[0][0]:
        return iso[0][1:]
    if r <= iso[-1][0]:
        return iso[-1][1:]
    for a, b in zip(iso[:-1], iso[1:]):
        if b[0] <= r <= a[0]:
            w = (r - b[0]) / (a[0] - b[0])
            return tuple(w * a[i] + (1 - w) * b[i] for i in (1, 2, 3))


def azimuth_params(rec, theta, r):
    """inverse-square interpolation in azimuth between the two quadrants around theta (ADCIRC spInterp)"""
    a = (theta - 45.0) % 360.0
    if a >= 360.0:                                        # -1e-16 % 360 = 360.0
        a = 0.0
    i = min(int(a // 90.0), 3)                            # quadrant at or before theta
    j = (i + 1) % 4
    d = a - 90.0 * i                                      # 0..90
    pi_, pj = quad_params(rec, i, r), quad_params(rec, j, r)
    if d < 1.0:
        return pi_
    if d > 89.0:
        return pj
    wi, wj = 1.0 / d ** 2, 1.0 / (90.0 - d) ** 2
    return tuple((wi * a + wj * b) / (wi + wj) for a, b in zip(pi_, pj))


def evaluate(rec, lonc, latc, ut, vt, lon, lat, model, rho, wrf=0.9, one2ten=0.8928, inflow=1, trans=1):
    """wind (east, north) [m/s] and pressure [Pa] of record rec centred at (lonc, latc) with translation (ut, vt)"""
    r, theta = dist_bearing(lonc, latc, lon, lat)
    if r < 1000.0:
        return 0.0, 0.0, rec["pc"]
    f = coriolis(latc)
    s = 1.0 if f >= 0 else -1.0
    dp = rec["pn"] - rec["pc"]
    if model == 1:
        rm = rec["rmw"]
        V = holland80(r, rec["vmaxg"], rm, rec["B80"], dp, rho, f)
        p = rec["pc"] + dp * math.exp(-(rm / r) ** rec["B80"])
        gam, sp = V / rec["vmaxg"], wrf * V
    elif model == 2:
        rm, Bg, phi = azimuth_params(rec, theta, r)
        V = gahm(r, rec["vmaxg"], rm, Bg, phi, f)
        p = rec["pc"] + dp * math.exp(-phi * (rm / r) ** Bg)
        gam, sp = V / rec["vmaxg"], wrf * V
    else:
        rm = rec["rmw"]
        V = holland10(r, rec["vms"], rm, rec["bs"], rec["xn"], rec["rn"])
        p = rec["pc"] + dp * math.exp(-(rm / r) ** rec["bs"])
        gam, sp = V / rec["vms"], V
    beta = fang(r, rm) if inflow == 1 else 0.0
    ang = math.radians((90.0 - theta) + s * (90.0 + beta))
    ue, vn = sp * math.cos(ang), sp * math.sin(ang)
    if trans == 1:
        ue, vn = ue + gam * ut, vn + gam * vt
    elif trans == 2:
        ue, vn = ue + ut, vn + vt
    return one2ten * ue, one2ten * vn, p


# ------------------------------------------------------------------ best track readers
def read_atcf(path):
    recs = {}
    for line in open(path):
        c = [x.strip() for x in line.split(",")]
        if len(c) < 11 or not c[2].isdigit():
            continue
        t = dt.datetime.strptime(c[2], "%Y%m%d%H") + dt.timedelta(minutes=int(c[3]) if c[3].isdigit() else 0)
        lat = int(c[6][:-1]) / 10.0 * (1 if c[6][-1] == "N" else -1)
        lon = int(c[7][:-1]) / 10.0 * (1 if c[7][-1] == "E" else -1)
        r = recs.setdefault(t, dict(time=t, lat=lat, lon=lon, vmax=float(c[8]), pc=float(c[9]) if c[9] else 0.0,
                                    pn=0.0, rmw=0.0, iso={}))
        rad = float(c[11]) if len(c) > 11 and c[11] else 0.0
        if rad in (34.0, 50.0, 64.0) and len(c) > 16:
            q = [float(x) if x else 0.0 for x in c[13:17]]
            if c[12] in ("NEQ", "") and any(v > 0 for v in q):
                r["iso"][rad] = q
        if len(c) > 17 and c[17]:
            r["pn"] = float(c[17])
        if len(c) > 19 and c[19]:
            r["rmw"] = float(c[19])
    return [recs[t] for t in sorted(recs)]


def read_ibtracs(path):
    import csv
    out = []
    rd = csv.DictReader(open(path))
    for row in rd:
        try:
            t = dt.datetime.strptime(row["ISO_TIME"].strip(), "%Y-%m-%d %H:%M:%S")
        except (ValueError, KeyError):
            continue                                       # the units line under the header
        num = lambda k: float(row[k]) if row.get(k, "").strip() not in ("", " ") else 0.0
        iso = {}
        for k in (34, 50, 64):
            q = [num("USA_R%d_%s" % (k, s)) for s in ("NE", "SE", "SW", "NW")]
            if any(v > 0 for v in q):
                iso[float(k)] = q
        out.append(dict(time=t, lat=num("LAT"), lon=num("LON"), vmax=num("USA_WIND"), pc=num("USA_PRES"),
                        pn=num("USA_POCI"), rmw=num("USA_RMW"), iso=iso))
    return out


# ------------------------------------------------------------------ fit
def translation(track):
    """translation (east, north) [m/s] from the positions, central differences"""
    n = len(track)
    for i, r in enumerate(track):
        a, b = track[max(i - 1, 0)], track[min(i + 1, n - 1)]
        dtt = (b["time"] - a["time"]).total_seconds()
        if dtt <= 0:
            r["ut"] = r["vt"] = 0.0
            continue
        d, brg = dist_bearing(a["lon"], a["lat"], b["lon"], b["lat"])
        sp = d / dtt
        r["ut"], r["vt"] = sp * math.sin(math.radians(brg)), sp * math.cos(math.radians(brg))


def fit_record(r, args):
    rho, wrf = args.rho, args.wrf
    f = coriolis(r["lat"])
    af = abs(f)
    s = 1.0 if f >= 0 else -1.0
    pc = (r["pc"] if r["pc"] > 0 else 1000.0) * 100.0
    pn = (r["pn"] if r["pn"] > r["pc"] else args.pn) * 100.0
    dp = max(pn - pc, 100.0)
    vt = math.hypot(r["ut"], r["vt"])
    vm = r["vmax"] * KT
    vmaxg = max(vm - vt, 0.5 * vm) / wrf
    vms = max(vm - vt, 0.5 * vm)
    rmw = r["rmw"] * NM if r["rmw"] > 0 else args.rmw * 1000.0
    out = dict(t=epoch(r["time"]), lon=r["lon"], lat=r["lat"], ut=r["ut"], vt=r["vt"], pc=pc, pn=pn,
               vmaxg=vmaxg, rmw=rmw, vms=vms)
    out["B80"] = min(max(rho * math.e * vmaxg ** 2 / dp, 1.0), 2.5)
    out["bs"] = min(max(rho * math.e * vms ** 2 / dp, 1.0), 2.5)
    out["Bg0"], out["phi0"] = gahm_bphi(vmaxg, rmw, dp, rho, f)
    # Holland 2010 outer exponent from the mean 34 kt radius
    r34 = [x for x in r["iso"].get(34.0, []) if x > 0]
    out["xn"], out["rn"] = 0.5, 0.0
    if r34:
        rn = float(np.mean(r34)) * NM
        vn = 34.0 * KT * (1 - vt / max(vm, 1e-6))
        if rn > rmw and vn < vms:
            g = (rmw / rn) ** out["bs"] * math.exp(1 - (rmw / rn) ** out["bs"])
            out["xn"] = min(max(math.log(vn / vms) / math.log(g), 0.3), 1.5)
            out["rn"] = rn
    # GAHM per quadrant and isotach
    R = [[0.0] * 3 for _ in range(4)]
    Rm = [[0.0] * 3 for _ in range(4)]
    Bg = [[0.0] * 3 for _ in range(4)]
    ph = [[0.0] * 3 for _ in range(4)]
    worst = 0.0
    for k, v in enumerate(ISO):
        if v * KT >= vm or v not in r["iso"]:
            continue
        for q in range(4):
            Rr = r["iso"][v][q] * NM
            if Rr <= 0:
                continue
            rm = min(rmw, 0.9 * Rr)
            for _ in range(20):
                beta = fang(Rr, rm) if args.inflow == 1 else 0.0
                ang = math.radians((90.0 - QUAD[q]) + s * (90.0 + beta))
                ax = wrf * math.cos(ang) + (r["ut"] / vmaxg if args.trans == 1 else 0.0)
                ay = wrf * math.sin(ang) + (r["vt"] / vmaxg if args.trans == 1 else 0.0)
                vg = v * KT / math.hypot(ax, ay)              # gradient wind at the isotach

                def h(x):
                    b_, p_ = gahm_bphi(vmaxg, x, dp, rho, f)
                    return gahm(Rr, vmaxg, x, b_, p_, f) - vg
                lo, hi = 1000.0, Rr * 0.999
                if vg >= vmaxg or h(lo) * h(hi) > 0:
                    rm_new = Rr if vg >= vmaxg or h(hi) < 0 else lo
                else:
                    from scipy.optimize import brentq
                    rm_new = brentq(h, lo, hi, xtol=0.01)
                if abs(rm_new - rm) < 0.1:
                    rm = rm_new
                    break
                rm = rm_new
            b_, p_ = gahm_bphi(vmaxg, rm, dp, rho, f)
            R[q][k], Rm[q][k], Bg[q][k], ph[q][k] = Rr, rm, b_, p_
    out.update(R=R, Rm=Rm, Bg=Bg, phi=ph)
    # check: the fitted vortex at each isotach (quadrant centre), 1-min surface wind as in the best track
    for k, v in enumerate(ISO):
        for q in range(4):
            if R[q][k] > 0:
                lon, lat = destination(r["lon"], r["lat"], QUAD[q], R[q][k])
                ue, vn, _ = evaluate(out, r["lon"], r["lat"], r["ut"], r["vt"], lon, lat, 2, rho, wrf, 1.0, args.inflow, args.trans)
                worst = max(worst, abs(math.hypot(ue, vn) - v * KT) / KT)
    out["fit_err_kt"] = worst
    return out


def destination(lon, lat, brg, d):
    p1, l1, b, a = math.radians(lat), math.radians(lon), math.radians(brg), d / R_EARTH
    p2 = math.asin(math.sin(p1) * math.cos(a) + math.cos(p1) * math.sin(a) * math.cos(b))
    l2 = l1 + math.atan2(math.sin(b) * math.sin(a) * math.cos(p1), math.cos(a) - math.sin(p1) * math.sin(p2))
    return math.degrees(l2), math.degrees(p2)


def write_trk(path, recs, rho, comment):
    with open(path, "w") as o:
        o.write("REEF3D-MET-TRACK 1 %d %.4f\n" % (len(recs), rho))
        o.write("# %s\n" % comment)
        o.write("# t lon lat ut vt pc pn vmaxg rmw B80 vms bs xn rn Bg0 phi0 | NE SE SW NW x 34 50 64 kt: R Rm Bg phi\n")
        for r in recs:
            v = [r["lon"], r["lat"], r["ut"], r["vt"], r["pc"], r["pn"], r["vmaxg"], r["rmw"], r["B80"], r["vms"],
                 r["bs"], r["xn"], r["rn"], r["Bg0"], r["phi0"]]
            for q in range(4):
                for k in range(3):
                    v += [r["R"][q][k], r["Rm"][q][k], r["Bg"][q][k], r["phi"][q][k]]
            o.write(stamp(EPOCH + dt.timedelta(seconds=r["t"])) + " " + " ".join("%.10g" % x for x in v) + "\n")


def read_trk(path):
    lines = [l for l in open(path) if l.strip() and not l.startswith("#")]
    h = lines[0].split()
    rho = float(h[3])
    recs = []
    for l in lines[1:]:
        c = l.split()
        x = [float(v) for v in c[1:]]
        r = dict(t=epoch(from_stamp(c[0])), lon=x[0], lat=x[1], ut=x[2], vt=x[3], pc=x[4], pn=x[5], vmaxg=x[6],
                 rmw=x[7], B80=x[8], vms=x[9], bs=x[10], xn=x[11], rn=x[12], Bg0=x[13], phi0=x[14])
        a = np.array(x[15:]).reshape(4, 3, 4)
        r["R"], r["Rm"], r["Bg"], r["phi"] = a[:, :, 0].tolist(), a[:, :, 1].tolist(), a[:, :, 2].tolist(), a[:, :, 3].tolist()
        recs.append(r)
    return recs, rho


def at_time(recs, t, lon, lat, model, rho, **kw):
    """time interpolation as REEF3D: centre and translation linear in time, the fields of the two records blended"""
    ts = [r["t"] for r in recs]
    if t <= ts[0] or len(recs) == 1:
        a = b = recs[0]; w = 0.0
    elif t >= ts[-1]:
        a = b = recs[-1]; w = 0.0
    else:
        k = int(np.searchsorted(ts, t, side="right")) - 1
        a, b = recs[k], recs[k + 1]
        w = (t - a["t"]) / (b["t"] - a["t"])
    lonc, latc = a["lon"] + w * (b["lon"] - a["lon"]), a["lat"] + w * (b["lat"] - a["lat"])
    ut, vt = a["ut"] + w * (b["ut"] - a["ut"]), a["vt"] + w * (b["vt"] - a["vt"])
    ua, va, pa = evaluate(a, lonc, latc, ut, vt, lon, lat, model, rho, **kw)
    if w == 0.0:
        return ua, va, pa
    ub, vb, pb = evaluate(b, lonc, latc, ut, vt, lon, lat, model, rho, **kw)
    return ua + w * (ub - ua), va + w * (vb - va), pa + w * (pb - pa)


# ------------------------------------------------------------------ commands
def cmd_fit(args):
    track = read_ibtracs(args.input) if args.input.lower().endswith(".csv") else read_atcf(args.input)
    track = [r for r in track if r["vmax"] > 0 and (args.start is None or stamp(r["time"]) >= args.start)
             and (args.end is None or stamp(r["time"]) <= args.end)]
    if not track:
        sys.exit("no records")
    translation(track)
    recs = [fit_record(r, args) for r in track]
    write_trk(args.output, recs, args.rho, "fitted from %s (tools/met_track.py, rho_a %g, Wrf %g, inflow %d, translation %d)"
              % (args.input, args.rho, args.wrf, args.inflow, args.trans))
    for r, f in zip(track, recs):
        print("%s  %6.2f %7.2f  Vmax %5.1f kt  pc %6.1f hPa  RMW %5.1f km  isotachs %-12s  |V - isotach| max %.3f kt"
              % (stamp(r["time"]), r["lat"], r["lon"], r["vmax"], f["pc"] / 100, f["rmw"] / 1000,
                 ",".join("%d" % k for k in sorted(r["iso"])), f["fit_err_kt"]))
    print("%s: %d records" % (args.output, len(recs)))


def cmd_synthetic(args):
    """straight track, constant intensity after a spin-up of 12 h; ATCF b-deck"""
    lines = []
    t0 = from_stamp(args.start)
    lon, lat = args.lon0, args.lat0
    for n in range(int(args.hours / 6) + 1):
        t = t0 + dt.timedelta(hours=6 * n)
        lon_n, lat_n = destination(args.lon0, args.lat0, args.heading, args.speed * 6 * 3600 * n)
        rad = {34: args.r34, 50: args.r50, 64: args.r64}
        for k in (34, 50, 64):
            if k > args.vmax or rad[k] is None:
                continue
            q = rad[k]
            lines.append("AL, 99, %s, , BEST,   0, %4d%s, %5d%s, %3d, %4d, HU, %3d, NEQ, %4d, %4d, %4d, %4d, %4d, %4d, %3d"
                         % (t.strftime("%Y%m%d%H"), round(abs(lat_n) * 10), "N" if lat_n >= 0 else "S",
                            round(abs(lon_n) * 10), "E" if lon_n >= 0 else "W", args.vmax, args.pc, k,
                            q[0], q[1], q[2], q[3], args.pn, 200, args.rmw))
    open(args.output, "w").write("\n".join(lines) + "\n")
    print("%s: %d lines" % (args.output, len(lines)))


def cmd_evaluate(args):
    recs, rho = read_trk(args.input)
    t = epoch(from_stamp(args.time))
    for l in open(args.points):
        if not l.strip() or l.startswith("#"):
            continue
        lon, lat = map(float, l.split()[:2])
        u, v, p = at_time(recs, t, lon, lat, args.model, rho, wrf=args.wrf, one2ten=args.one2ten, inflow=args.inflow, trans=args.trans)
        print("%.10f %.10f %.12e %.12e %.12e" % (lon, lat, u, v, p))


def cmd_info(args):
    recs, rho = read_trk(args.input)
    for r in recs:
        print("%s  %7.2f %6.2f  Vg %5.1f m/s  pc %6.1f hPa  RMW %5.1f km  B %.2f  Bg0 %.2f phi0 %.3f  xn %.2f"
              % (stamp(EPOCH + dt.timedelta(seconds=r["t"])), r["lon"], r["lat"], r["vmaxg"], r["pc"] / 100,
                 r["rmw"] / 1000, r["B80"], r["Bg0"], r["phi0"], r["xn"]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def phys(p):
        p.add_argument("--wrf", type=float, default=0.9, help="gradient to 10 m wind reduction (W 216)")
        p.add_argument("--inflow", type=int, default=1, help="inflow angle of ADCIRC 1, none 0 (W 218)")
        p.add_argument("--trans", type=int, default=1, help="translation: 1 scaled by Vg/Vmax, 2 full, 0 none (W 219)")

    f = sub.add_parser("fit")
    f.add_argument("input")
    f.add_argument("-o", "--output", default="met-vortex.trk")
    f.add_argument("--rho", type=float, default=1.15, help="air density of the vortex [kg/m3]")
    f.add_argument("--pn", type=float, default=1013.0, help="outer pressure [hPa] where the track has none")
    f.add_argument("--rmw", type=float, default=40.0, help="radius of maximum winds [km] where the track has none")
    f.add_argument("--start")
    f.add_argument("--end")
    phys(f)
    f.set_defaults(func=cmd_fit)

    s = sub.add_parser("synthetic")
    s.add_argument("-o", "--output", default="synthetic.dat")
    s.add_argument("--start", default="20181010.000000")
    s.add_argument("--hours", type=float, default=48.0)
    s.add_argument("--lon0", type=float, default=-86.0)
    s.add_argument("--lat0", type=float, default=26.0)
    s.add_argument("--heading", type=float, default=0.0, help="deg clockwise from north")
    s.add_argument("--speed", type=float, default=5.0, help="m/s")
    s.add_argument("--vmax", type=int, default=110, help="kt")
    s.add_argument("--pc", type=int, default=950, help="hPa")
    s.add_argument("--pn", type=int, default=1012, help="hPa")
    s.add_argument("--rmw", type=int, default=20, help="nm")
    s.add_argument("--r34", type=int, nargs=4, default=[150, 130, 90, 110], metavar=("NE", "SE", "SW", "NW"))
    s.add_argument("--r50", type=int, nargs=4, default=[80, 70, 50, 60])
    s.add_argument("--r64", type=int, nargs=4, default=[45, 40, 30, 35])
    s.set_defaults(func=cmd_synthetic)

    e = sub.add_parser("evaluate")
    e.add_argument("input")
    e.add_argument("--model", type=int, default=2)
    e.add_argument("--time", required=True)
    e.add_argument("--points", required=True, help="file with lines: lon lat")
    e.add_argument("--one2ten", type=float, default=0.8928)
    phys(e)
    e.set_defaults(func=cmd_evaluate)

    i = sub.add_parser("info")
    i.add_argument("input")
    i.set_defaults(func=cmd_info)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
