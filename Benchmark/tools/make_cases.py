#!/usr/bin/env python3
# Architect: Hans Bihs
"""
Writes the benchmark case directories (cases/<name>/case.json, control.txt, ctrl.txt).

The case files in cases/ are the source of truth and can be edited by hand; this script documents
how they were set up (gauge lists, geometry) and can regenerate them:  ./tools/make_cases.py
"""

import json
import math
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASES = os.path.join(ROOT, "cases")
TUT = os.environ.get("REEF3D_TUTORIALS", os.path.join(ROOT, "..", "Tutorials"))


def write_case(name, control, ctrl, case, files=None):
    d = os.path.join(CASES, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "control.txt"), "w") as f:
        f.write(control.strip() + "\n")
    with open(os.path.join(d, "ctrl.txt"), "w") as f:
        f.write(ctrl.strip() + "\n")
    with open(os.path.join(d, "case.json"), "w") as f:
        json.dump(case, f, indent=1)
        f.write("\n")
    for src, dst in (files or {}).items():
        shutil.copy(src, os.path.join(d, dst))


def gauges(xs, y, fmt="P 51 %.4f %.4f"):
    return "\n".join(fmt % (x, y) for x in xs)


def frange(a, b, h):
    n = int(round((b - a) / h))
    return [a + i * h for i in range(n + 1)]


# references (short label, full citation)
REF = {
    "beji": ("Beji & Battjes (1993)",
             "Beji, S. & Battjes, J.A. (1993) Experimental investigation of wave propagation over a bar. Coastal "
             "Engineering 19, 151-162; Luth, Klopman & Kitou (1994) / Dingemans (1994). Measured series: secondary, "
             "digitised (Basilisk test 'bar', gauge-4..11), see refdata/waves/PROVENANCE.md"),
    "berkhoff": ("Berkhoff, Booy & Radder (1982)",
                 "Berkhoff, J.C.W., Booy, N. & Radder, A.C. (1982) Verification computations with linear wave "
                 "propagation models. Coastal Engineering 6, 255-279. Sections 2, 3, 5, 7: secondary, digitised "
                 "(Basilisk example 'shoal'), see refdata/waves/PROVENANCE.md"),
    "synolakis": ("Synolakis (1987)",
                  "Synolakis, C.E. (1987) The runup of solitary waves. J. Fluid Mech. 185, 523-545; run-up law as "
                  "quoted in Synolakis et al. (2008) PAGEOPH 165 (NOAA NCTR benchmark 'canonical bathymetry')"),
    "ting": ("Ting & Kirby (1995)",
             "Ting, F.C.K. & Kirby, J.T. (1995) Dynamics of surf-zone turbulence in a strong plunging breaker. "
             "Coastal Engineering 24, 177-204 (plunging case H = 0.128 m, T = 5 s, 1:35 slope); breaking point "
             "x_b = 7.795 m from Derakhti et al. (2016) Table 1, see refdata/waves/ting_kirby_breaking/setup.txt"),
    "fenton": ("Fenton (1985) 5th-order Stokes",
               "Fenton, J.D. (1985) A fifth-order Stokes theory for steady waves. J. Waterway, Port, Coastal and "
               "Ocean Eng. 111(2), 216-234 (the generating theory, B 92 5, recorded with P 50)"),
    "stokes2": ("2nd-order Stokes theory",
                "Stokes second-order wave theory (e.g. Dean & Dalrymple 1991), recorded with P 50"),
    "ritter": ("Ritter (1892), SWASHES",
               "Ritter, A. (1892) Die Fortpflanzung der Wasserwellen. Z. Verein Deutscher Ing. 36; exact solution "
               "as in Delestre et al. (2013) SWASHES, Int. J. Numer. Meth. Fluids 72, 269-300"),
    "stoker": ("Stoker (1957), SWASHES",
               "Stoker, J.J. (1957) Water Waves. Interscience; exact solution as in Delestre et al. (2013) SWASHES"),
    "solitary": ("Solitary wave celerity",
                 "Boussinesq solitary wave, c = sqrt(g (d + H)) (e.g. Synolakis 1987; Laitone 1960)"),
    "martin": ("Martin & Moyce (1952)",
               "Martin, J.C. & Moyce, W.J. (1952) An experimental study of the collapse of liquid columns on a "
               "rigid horizontal plane. Phil. Trans. R. Soc. A 244, 312-324; n^2 = 2, a = 2.25 in, surge front "
               "Z(T) digitised from Fig. 3 (secondary: PySPH db_exp_data.py), see refdata/dambreak/PROVENANCE.md"),
    "kleefsman": ("Kleefsman et al. (2005)",
                  "Kleefsman, K.M.T. et al. (2005) A volume-of-fluid based simulation method for wave impact "
                  "problems. J. Comput. Phys. 206, 363-393 (MARIN, SPHERIC Test 2). P1/P3 pressures: secondary, "
                  "digitised (PySPH), see refdata/dambreak/PROVENANCE.md"),
    "cylinder": ("Williamson St(Re); Qu et al. (2013)",
                 "Williamson (1988) St-Re fit as given by Beaudan & Moin (1994) Eq. 1 (50 < Re < 180); mean C_D at "
                 "Re = 100: Park et al. (1998) 1.33, Posdziech & Grundmann (2007) 1.325, Qu et al. (2013) 1.319 "
                 "(Qu et al. 2013 J. Fluids Struct. 39, Table 1), see refdata/misc/cylinder_re100_unconfined.txt"),
    "chen": ("Chen et al. (2014)",
             "Chen, L.F., Zang, J., Hillis, A.J., Morgan, G.C.J. & Plummer, A.R. (2014) Numerical investigation of "
             "wave-structure interaction using OpenFOAM. Ocean Engineering 88, 91-109; measured data from the "
             "REEF3D tutorial 11_14 (Bihs et al. 2016, Computers & Fluids 140)"),
    "irschik": ("Irschik et al. (2002)",
                "Irschik, K., Sparboom, U. & Oumeraci, H. (2002) Breaking wave characteristics for the loading of "
                "a slender pile. Proc. 28th ICCE; GWK Hannover; measured data from the REEF3D tutorial 11_15 "
                "(Kamath et al. 2016, Ocean Engineering)"),
    "hulme": ("Hulme (1982) linear theory",
              "Hulme, A. (1982) The wave forces acting on a floating hemisphere undergoing forced periodic "
              "oscillations. J. Fluid Mech. 121, 443-463 (added mass and damping, table via LHEEA "
              "Hulme_Heaving_Sphere), see refdata/misc/hulme_hemisphere_heave.txt; equilibrium: Archimedes"),
    "scour": ("HEC-18 / Sumer et al. (1992)",
              "Arneson et al. (2012) FHWA HEC-18 (CSU pier scour equation); Sumer, Fredsoe & Christiansen (1992) "
              "J. Waterway, Port, Coastal, Ocean Eng. 118, S/D = 1.3 (sigma 0.7) for live-bed current scour; see "
              "refdata/misc/pier_scour.txt"),
}


def ref(k):
    return {"short": REF[k][0], "citation": REF[k][1]}


# ==============================================================================================
# SFLOW
# ==============================================================================================

def sflow_dambreak(name, hr, key):
    """1D dam break, flat frictionless bed, hydrostatic SFLOW. hl = 0.5 m, dam at x0 = 20 m."""
    hl, x0 = 0.5, 20.0
    gx = [12.0, 16.0, 19.0, 21.0, 24.0, 28.0, 32.0, 36.0]
    control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 0.01
B 10 0.0 40.0 0.0 0.01 0.0 1.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 210 3
A 211 4
A 220 0
A 217 1
F 60 %(f60)g
F 72 0.0 %(x0)g -1.0 1.0 %(hl)g
N 41 4.0
N 47 0.25
M 10 1
P 12 100
W 22 -9.81
%(g)s
""" % {"f60": max(hr, 1.0e-4), "hl": hl, "x0": x0, "g": gauges(gx, 0.005)}
    case = {
        "description": "1D dam break on a flat frictionless bed, h_l = %g m, h_r = %g m, dam at x = %g m, hydrostatic "
                       "SFLOW (A 220 0), WENO reconstruction; depth time series at 8 gauges for 0 < t < 4 s against "
                       "the exact %s solution." % (hl, hr, x0, "Ritter (dry bed)" if hr == 0 else "Stoker (wet bed)"),
        "reference": ref(key),
        "tags": ["sflow", "analytical", "1d", "quick"],
        "time": 4.0,
        "check": {"type": "swe_dambreak", "x0": x0, "hl": hl, "hr": hr, "z_bed": -max(hr, 1.0e-4), "datum": 0.0,
                  "t_start": 0.0, "t_end": 4.0, "tol": 0.025 if hr == 0 else 0.002},
        "levels": {
            "nightly": {"control_set": {"B 1": "0.05", "B 10": "0.0 40.0 0.0 0.05 0.0 1.0"},
                        "ctrl_set": {"P 51": [("%.4f 0.025" % x) for x in gx]},
                        "check_set": {"tol": 0.03 if hr == 0 else 0.005}},
            "release": {"np": 2}},
        "baseline": {"nightly": 0.0212 if hr == 0 else 0.0026, "release": 0.0198 if hr == 0 else 0.0005,
                     "commit": "ff1bb68"},
    }
    write_case(name, control, ctrl, case)


def sflow_solitary():
    """Tutorial 8_3: solitary wave H/d = 0.1 on d = 0.5 m, non-hydrostatic SFLOW."""
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.01
B 10 0.0 40.0 0.0 0.01 0.0 1.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 210 3
A 211 4
A 216 2
A 220 2
A 223 0.5
A 240 1
A 241 1
A 243 0
A 246 0
B 90 1
B 92 9
B 91 0.05 4.0
B 96 4.0 8.0
B 98 2
B 99 1
F 60 0.5
N 41 18.0
N 47 0.2
M 10 1
P 12 100
P 51 8.0 0.005
P 51 16.0 0.005
P 51 24.0 0.005
W 22 -9.81
"""
    case = {
        "description": "Solitary wave H/d = 0.1 (d = 0.5 m) propagating over constant depth with the non-hydrostatic "
                       "SFLOW (A 220 2), set-up of tutorial 8_3 extended to 40 m: celerity between x = 8 and 24 m "
                       "against c = sqrt(g(d+H)) and crest height at x = 24 m.",
        "reference": ref("solitary"),
        "tags": ["sflow", "analytical", "1d", "quick"],
        "time": 18.0,
        "check": {"type": "solitary", "x1": 8.0, "x2": 24.0, "d": 0.5, "H": 0.05, "tol_c": 0.01, "tol_H": 0.03},
        "levels": {
            "nightly": {"control_set": {"B 1": "0.02", "B 10": "0.0 40.0 0.0 0.02 0.0 1.0"},
                        "ctrl_set": {"P 51": ["8.0 0.01", "16.0 0.01", "24.0 0.01"]},
                        "check_set": {"tol_c": 0.01, "tol_H": 0.03}},
            "release": {"np": 2}},
    }
    write_case("sflow_solitary_celerity", control, ctrl, case)


# Beji & Battjes: measured gauges WG4..WG11 in the Basilisk frame (bar toe at x = 6 m)
BB_GAUGES = [("WG4", 10.5), ("WG5", 12.5), ("WG6", 13.5), ("WG7", 14.5), ("WG8", 15.7), ("WG9", 17.3),
             ("WG10", 19.0), ("WG11", 21.0)]


def bb_signals(shift):
    return [{"name": n, "source": "gauge", "index": i + 1, "data": "waves/beji_battjes_bar/%s_eta_timeseries.txt" % n,
             "data_y_factor": 0.01} for i, (n, x) in enumerate(BB_GAUGES)]


def bb_check(tol_rms, tol_h):
    return {"type": "timeseries", "signals": bb_signals(0), "align": [n for n, _ in BB_GAUGES], "lag_range": [-1.2, 1.2],
            "lag_step": 0.005, "local_lag": 0.1, "tol_rms": tol_rms, "tol_height": tol_h,
            "note": "WG4..WG11 at x = 10.5..21 m from the wave board (bar toe at 6 m); the measured time origin "
                    "is arbitrary, one common lag is fitted on all gauges, plus at most +-0.1 s (T/20) per gauge for the gauge "
                    "position and digitising uncertainty"}


def sflow_beji():
    shift = 5.0   # tutorial 8_4 geometry: bar toe at x = 11 m
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.01
B 10 0.0 40.0 0.0 0.01 0.0 0.8

S 61 11.0 17.0 0.0 0.01 0.0 0.3
S 10 17.0 19.0 0.0 0.01 0.0 0.3
S 61 19.0 22.0 0.0 0.01 0.3 0.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 220 3
B 90 1
B 92 4
B 93 0.02 2.02
B 96 5.0 9.0
B 98 2
B 99 1
F 60 0.40
N 41 42.0
N 47 0.2
M 10 1
P 12 100
W 10 0.0
W 22 -9.81
%s
""" % gauges([x + shift for _, x in BB_GAUGES], 0.005)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar, non-hydrostatic "
                       "SFLOW with the improved-dispersion quadratic pressure (A 220 3), geometry of tutorial 8_4 (bar toe at x = 11 m, i.e. +5 m against the experiment). "
                       "Surface elevation at WG4..WG11 against the measured series (shape incl. the released higher "
                       "harmonics on the lee side).",
        "reference": ref("beji"),
        "tags": ["sflow", "experiment", "1d", "dispersion"],
        "time": 42.0,
        "check": bb_check(0.35, 0.2),
        "levels": {
            "nightly": {"control_set": {"B 1": "0.025", "B 10": "0.0 40.0 0.0 0.025 0.0 0.8",
                                        "S 61": ["11.0 17.0 0.0 0.025 0.0 0.3", "19.0 22.0 0.0 0.025 0.3 0.0"],
                                        "S 10": "17.0 19.0 0.0 0.025 0.0 0.3"},
                        "ctrl_set": {"P 51": ["%.4f 0.0125" % (x + shift) for _, x in BB_GAUGES]},
                        "check_set": {"tol_rms": 0.35, "tol_height": 0.2}},
            "release": {"np": 2}},
    }
    write_case("sflow_beji_battjes_bar", control, ctrl, case)


# ==============================================================================================
# NHFLOW
# ==============================================================================================

def nhflow_beji():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 1500 1 10
B 10 0.0 35.0 0.0 0.01 0.0 1.0

S 61 6.0 12.0 0.0 0.01 0.0 0.3
S 10 12.0 14.0 0.0 0.01 0.0 0.3
S 61 14.0 17.0 0.0 0.01 0.3 0.0

B 103 5
B 113 2.5
B 116 1.0

M 10 4
"""
    ctrl = """
A 10 5
B 90 1
B 92 4
B 93 0.02 2.02
B 96 3.73 8.73
B 98 2
B 99 2
F 60 0.4
N 41 42.0
N 47 0.5
M 10 4
P 12 100
W 22 -9.81
%s
""" % gauges([x for _, x in BB_GAUGES], 0.005)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar with NHFLOW, "
                       "geometry of tutorial NHFLOW/2 (bar toe at x = 6 m as in the experiment). Surface elevation "
                       "at WG4..WG11 against the measured series.",
        "reference": ref("beji"),
        "tags": ["nhflow", "experiment", "2dv", "dispersion"],
        "np": 4,
        "time": 42.0,
        "check": bb_check(0.35, 0.3),
        "levels": {
            "nightly": {"np": 2, "control_set": {"B 2": "1000 1 8", "M 10": "2"},
                        "check_set": {"tol_rms": 0.45, "tol_height": 0.3}},
            "release": {}},
    }
    write_case("nhflow_beji_battjes_bar", control, ctrl, case)


def nhflow_stokes5():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 800 1 10
B 10 0.0 200.0 0.0 0.05 0.0 1.0

B 103 5
B 113 2.5
B 116 1.0

M 10 4
M 20 2
"""
    ctrl = """
A 10 5
B 90 1
B 92 5
B 93 1.0 4.5
B 96 25.0 50.0
B 98 2
B 99 1
F 60 4.01
I 12 1
N 41 80.0
N 47 0.5
M 10 4
P 12 100
P 50 50.0 0.025
P 50 75.0 0.025
P 50 100.0 0.025
P 50 125.0 0.025
P 51 50.0 0.025
P 51 75.0 0.025
P 51 100.0 0.025
P 51 125.0 0.025
W 22 -9.81
"""
    case = {
        "description": "5th-order Stokes waves (Fenton), H = 1 m, T = 4.5 s, d = 4.01 m (kd ~ 1, H/L ~ 0.04) "
                       "propagating 100 m in NHFLOW (tutorial NHFLOW/1 extended). The simulated surface elevation "
                       "must follow the generating theory (P 50) at x = 50..125 m: no loss of height, no phase "
                       "drift, the crest/trough asymmetry kept.",
        "reference": ref("fenton"),
        "tags": ["nhflow", "analytical", "2dv", "quick"],
        "np": 4,
        "time": 80.0,
        "check": {"type": "theory_gauges", "t_start": 50.0, "t_end": 80.0, "tol_rms": 0.10, "tol_height": 0.03},
        "levels": {
            "nightly": {"np": 2, "time": 60.0, "control_set": {"B 2": "400 1 8", "M 10": "2"},
                        "check_set": {"t_start": 45.0, "t_end": 60.0, "tol_rms": 0.12, "tol_height": 0.04}},
            "release": {}},
    }
    write_case("nhflow_stokes5_propagation", control, ctrl, case)


def synolakis_geometry(d, Hd, dx, x_toe, x_end, z_top):
    """plane beach 1:19.85, toe at x_toe, still water depth d; returns (x_shore, S 61 line)"""
    cot = 19.85
    x_shore = x_toe + cot * d
    z_end = (x_end - x_toe) / cot
    return x_shore, "S 61 %.3f %.3f 0.0 %g 0.0 %.4f" % (x_toe, x_end, dx, z_end)


def nhflow_synolakis():
    d, Hd = 1.0, 0.0185
    x_toe, x_end = 50.0, 79.0
    x_shore, s61 = synolakis_geometry(d, Hd, 0.05, x_toe, x_end, 1.6)
    control = """
C 11 6
C 12 3
C 13 3
C 14 21
C 15 21
C 16 3

B 1 0.05
B 10 0.0 %(xe)g 0.0 0.05 0.0 1.6
%(s61)s

M 10 4
""" % {"xe": x_end, "s61": s61}
    ctrl = """
A 10 5
B 90 1
B 92 9
B 93 %(H)g 10.0
B 96 10.0 0.0
B 98 2
B 99 0
F 60 %(d)g
N 41 45.0
N 47 0.5
M 10 4
P 12 100
P 51 20.0 0.025
P 51 %(xt)g 0.025
P 52 0.025
P 55 0.05
P 133 0.025
P 134 0.025
W 22 -9.81
""" % {"H": Hd * d, "d": d, "xt": x_toe}
    case = {
        "description": "Non-breaking solitary wave H/d = 0.0185 (d = 1 m) running up a 1:19.85 plane beach "
                       "(Synolakis 1987 canonical case), NHFLOW with wetting and drying. Maximum run-up (P 134) "
                       "against the run-up law R/d = 2.831 sqrt(cot b)(H/d)^1.25 = 0.0862 (valid for non-breaking "
                       "waves; breaking starts at H/d = 0.045 on this slope).",
        "reference": ref("synolakis"),
        "tags": ["nhflow", "analytical", "2dv", "runup"],
        "np": 4,
        "time": 45.0,
        "check": {"type": "runup_law", "d": d, "H_over_d": Hd, "cot_beta": 19.85, "swl": d, "tol": 0.08},
        "extra_checks": [{"type": "profiles", "d": d, "H_over_d": Hd, "swl": d, "x_shore": x_shore, "x_sign": -1.0,
                          "profiles": [[t, "waves/synolakis_runup_funwave/synolakis_H0.0185_t%d.dat" % t]
                                       for t in (30, 40, 50, 60, 70)], "tol": 0.4,
                          # ff1bb68: rms/H 0.09-0.16 for t* = 30-60 and 0.32 (nightly) / 0.33 (release) in the
                          # run-down at t* = 70, i.e. independent of the grid
                          "note": "measured profiles (NTHMP BP4 via FUNWAVE-TVD), x/d positive offshore"}],
        "levels": {
            "nightly": {"np": 2, "control_set": {"B 1": "0.1", "B 10": "0.0 %g 0.0 0.1 0.0 1.6" % x_end,
                                                 "S 61": "%.3f %.3f 0.0 0.1 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85),
                                                 "M 10": "2"},
                        "ctrl_set": {"P 51": ["20.0 0.05", "%g 0.05" % x_toe], "P 52": "0.05", "P 133": "0.05",
                                     "P 134": "0.05"},
                        "check_set": {"tol": 0.12}},
            "release": {"control_set": {"B 1": "0.025", "B 10": "0.0 %g 0.0 0.025 0.0 1.6" % x_end,
                                        "S 61": "%.3f %.3f 0.0 0.025 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85)},
                        "ctrl_set": {"P 51": ["20.0 0.0125", "%g 0.0125" % x_toe], "P 52": "0.0125",
                                     "P 133": "0.0125", "P 134": "0.0125"}}},
    }
    write_case("nhflow_synolakis_runup_nonbreaking", control, ctrl, case)

    # breaking case H/d = 0.28: profiles
    Hd = 0.28
    ctrl2 = """
A 10 5
A 550 1
B 90 1
B 92 9
B 93 %(H)g 10.0
B 96 10.0 0.0
B 98 2
B 99 0
F 60 %(d)g
N 41 40.0
N 47 0.5
M 10 4
P 12 100
P 52 0.025
P 55 0.02
P 133 0.025
P 134 0.025
W 22 -9.81
""" % {"H": Hd * d, "d": d}
    profiles = [[t, "waves/synolakis_runup/profile_Hd0.28_t%d.txt" % t] for t in range(10, 70, 5)]
    case2 = {
        "description": "Breaking solitary wave H/d = 0.28 (d = 1 m) on the 1:19.85 beach (Synolakis 1987), NHFLOW "
                       "with breaking (A 550 1) and wetting/drying: free-surface profiles at t sqrt(g/d) = 10..65 "
                       "(shoaling, breaking, run-up and run-down) against the measured profiles. The time origin is "
                       "fitted on the first profile.",
        "reference": ref("synolakis"),
        "tags": ["nhflow", "experiment", "2dv", "runup", "breaking"],
        "np": 4,
        "time": 40.0,
        "check": {"type": "profiles", "d": d, "H_over_d": Hd, "swl": d, "x_shore": x_shore, "profiles": profiles,
                  "tol": 0.22},
        "levels": {
            "nightly": {"np": 2, "control_set": {"B 1": "0.1", "B 10": "0.0 %g 0.0 0.1 0.0 1.6" % x_end,
                                                 "S 61": "%.3f %.3f 0.0 0.1 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85),
                                                 "M 10": "2"},
                        "ctrl_set": {"P 52": "0.05", "P 133": "0.05", "P 134": "0.05"},
                        "check_set": {"tol": 0.2}},
            "release": {"control_set": {"B 1": "0.025", "B 10": "0.0 %g 0.0 0.025 0.0 1.6" % x_end,
                                        "S 61": "%.3f %.3f 0.0 0.025 0.0 %.4f" % (x_toe, x_end, (x_end - x_toe) / 19.85)},
                        "ctrl_set": {"P 52": "0.0125", "P 133": "0.0125", "P 134": "0.0125"}}},
    }
    write_case("nhflow_synolakis_breaking_profiles", control, ctrl2, case2)


BERKHOFF_SECTIONS = [
    {"name": "section 2", "x": 3.0, "sign": -1.0, "data": "waves/berkhoff_shoal/section2_H_over_H0.txt"},
    {"name": "section 3", "x": 5.0, "sign": -1.0, "data": "waves/berkhoff_shoal/section3_H_over_H0.txt"},
    {"name": "section 5", "x": 9.0, "sign": -1.0, "data": "waves/berkhoff_shoal/section5_H_over_H0.txt"},
    {"name": "section 7", "y": 0.0, "sign": 1.0, "data": "waves/berkhoff_shoal/section7_H_over_H0.txt"},
]


def berkhoff_gauges():
    """P 51 lines at the measured points; domain frame = shoal frame + (10, 10)"""
    out = []
    for sec in BERKHOFF_SECTIONS:
        rows = [l.split() for l in open(os.path.join(ROOT, "refdata", sec["data"])) if l.strip() and l[0] != "#"]
        for r in rows:
            s = float(r[0])
            x = sec["x"] if "x" in sec else sec["sign"] * s
            y = sec["y"] if "y" in sec else sec["sign"] * s
            out.append("P 51 %.4f %.4f" % (x + 10.0, y + 10.0))
    return "\n".join(out)


def berkhoff(model):
    """elliptic shoal, 3D: domain x' = -10..15, y' = -10..10 (shoal frame) -> 0..25, 0..20"""
    control = """
C 11 6
C 12 21
C 13 21
C 14 7
C 15 21
C 16 3

B 1 0.04
B 10 0.0 25.0 0.0 20.0 0.0 0.6
B 103 5
B 113 2.0
B 116 0.6

G 10 1
G 15 2
G 31 0

M 10 16
M 20 2
"""
    if model == "nhflow":
        mk = "A 10 5\nB 92 4\nN 47 0.5"
        sig = "B 2 625 500 10"
    else:
        mk = "A 10 3\nA 310 3\nA 311 4\nA 343 0\nB 92 4\nN 47 1.0"
        sig = "B 2 625 500 10"
    control = control.replace("B 1 0.04", sig)
    ctrl = """
%(mk)s
B 90 1
B 93 0.0464 1.0
B 96 1.5 3.0
B 98 2
B 99 1
F 60 0.45
N 41 30.0
M 10 16
P 12 100
W 22 -9.81
%(g)s
""" % {"mk": mk, "g": berkhoff8_gauges()}
    case = {
        "description": "Berkhoff et al. (1982) elliptic shoal on a 1:50 slope rotated by 20 deg, regular waves "
                       "T = 1 s, H = 4.64 cm (3D refraction, diffraction, nonlinear focusing behind the shoal) with "
                       "%s. Wave heights H/H0 along sections 1-5 (across the tank at x = 1, 3, 5, 7, 9 m behind the "
                       "shoal centre) and 6-8 (along the tank at y = -2, 0, 2 m) against the measurements. Bathymetry from geo.dat written "
                       "by tools/berkhoff_geo.py (Basilisk shoal.c formula, minimum depth 0.07 m)." % model.upper(),
        "reference": ref("berkhoff8"),
        "tags": [model, "experiment", "3d", "refraction"],
        "np": 16,
        "time": 30.0,
        "generate": {"script": "berkhoff_geo.py", "args": [0.04]},
        "check": {"type": "section_heights", "sections": BERKHOFF8, "dx": 10.0, "dy": 10.0, "H0": 0.0464,
                  "t_start": 22.0, "t_end": 30.0, "tol": 0.15,
                  "note": "sections 1-8 from FUNWAVE-TVD (amplitudes, H/H0 = a/23.2 mm); y in the FUNWAVE frame, "
                          "which matches the Basilisk shoal.c frame used for the bathymetry (checked against the "
                          "Basilisk copy of section 5: rms difference 0.03)"},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 2": "250 200 8", "M 10": "4"},
                        "ctrl_set": {"M 10": "4"}, "generate": {"script": "berkhoff_geo.py", "args": [0.1]},
                        "check_set": {"tol": 0.3}},
            "release": {}},
    }
    if model == "fnpf":
        # nightly run (ff1bb68, 2 cores): rms 0.07-0.27 on sections 1-6 and 8, 0.31 on the centreline
        # section 7 (focus peak H/H0 1.85 vs 2.02, slightly downstream)
        case["levels"]["nightly"]["check_set"]["tol"] = 0.35
        # release run (dx 0.04): rms 0.09-0.18 on sections 1-6 and 8, 0.21 on section 7 (focus H/H0 2.17 vs 2.02)
        case["levels"]["release"]["check_set"] = {"tol": 0.25}
    if model == "nhflow":
        case["levels"]["nightly"]["control_set"]["B 2"] = "250 200 5"
        case["levels"]["nightly"]["xfail"] = (
            "hans_dev ff1bb68, dx = 0.1 m (L/15): the wave heights are about 60 % of the measured ones on all "
            "sections, already at the first point of section 7 (H/H0 0.67 vs 1.07); 5 or 8 sigma layers give "
            "the same (rms 0.34-0.72 on sections 1-7, 0.20 on section 8). Not a 3D effect: a 2D flume with the same "
            "wave at dx = 0.1 m loses 16 % of the height over 14 m with the default WENO-JS reconstruction "
            "(dx 0.05: 2 %, dx 0.025: none). WENO-Z (A 527 1) keeps 0.95 along the flume and reduces the "
            "Berkhoff error to 0.09-0.49 (overall 0.31); the rest is the coarse grid over the shoal (L/10)")
    write_case("%s_berkhoff_shoal" % model, control, ctrl, case)


# ==============================================================================================
# FNPF
# ==============================================================================================

def fnpf_beji():
    shift = 5.0   # tutorial 9_2 geometry: bar toe at x = 11 m
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 2000 1 10
B 10 0.0 40.0 0.0 0.02 0.0 1.0

S 61 11.0 17.0 0.0 0.02 0.0 0.3
S 10 17.0 19.0 0.0 0.02 0.0 0.3
S 61 19.0 22.0 0.0 0.02 0.3 0.0

B 103 5
B 113 3.0
B 116 1.0

M 10 2
M 20 2
"""
    ctrl = """
A 10 3
A 310 3
A 311 4
A 343 0
B 90 1
B 92 4
B 93 0.02 2.02
B 96 5.0 10.0
B 98 2
B 99 1
F 60 0.4
N 41 42.0
N 47 1.0
M 10 2
P 12 100
W 22 -9.81
%s
""" % gauges([x + shift for _, x in BB_GAUGES], 0.01)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar with FNPF, "
                       "geometry of tutorial 9_2 (bar toe at x = 11 m, i.e. +5 m against the experiment). Surface "
                       "elevation at WG4..WG11 against the measured series.",
        "reference": ref("beji"),
        "tags": ["fnpf", "experiment", "2dv", "dispersion", "quick"],
        "np": 2,
        "time": 42.0,
        "check": bb_check(0.3, 0.25),
        "levels": {
            "nightly": {"np": 1, "control_set": {"B 2": "800 1 8", "M 10": "1", "B 10": "0.0 40.0 0.0 0.05 0.0 1.0",
                                                 "S 61": ["11.0 17.0 0.0 0.05 0.0 0.3", "19.0 22.0 0.0 0.05 0.3 0.0"],
                                                 "S 10": "17.0 19.0 0.0 0.05 0.0 0.3"},
                        "ctrl_set": {"M 10": "1", "P 51": ["%.4f 0.025" % (x + shift) for _, x in BB_GAUGES]},
                        "check_set": {"tol_rms": 0.3, "tol_height": 0.2}},
            "release": {}},
    }
    write_case("fnpf_beji_battjes_bar", control, ctrl, case)


def fnpf_stokes2():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 800 1 10
B 10 0.0 40.0 0.0 0.05 0.0 2.0
B 103 5
B 113 2.5
B 116 2.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 3
A 310 3
A 311 4
A 343 0
B 90 1
B 92 4
B 91 0.04 4.0
B 96 4.0 8.0
B 98 2
B 99 1
F 60 2.0
N 41 50.0
N 47 1.0
M 10 1
P 12 100
P 50 15.0 0.025
P 50 20.0 0.025
P 50 25.0 0.025
P 51 15.0 0.025
P 51 20.0 0.025
P 51 25.0 0.025
W 22 -9.81
"""
    case = {
        "description": "2nd-order Stokes waves H = 0.04 m, L = 4 m in deep water (d = 2 m), FNPF, tutorial 9_1: "
                       "the simulated surface elevation must follow the generating theory (P 50) at x = 15..25 m.",
        "reference": ref("stokes2"),
        "tags": ["fnpf", "analytical", "2dv", "quick"],
        "check": {"type": "theory_gauges", "t_start": 35.0, "t_end": 50.0, "tol_rms": 0.15, "tol_height": 0.02,
                  "note": "hans_dev ff1bb68: the simulated wave is about 0.4 % faster than the theory (phase lead "
                          "0.017 s at x = 15 m, 0.032 s at 25 m, at both levels), height within 1 %"},
        "time": 50.0,
        "levels": {
            "nightly": {"time": 40.0, "control_set": {"B 2": "400 1 8"},
                        "check_set": {"t_start": 30.0, "t_end": 40.0, "tol_rms": 0.17, "tol_height": 0.03}},
            "release": {}},
    }
    write_case("fnpf_stokes2_propagation", control, ctrl, case)


def fnpf_ting_kirby():
    # tutorial 9_3: slope from x = 5.8 m (depth 0.4 m), depth 0.38 m at x = 6.5 m (Ting & Kirby x = 0)
    x0 = 5.8 + 0.02 * 35.0
    xs = frange(8.0, 18.0, 0.1)
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 600 1 10
B 10 0.0 30.0 0.0 0.025 0.0 0.748
B 103 5
B 113 1.0
B 116 0.4

S 61 5.8 32.0 0.0 0.025 0.0 0.748

M 10 2
M 20 2
"""
    ctrl = """
A 10 3
A 341 2.0
A 343 1
A 350 1
A 351 3
A 352 3
A 365 0.0025
B 90 1
B 92 8
B 93 0.128 5.0
B 96 9.5 9.0
B 98 3
B 99 1
F 60 0.4
N 41 100.0
N 47 1.0
M 10 2
P 12 100
W 22 -9.81
%s
""" % gauges(xs, 0.0125)
    case = {
        "description": "Plunging breaker of Ting & Kirby (1995): cnoidal waves H = 0.128 m, T = 5 s on a 1:35 slope, "
                       "FNPF with the viscosity breaking model (tutorial 9_3 at its resolution dx = 0.05 m, Dirichlet generation). Breaking point "
                       "= location of the maximum wave height from 101 gauges between x = 8 and 18 m, against the "
                       "measured breaking point x_b = 7.795 m from the x = 0 station (d = 0.38 m, here x = %.2f m)."
                       % x0,
        "reference": ref("ting"),
        "tags": ["fnpf", "experiment", "2dv", "breaking"],
        "np": 2,
        "time": 100.0,
        "check": {"type": "breaking_point", "x_b": 7.795, "x_offset": x0, "t_start": 70.0, "t_end": 100.0,
                  "tol": 0.4},
        "levels": {
            "nightly": {"np": 1, "time": 70.0, "control_set": {"M 10": "1"},
                        "ctrl_set": {"M 10": "1"}, "check_set": {"t_start": 45.0, "t_end": 70.0, "tol": 0.8}},
            "release": {}},
        "xfail": "on hans_dev ff1bb68 (dx = 0.05 m, tutorial resolution) the maximum wave height is reached at "
                 "x = 12.4 m (d = 0.21 m, H = 0.20 m), 1.9 m before the measured breaking point (d = 0.157 m); "
                 "with dx = 0.025 m the run stops at t = 91 s (N 61 velocity limit) after the wave train becomes "
                 "irregular at t ~ 75 s",
    }
    write_case("fnpf_ting_kirby_plunging", control, ctrl, case)


# ==============================================================================================
# CFD
# ==============================================================================================

def cfd_martin_moyce():
    a = 0.05715  # 2.25 in
    Zs = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    L = 14 * a
    control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 %(dx)g
B 10 0.0 %(L).5f 0.0 %(dx)g 0.0 %(H).5f

M 10 2
M 20 2
""" % {"dx": a / 40, "L": L, "H": 3 * a}
    ctrl = """
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 54 %(a).5f
F 56 %(h).5f
N 40 3
N 41 0.5
N 47 0.2
M 10 2
P 12 100
T 10 0
W 22 -9.81
%(g)s
""" % {"a": a, "h": 2 * a, "g": gauges([Z * a for Z in Zs], a / 80)}
    case = {
        "description": "2D collapse of a water column a x 2a (n^2 = 2, a = 2.25 in = 0.05715 m) in a closed tank, "
                       "two-phase CFD with level set (set-up of tutorial 11_1 scaled to the experiment). Surge-front "
                       "arrival at gauges Z = x/a = 1.5..10 as T = t sqrt(2g/a) against Martin & Moyce (1952). "
                       "The front has arrived when the surface is 0.05 a above the bed.",
        "reference": ref("martin"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "quick"],
        "np": 2,
        "time": 0.5,
        "check": {"type": "front_arrival", "a": a, "x_wall": 0.0, "z_bed": 0.0, "h_thr": 0.05 * a, "datum": 0.0,
                  "data": "dambreak/martin_moyce_1952_surge_n2_2_a2.25in_pysph.dat", "tol_max": 0.25,
                  "tol_mean": 0.08,
                  "note": "the largest deviation is at Z = 1.5 (front 20 % early): the experiment includes the "
                          "opening of the column, the simulation starts from rest without it"},
        "levels": {
            "nightly": {"np": 1, "control_set": {"B 1": "%g" % (a / 20), "B 10": "0.0 %.5f 0.0 %g 0.0 %.5f" % (L, a / 20, 3 * a),
                                                 "M 10": "1"},
                        "ctrl_set": {"M 10": "1", "P 51": ["%.5f %g" % (Z * a, a / 40) for Z in Zs]},
                        "check_set": {"tol_max": 0.25, "tol_mean": 0.1}},
            "release": {}},
    }
    write_case("cfd_dambreak_2d_martin_moyce", control, ctrl, case)


def cfd_kleefsman():
    H = 0.55
    tf = 1.0 / math.sqrt(9.81 / H)
    pf = 1000.0 * 9.81 * H
    control = """
C 11 21
C 12 21
C 13 21
C 14 21
C 15 21
C 16 3

B 1 0.01
B 10 0.0 3.22 0.0 1.0 0.0 1.0
S 10 0.67 0.83 0.3 0.7 0.0 0.16

M 10 16
M 20 2
"""
    ctrl = """
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 51 2.0
F 56 0.55
N 40 3
N 41 4.0
N 47 0.2
M 10 16
P 12 100
T 10 0
W 1 1000.0
W 22 -9.81
P 64 0.835 0.5 0.021
P 64 0.835 0.5 0.101
P 51 2.725 0.5
P 51 1.0 0.5
"""
    sig = lambda n, f: {"name": "P%d" % n, "source": "pressure", "n": 1 if n == 1 else 2,
                        "data": "dambreak/kleefsman_2005_P%d_pysph.dat" % n,
                        "data_t_factor": tf, "data_y_factor": pf}
    case = {
        "description": "MARIN 3D dam break with a box obstacle (Kleefsman et al. 2005, SPHERIC Test 2): tank "
                       "3.22 x 1 x 1 m, water column 1.22 x 1 x 0.55 m behind the door at x = 2 m, box 0.16 x 0.40 "
                       "x 0.16 m at x = 0.67..0.83 m. Pressure at P1 (z = 0.021 m) and P3 (z = 0.101 m) on the "
                       "water-facing box face against the measurements: impact peak value and time, rms error. "
                       "Sensor heights follow PySPH/LS-DYNA (0.021/0.101 m; the ComFLOW page gives 0.025/0.099 m).",
        "reference": ref("kleefsman"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "impact"],
        "np": 16,
        "time": 4.0,
        "check": {"type": "timeseries", "signals": [sig(1, 1), sig(3, 2)], "align": [], "lag": 0.0,
                  "tol_rms": 0.35, "tol_height": None, "rms_about_mean": False,
                  "peak": {"signals": ["P1"], "window": [0.3, 0.8], "tol_value": 0.30, "tol_time": 0.03}},
        "levels": {
            "nightly": {"np": 4, "time": 2.0, "control_set": {"B 1": "0.025", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 64": ["0.85 0.5 0.021", "0.85 0.5 0.101"]},
                        "check_set": {"tol_rms": 0.5, "window": [0.0, 2.0],
                                      "peak": {"signals": ["P1"], "window": [0.3, 0.8], "tol_value": 0.3,
                                               "tol_time": 0.1}}},
            "release": {}},
    }
    write_case("cfd_dambreak_3d_kleefsman", control, ctrl, case)


def cfd_cylinder():
    """laminar 2D cylinder Re = 100, D = 0.1 m, U = 1 m/s, nu = 1e-3; domain 30D x 16D"""
    D, U, nu = 0.1, 1.0, 1.0e-3
    Hz = 1.6
    control = """
C 11 1
C 12 3
C 13 3
C 14 2
C 15 3
C 16 3

B 1 0.002
B 10 0.0 3.0 0.0 0.002 0.0 1.6
B 101 11
B 127 0.002 0.03 0.9 0.4 1.05
B 103 11
B 129 0.002 0.03 0.79 0.3 1.05

S 32 0.8 0.79 0.05

M 10 4
M 20 2
"""
    ctrl = """
B 10 0
B 20 2
B 60 1
B 61 1
D 10 4
D 20 2
D 22 2
D 30 1
F 30 0
F 40 0
I 11 1
N 40 3
N 41 50.0
N 47 0.3
M 10 4
P 12 100
P 81 0.7 0.9 -1.0 1.0 0.69 0.89
T 10 0
W 1 1000.0
W 2 %(nu)g
W 10 %(Q)g
"""
    case = {
        "description": "Laminar flow past a circular cylinder at Re = U D / nu = 100 (D = 0.1 m, U = 1 m/s, nu = "
                       "1e-3 m2/s), single-phase CFD in 2D, domain 30 D x 16 D (8 D upstream, blockage 6 %, cylinder 0.1 D below "
                       "the centre line so that the shedding starts), "
                       "cell-size stretching around the cylinder (D/50 release, D/25 nightly), implicit diffusion with the wall ghost values at the solid (D 22 2: with the default D 22 1 the cylinder wall is free of shear in the implicit diffusion, C_D drops by a third and the wake does not shed). Strouhal number from "
                       "the lift and mean drag coefficient after the shedding has settled.",
        "reference": ref("cylinder"),
        "tags": ["cfd", "analytical", "2d", "single-phase", "quick"],
        "np": 4,
        "time": 50.0,
        "check": {"type": "vortex_shedding", "U": U, "D": D, "W": 0.002, "nu": nu, "rho": 1000.0, "lift": "Fz",
                  "Cd_ref": 1.33, "t_start": 30.0, "t_end": 50.0, "tol_St": 0.05, "tol_Cd": 0.10,
                  "note": "blockage D/H = 1/16 raises St and C_D by a few per cent against the unconfined values"},
        "levels": {
            "nightly": {"np": 2, "time": 40.0,
                        "control_set": {"B 1": "0.004", "B 10": "0.0 3.0 0.0 0.004 0.0 1.6",
                                        "B 127": "0.004 0.04 0.9 0.4 1.08", "B 129": "0.004 0.04 0.79 0.3 1.08",
                                        "M 10": "2"},
                        "ctrl_set": {"M 10": "2", "W 10": "%g" % (U * Hz * 0.004)},
                        "check": {"type": "vortex_shedding", "U": U, "D": D, "W": 0.004, "nu": nu, "rho": 1000.0,
                                  "lift": "Fz", "Cd_ref": 1.33, "t_start": 25.0, "t_end": 40.0, "tol_St": 0.08,
                                  "tol_Cd": 0.15},
                        # ff1bb68: C_D 1.01 from the sign error in the P 81 wall shear (patch 0002,
                        # in hans_dev since d80a861): C_D 1.33, St 0.1745, C_L 0.31 on this grid
                        },
            "release": {}},
    }
    ctrl = ctrl % {"nu": nu, "Q": U * Hz * 0.002}
    write_case("cfd_cylinder_re100", control, ctrl, case)


def cfd_beji():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.01
B 10 0.0 30.0 0.0 0.01 0.0 0.6

S 61 6.0 12.0 0.0 0.01 0.0 0.3
S 10 12.0 14.0 0.0 0.01 0.0 0.3
S 61 14.0 17.0 0.0 0.01 0.3 0.0

M 10 8
M 20 2
"""
    ctrl = """
B 10 1
B 50 0.0001
B 90 1
B 92 2
B 93 0.02 2.02
B 96 3.73 5.0
B 98 2
B 99 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 0.6
F 60 0.4
I 12 1
N 40 3
N 41 42.0
N 47 0.2
M 10 8
P 12 100
W 22 -9.81
%s
""" % gauges([x for _, x in BB_GAUGES], 0.005)
    case = {
        "description": "Regular waves (H = 2 cm, T = 2.02 s) over the Beji & Battjes submerged bar, two-phase CFD in "
                       "2D (tutorial 11_10 with the flume extended to 30 m so that all gauges lie outside the "
                       "numerical beach). Surface elevation at WG4..WG11 against the measured series.",
        "reference": ref("beji"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "dispersion"],
        "np": 8,
        "time": 42.0,
        "check": bb_check(0.35, 0.2),
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.02", "B 10": "0.0 30.0 0.0 0.02 0.0 0.6",
                                                 "S 61": ["6.0 12.0 0.0 0.02 0.0 0.3", "14.0 17.0 0.0 0.02 0.3 0.0"],
                                                 "S 10": "12.0 14.0 0.0 0.02 0.0 0.3", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 51": ["%.4f 0.01" % x for _, x in BB_GAUGES]},
                        "check_set": {"tol_rms": 0.4, "tol_height": 0.2}},
            "release": {}},
    }
    write_case("cfd_beji_battjes_bar", control, ctrl, case)


def cfd_ting_kirby():
    # tutorial 11_11: slope from x = 13.8 m, depth 0.38 m at x = 14.5 m
    x0 = 13.8 + 0.02 * 35.0
    xs = frange(18.0, 25.0, 0.1)
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.005
B 10 0.0 40.0 0.0 0.005 0.0 1.00
S 61 13.8 40.0 0.0 0.005 0.0 0.748

M 10 32
M 20 2
"""
    ctrl = """
B 10 1
B 11 1
B 50 0.0001
B 90 1
B 92 8
B 93 0.128 5.0
B 96 9.8 0.0
B 98 2
B 99 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 1.0
F 60 0.4
I 12 1
N 40 3
N 41 60.0
N 47 0.1
M 10 32
P 12 100
T 10 2
T 36 2
W 10 0.0
W 22 -9.81
%s
""" % gauges(xs, 0.0025)
    case = {
        "description": "Plunging breaker of Ting & Kirby (1995), cnoidal waves H = 0.128 m, T = 5 s on a 1:35 "
                       "slope, two-phase CFD with k-omega (tutorial 11_11). Breaking point = location of the "
                       "maximum wave height from gauges every 0.1 m between x = 18 and 25 m, against the measured "
                       "x_b = 7.795 m from the d = 0.38 m station (here x = %.2f m)." % x0,
        "reference": ref("ting"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "breaking"],
        "np": 32,
        "time": 60.0,
        "check": {"type": "breaking_point", "x_b": 7.795, "x_offset": x0, "t_start": 40.0, "t_end": 60.0,
                  "tol": 0.4},
        "levels": {
            "nightly": {"np": 4, "time": 45.0,
                        "control_set": {"B 1": "0.02", "B 10": "0.0 40.0 0.0 0.02 0.0 1.00",
                                        "S 61": "13.8 40.0 0.0 0.02 0.0 0.748", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 51": ["%.4f 0.01" % x for x in xs]},
                        "check_set": {"t_start": 30.0, "t_end": 45.0, "tol": 1.0}},
            "release": {}},
    }
    write_case("cfd_ting_kirby_plunging", control, ctrl, case)


def cfd_chen():
    t = os.path.join(TUT, "REEF3D_CFD", "11_14 Non-Breaking Wave Forces")
    control = """
C 11 6
C 12 21
C 13 3
C 14 7
C 15 21
C 16 3

B 1 0.025
B 10 0.0 18.0 1.5 3.0 0.0 1.0

S 33 7.50 1.50 0.125

M 10 16
M 20 2
"""
    ctrl = """
B 10 1
B 11 1
B 50 0.0001
B 90 1
B 92 4
B 93 0.07 1.22
B 96 2.11 4.22
B 98 2
B 99 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 1.0
F 60 0.505
I 12 1
N 40 3
N 41 33.0
N 47 0.1
M 10 16
P 12 100
T 10 2
T 36 2
W 22 -9.81
P 51 5.00 1.5125
P 81 7.325 7.675 1.4 1.675 0.0 1.0
"""
    case = {
        "description": "Regular waves (2nd-order Stokes, H = 0.07 m, T = 1.22 s, d = 0.505 m) on a vertical "
                       "cylinder D = 0.25 m (Chen et al. 2014), two-phase CFD with k-omega, tutorial 11_14 on the "
                       "half domain y = 1.5..3 m (symmetry plane through the cylinder axis, force x 2). Surface "
                       "elevation at x = 5 m and inline force against the measurements; the time lag is fitted on "
                       "the gauge.",
        "reference": ref("chen"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "wave-force"],
        "np": 16,
        "time": 33.0,
        "check": {"type": "timeseries",
                  "signals": [{"name": "eta x=5 m", "source": "gauge", "index": 1,
                               "data": "wave_forces/chen2014_eta_x5.txt"},
                              {"name": "Fx", "source": "force", "n": 1, "component": "Fx", "factor": 2.0,
                               "data": "wave_forces/chen2014_force.txt"}],
                  "align": ["eta x=5 m"], "lag_range": [-1.3, 1.3], "lag_step": 0.005,
                  "tol_rms": 0.25, "tol_height": 0.10},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.05", "M 10": "4"}, "ctrl_set": {"M 10": "4",
                        "P 51": "5.0 1.525"},
                        "check_set": {"tol_rms": 0.4, "tol_height": 0.2}},
            "release": {}},
    }
    write_case("cfd_wave_force_chen2014", control, ctrl, case)


def cfd_irschik():
    control = """
C 11 6
C 12 21
C 13 3
C 14 8
C 15 21
C 16 3

B 1 0.05
B 10 0.0 54.0 2.5 5.0 0.0 7.0

S 33 44.0 2.50 0.35
S 61 21.0 44.0 2.5 5.0 0.0 2.30
S 10 44.0 54.0 2.5 5.0 0.0 2.30

M 10 64
M 20 2
"""
    ctrl = """
B 10 1
B 50 0.0001
B 90 1
B 92 5
B 93 1.3 4.0
B 96 21.0 0.0
B 98 2
B 99 3
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 42 7.0
F 60 3.80
I 12 1
N 40 3
N 41 27.0
N 47 0.1
M 10 64
P 12 100
T 10 2
T 36 2
W 22 -9.81
P 51 43.65 2.55
P 81 43.30 44.70 2.5 3.20 2.30 8.0
"""
    case = {
        "description": "Breaking wave impact on a vertical pile D = 0.7 m in the GWK (Irschik et al. 2002): 5th-order "
                       "Stokes waves H = 1.3 m, T = 4 s on d = 3.8 m, 1:10 slope, the wave breaks at the pile "
                       "(tutorial 11_15, half domain with a symmetry plane, force x 2). Release level only: the "
                       "slamming force (peak and rms over the measured impact) against the measured force; the "
                       "lag is fitted on the force itself because the position of the measured gauge in the tutorial "
                       "data is not documented.",
        "reference": ref("irschik"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "wave-force", "breaking"],
        "np": 64,
        "time": 27.0,
        "check": {"type": "timeseries",
                  "signals": [{"name": "Fx", "source": "force", "n": 1, "component": "Fx", "factor": 2.0,
                               "data": "wave_forces/irschik2002_force.txt"}],
                  "align": ["Fx"], "lag_range": [-4.0, 4.0], "lag_step": 0.005,
                  "tol_rms": 0.4, "tol_height": 0.25},
        "levels": {"release": {}},
    }
    write_case("cfd_breaking_wave_force_irschik2002", control, ctrl, case)


def hulme_values():
    import sys
    sys.path.insert(0, ROOT)
    import benchmark
    T, zeta, kR = benchmark.hulme_prediction(benchmark.read_refdata("misc/hulme_hemisphere_heave.txt"), 1.0, 0.5)
    return T, zeta


def cfd_sphere():
    t = os.path.join(TUT, "REEF3D_CFD", "11_16 Heave Decay of Sphere")
    control = """
C 11 21
C 12 21
C 13 21
C 14 21
C 15 21
C 16 3

B 1 0.025
B 10 0.0 6.0 0.0 6.0 0.0 4.0

B 101 11
B 127 0.025 0.2 3.0 1.5 1.08
B 102 11
B 128 0.025 0.2 3.0 1.5 1.08
B 103 11
B 129 0.025 0.2 2.0 1.5 1.08

M 10 32
M 20 2
"""
    ctrl = """
B 10 1
B 50 0.00001
B 90 1
B 99 2
B 107 0.0 0.0 0.0 6.0 0.5
B 107 6.0 6.0 0.0 6.0 0.5
B 107 0.0 6.0 0.0 0.0 0.5
B 107 0.0 6.0 6.0 6.0 0.5
D 10 4
D 20 2
D 30 1
F 30 3
F 35 5
F 40 3
F 60 2.0
I 10 1
N 40 4
N 41 8.0
N 47 0.3
M 10 32
P 12 100
W 1 1000.0
W 22 -9.81
X 10 1
X 11 0 0 1 0 0 0
X 21 500.0
X 180 1
X 182 3.0 3.0 2.1
"""
    case = {
        "description": "Free heave decay of a sphere R = 0.5 m of half the water density (tutorial 11_16, heave "
                       "only, initial offset 0.1 m = 0.2 R, damping zones at the walls). Equilibrium from Archimedes "
                       "(centre at the still water level), damped period and damping ratio of the first periods "
                       "against linear potential theory with Hulme's (1982) hemisphere coefficients "
                       "(T = %.3f s, zeta = %.3f). The offset is moderate, so a few per cent nonlinear deviation is "
                       "expected and allowed in the tolerances." % hulme_values(),
        "reference": ref("hulme"),
        "tags": ["cfd", "semi-analytical", "3d", "two-phase", "6dof"],
        "np": 32,
        "time": 8.0,
        "files": ["floating.stl"],
        "check": {"type": "heave_decay", "R": 0.5, "mass_ratio": 1.0, "z_eq": 2.0, "n_periods": 2,
                  "hulme": "misc/hulme_hemisphere_heave.txt", "tol_T": 0.05, "tol_zeta": 0.25, "tol_eq": 0.02,
                  "eq_window": 1.0},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.05", "B 127": "0.05 0.3 3.0 1.5 1.1",
                                                 "B 128": "0.05 0.3 3.0 1.5 1.1", "B 129": "0.05 0.3 2.0 1.5 1.1",
                                                 "M 10": "4"},
                        "ctrl_set": {"M 10": "4"},
                        "check_set": {"tol_T": 0.06, "tol_zeta": 0.25, "tol_eq": 0.02}},
            "release": {}},
    }
    write_case("cfd_sphere_heave_decay", control, ctrl, case, {os.path.join(t, "floating.stl"): "floating.stl"})


def cfd_pier_scour():
    D, h, U = 0.2, 0.3, 0.3
    Fr = U / math.sqrt(9.81 * h)
    hec = 2.0 * 1.0 * 1.0 * 1.1 * (h / D) ** 0.35 * Fr ** 0.43
    control = """
C 11 1
C 12 21
C 13 21
C 14 2
C 15 21
C 16 3

B 1 0.02
B 10 1.0 5.5 0.0 2.0 -0.3 0.5

S 33 2.5 1.0 0.1

M 10 32
M 13 0
"""
    ctrl = """
B 10 1
B 11 1
B 50 0.006
B 60 1
B 61 2
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 60 0.3
I 10 1
N 40 2
N 47 0.3
M 10 32
P 12 100
P 122 1
T 10 2
T 36 1
W 10 0.18
W 22 -9.81
S 10 1
S 11 3
S 12 0
S 13 5.0
S 14 0.3
S 16 4
S 19 75600
S 20 0.00097
S 21 5.0
S 22 2570.0
S 24 0.75
S 30 0.035
S 41 1
S 42 1
S 43 100
S 44 1
S 50 1
S 57 0.0
S 60 0
S 73 0.00 0.5 0.0 1.0 0.0
S 80 4
S 81 35.0
S 82 10.0
S 90 1
S 91 10
"""
    case = {
        "description": "Local scour at a vertical circular pier D = 0.2 m in a steady current (h = 0.3 m, U = 0.3 "
                       "m/s, d50 = 0.97 mm), CFD with k-omega and the Exner model (tutorial 11_18). Release: "
                       "21 h of sediment time (S 19 75600); the equilibrium scour depth S/D must lie in the range "
                       "of established estimates: HEC-18/CSU %.2f (K1 = 1, K2 = 1, K3 = 1.1), Sumer et al. 1.3 "
                       "(sigma 0.7). Nightly: 1 h of sediment time on a coarse grid as a smoke test (scour must "
                       "start, S/D between 0.05 and 2)." % hec,
        "reference": ref("scour"),
        "tags": ["cfd", "empirical", "3d", "sediment"],
        "np": 32,
        "check": {"type": "scour_depth", "D": D, "band": [0.9, 2.0], "hec18": round(hec, 3), "bed0": 0.0},
        "levels": {
            "nightly": {"np": 4, "control_set": {"B 1": "0.04", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "S 19": "3600", "S 17": "5.0"},
                        "check_set": {"band": [0.05, 2.0]}},
            "release": {}},
    }
    write_case("cfd_pier_scour", control, ctrl, case)


# ==============================================================================================
# second round: data from public model test suites (FUNWAVE-TVD, NHWAVE, AQUAgpusph, olaFlow)
# ==============================================================================================

REF.update({
    "spheric2": ("Kleefsman et al. (2005), SPHERIC Test 2",
                 "Kleefsman, K.M.T. et al. (2005) A volume-of-fluid based simulation method for wave impact "
                 "problems. J. Comput. Phys. 206, 363-393 (MARIN, SPHERIC Test 2). Measured P1-P8 and H1-H4 "
                 "(instrument records, dt 1 ms) from the AQUAgpusph repository, see refdata/PROVENANCE_git_sources.md"),
    "berkhoff8": ("Berkhoff, Booy & Radder (1982)",
                  "Berkhoff, J.C.W., Booy, N. & Radder, A.C. (1982) Verification computations with linear wave "
                  "propagation models. Coastal Engineering 6, 255-279. Sections 1-8 (amplitudes) from the "
                  "FUNWAVE-TVD benchmark car_berkhoff_2d, see refdata/PROVENANCE_git_sources.md"),
    "briggs": ("Briggs et al. (1995) conical island",
               "Briggs, M.J., Synolakis, C.E., Harkins, G.S. & Green, D.R. (1995) Laboratory experiments of tsunami "
               "runup on a circular island. Pure Appl. Geophys. 144, 569-593 (NTHMP benchmark BP6). Gauges 6, 9, 16, "
               "22 (FUNWAVE-TVD), run-up around the island (NHWAVE), see refdata/PROVENANCE_git_sources.md"),
    "mase": ("Mase & Kirby (1992)",
             "Mase, H. & Kirby, J.T. (1992) Hybrid frequency-domain KdV equation for random wave transformation. "
             "Proc. 23rd ICCE, 474-487; measured records (dt 0.05 s, 12 gauges on a 1:20 slope) from the "
             "FUNWAVE-TVD benchmark car_mase_kirby, see refdata/PROVENANCE_git_sources.md"),
    "lin1998": ("Lin (1998), Liu et al. (1999)",
                "Lin, P. (1998) Numerical modeling of breaking waves. PhD thesis, Cornell University; Liu, P.L.-F., "
                "Lin, P., Chang, K.-A. & Sakakiyama, T. (1999) Numerical modeling of wave interaction with porous "
                "structures. J. Waterway, Port, Coastal, Ocean Eng. 125, 322-330. Digitised profiles from the "
                "olaFlow tutorial CR35_dambreak, see refdata/PROVENANCE_git_sources.md"),
    "thacker": ("Thacker (1981), SWASHES",
                "Thacker, W.C. (1981) Some exact solutions to the nonlinear shallow-water wave equations. J. Fluid "
                "Mech. 107, 499-508; planar surface in a parabola as in Delestre et al. (2013) SWASHES 4.2.1"),
    "jonswap": ("JONSWAP spectrum (DNV-RP-C205)",
                "Hasselmann et al. (1973) JONSWAP spectrum in the form of DNV-RP-C205 (2010) Sec. 3.5.5, "
                "A_gamma = 1 - 0.287 ln(gamma)"),
    "sloshing": ("Linear sloshing theory",
                 "Free oscillation of the first sloshing mode, linear potential theory T = 2 pi / sqrt(g k tanh kh), "
                 "k = pi/L (e.g. Faltinsen & Timokha 2009, Sloshing, Ch. 4)"),
})


def cfd_kleefsman2():
    """MARIN dam break, SPHERIC Test 2 geometry (AQUAgpusph frame shifted by +1.992 m in x)"""
    H = 0.55
    xf = 0.8245            # water-facing box face
    def signals(tp_front, tp_top, th, th_behind):
        """per-signal rms tolerances: box front P1-P4, box top P5-P8, water heights, height behind the box"""
        sig = []
        for n in range(1, 9):
            sig.append({"name": "P%d" % n, "source": "pressure", "n": n,
                        "data": "dambreak/spheric_test2/spheric_test2_P1-P8_H1-H4.dat", "data_col": n,
                        "tol_rms": tp_front if n <= 4 else tp_top})
        for k, (nm, x) in enumerate([("H x=1.456", 1.456), ("H x=0.960", 0.960), ("H x=0.464", 0.464),
                                      ("H x=2.606", 2.606)]):
            sig.append({"name": nm, "source": "gauge", "index": k + 1, "datum": 0.0, "clip_min": 0.0,
                        "data": "dambreak/spheric_test2/spheric_test2_P1-P8_H1-H4.dat", "data_col": 9 + k,
                        "tol_rms": th_behind if x < 0.8 else th})
        return sig
    # The top probes P5-P8 see a thin film / splash the grid does not resolve (rms 0.83-0.95 at the nightly
    # grid up to t = 1.7 s) and the flow behind the box arrives late (H x=0.464 rms 0.46): looser tolerances
    # there. Release values are provisional (not run yet).
    sig = signals(0.6, 0.9, 0.35, 0.5)
    sig_n = signals(0.7, 1.0, 0.35, 0.6)
    control = """
C 11 21
C 12 21
C 13 21
C 14 21
C 15 21
C 16 3

B 1 0.01
B 10 0.0 3.22 0.0 1.0 0.0 1.0
S 10 0.6635 0.8245 0.2985 0.7015 0.0 0.161

M 10 16
M 20 2
"""
    probes = ["%.4f 0.5 %.3f" % (xf + 0.005, 0.021 + 0.04 * i) for i in range(4)] + \
             ["%.4f 0.5 0.166" % (xf - 0.021 - 0.04 * i) for i in range(4)]
    probes_n = ["0.8375 0.5 %.3f" % (0.021 + 0.04 * i) for i in range(4)] + \
               ["%.4f 0.5 0.1875" % (xf - 0.021 - 0.04 * i) for i in range(4)]
    ctrl = """
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 46 3
F 47 10
F 51 1.992
F 56 0.55
N 40 3
N 41 6.0
N 47 0.2
M 10 16
P 12 100
T 10 0
W 1 1000.0
W 22 -9.81
%s
P 51 1.456 0.5
P 51 0.960 0.5
P 51 0.464 0.5
P 51 2.606 0.5
""" % "\n".join("P 64 " + p for p in probes)
    case = {
        "description": "MARIN 3D dam break with a box obstacle (Kleefsman et al. 2005, SPHERIC Test 2): tank 3.22 x 1 x "
                       "1 m, water column 1.228 x 1 x 0.55 m behind the door at x = 1.992 m, box 0.161 x 0.403 x 0.161 "
                       "m with its water-facing side at x = 0.8245 m (SPHERIC/AQUAgpusph dimensions). Pressures P1-P4 "
                       "on the box front (z = 0.021-0.141 m) and P5-P8 on its top, water heights at x = 1.456, 0.960, "
                       "0.464 (behind the box) and 2.606 m (reservoir) against the measured records: rms errors, peak "
                       "value and time of P1. The probes sit 5 mm off the box surface (one cell centre at the nightly "
                       "grid).",
        "reference": ref("spheric2"),
        "tags": ["cfd", "experiment", "3d", "two-phase", "impact"],
        "np": 16,
        "time": 6.0,
        "check": {"type": "timeseries", "signals": sig, "align": [], "lag": 0.0, "tol_rms": 0.6, "tol_height": None,
                  "rms_about_mean": False,
                  "note_mass": "F 46 3 (level-set volume correction, works with N 40 3 since hans_dev d80a861): "
                               "without it the volume drops by 18 % and the run stops at t = 1.77 s",
                  "peak": {"signals": ["P1"], "window": [0.3, 0.8], "tol_value": 0.3, "tol_time": 0.03},
                  "note": "the water-height columns follow AQUAgpusph's reading (column 10 = gauge nearest the "
                          "reservoir); the labels H1-H4 of the data file are therefore not used"},
        "levels": {
            "nightly": {"np": 4, "time": 2.0, "control_set": {"B 1": "0.025", "M 10": "4"},
                        "ctrl_set": {"M 10": "4", "P 64": probes_n},
                        "check_set": {"tol_rms": 0.7, "window": [0.0, 2.0], "signals": sig_n,
                                      "peak": {"signals": ["P1"], "window": [0.3, 0.8], "tol_value": 0.3,
                                               "tol_time": 0.1}}},
            "release": {}},
    }
    write_case("cfd_dambreak_3d_kleefsman", control, ctrl, case)


BERKHOFF8 = ([{"name": "section %d" % n, "x": x, "sign": 1.0, "col": 1, "y_factor": 1.0 / 23.2,
               "data": "waves/berkhoff_shoal_funwave/section%d.dat" % n}
              for n, x in [(1, 1.0), (2, 3.0), (3, 5.0), (4, 7.0), (5, 9.0)]] +
             [{"name": "section %d" % n, "y": y, "sign": -1.0, "col": k, "y_factor": 1.0 / 23.2,
               "data": "waves/berkhoff_shoal_funwave/section678.dat"}
              for n, y, k in [(6, -2.0, 1), (7, 0.0, 2), (8, 2.0, 3)]])


def berkhoff8_gauges():
    out = []
    for sec in BERKHOFF8:
        rows = [l.split() for l in open(os.path.join(ROOT, "refdata", sec["data"])) if l.strip() and l[0] != "#"]
        for r in rows:
            s = float(r[0])
            x = sec["x"] if "x" in sec else sec["sign"] * s
            y = sec["y"] if "y" in sec else sec["sign"] * s
            out.append("P 51 %.4f %.4f" % (x + 10.0, y + 10.0))
    return "\n".join(sorted(set(out), key=out.index))


def conical_geo_case(model, case, H_over_d, dx_rel, dx_night):
    """Briggs et al. (1995): conical island, solitary wave, flat depth 0.32 m"""
    d, xc, yc = 0.32, 15.0, 13.8
    radii = [round(1.40 + 0.03 * i, 3) for i in range(36)]
    angles = [0.0, 45.0, 67.5, 90.0, 112.5, 135.0, 180.0, 225.0, 247.5, 270.0, 292.5, 315.0]
    run = {"type": "runup_angles", "data": "waves/conical_island/briggs_case%s_runup.dat" % case,
           "data_y_factor": 0.01, "angles": angles, "radii": radii, "xc": xc, "yc": yc, "r_toe": 3.6,
           "slope": 0.25, "height": 0.625, "depth": d, "theta_front": 270.0, "thr": 0.003,
           "tol": 0.08, "tol_max": 0.12,
           "note": "run-up angle convention taken from the data: 270 deg faces the incident wave (largest run-up), "
                   "90 deg is the lee side; not stated in the source"}
    pts = []
    for th in angles:
        phi = math.radians(th - 270.0 + 180.0)
        for r in radii:
            pts.append((xc + r * math.cos(phi), yc + r * math.sin(phi)))
    # wave gauges 6, 9, 16, 22: distances from the FUNWAVE grid indices (cone centre i = 461, j = 277)
    g4 = [("gauge 6", xc - 3.60, yc), ("gauge 9", xc - 2.45, yc), ("gauge 16", xc, yc + 2.55), ("gauge 22", xc + 2.55, yc)]
    sig = [{"name": n, "source": "gauge", "index": k + 1, "data": "waves/conical_island/briggs_case%s_gauges_6_9_16_22.dat" % case,
            "data_col": k + 1} for k, (n, _, _) in enumerate(g4)]
    gauges = "\n".join("P 51 %.4f %.4f" % (x, y) for _, x, y in g4) + "\n" + \
             "\n".join("P 51 %.4f %.4f" % p for p in pts)
    control = """
C 11 6
C 12 21
C 13 21
C 14 7
C 15 21
C 16 3

B 1 %(dx)g
B 10 0.0 26.0 0.0 27.6 0.0 1.0

G 10 1
G 15 2
G 31 0

M 10 8
M 20 2
""" % {"dx": dx_rel}
    if model == "sflow":
        mk = "A 10 2\nA 246 1\nN 47 0.3"
    else:
        mk = "A 10 5\nA 550 1\nN 47 0.5"
        control = control.replace("B 1 %g" % dx_rel, "B 2 %d %d 5" % (round(26 / dx_rel), round(27.6 / dx_rel)))
    ctrl = """
%(mk)s
B 90 1
B 92 9
B 93 %(H)g 10.0
B 96 2.0 3.0
B 98 2
B 99 1
F 60 0.32
N 41 25.0
M 10 8
P 12 100
W 22 -9.81
%(g)s
""" % {"mk": mk, "H": H_over_d * d, "g": gauges}
    lab_peak = {"A": 31.0, "B": 29.8, "C": 28.76}[case]   # measured peak at gauge 6 (lab time)
    nightly_cs = {"B 1": "%g" % dx_night, "M 10": "4"} if model == "sflow" else \
        {"B 2": "%d %d 5" % (round(26 / dx_night), round(27.6 / dx_night)), "M 10": "4"}
    case_d = {
        "description": "Solitary wave H/d = %g (case %s, d = 0.32 m) on the conical island of Briggs et al. (1995) "
                       "with %s: island of 7.2 m toe diameter, 1:4 slope, 0.625 m high, centred at (15, 13.8) m; "
                       "bathymetry from geo.dat (tools/conical_geo.py). Surface elevation at gauges 6 (toe, front), "
                       "9 (front slope), 16 (side) and 22 (lee) against the measured records (common lag fitted), and "
                       "the maximum run-up on 12 radial lines (gauges every 0.03 m along the slope) against the "
                       "measured run-up. Gauge radii from the FUNWAVE-TVD grid indices (3.60, 2.45, 2.55, 2.55 m)."
                       % (H_over_d, case, model.upper()),
        "reference": ref("briggs"),
        "tags": [model, "experiment", "3d", "runup"],
        "np": 8,
        "time": 25.0,
        "generate": {"script": "conical_geo.py", "args": [min(dx_rel, 0.05)]},
        "check": {"type": "timeseries", "signals": sig, "align": [n for n, _, _ in g4], "lag_range": [-40.0, 0.0],
                  "lag_step": 0.01, "local_lag": 0.2, "window": [lab_peak - 4.0, lab_peak + 12.0],
                  "tol_rms": 0.5, "tol_height": 0.25},
        "extra_checks": [run],
        "levels": {
            "nightly": {"np": 4, "control_set": nightly_cs, "ctrl_set": {"M 10": "4"},
                        "generate": {"script": "conical_geo.py", "args": [0.05]},
                        "check_set": {"tol_rms": 0.6, "tol_height": 0.4},
                        "extra_checks": [dict(run, tol=0.1, tol_max=0.15)]},
            "release": {}},
    }
    write_case("%s_conical_island_%s" % (model, case), control, ctrl, case_d)


MK_DEPTHS = [470, 350, 300, 250, 200, 175, 150, 125, 100, 75, 50, 25]


def mase_kirby(model):
    """irregular waves on a 1:20 slope; linear wave components from the measured record at the slope
    toe (h = 47 cm), generated in a relaxation zone; wave origin (B 105) at the position of that gauge"""
    x_g, x_toe, x_end = 2.75, 3.0, 14.5
    gx = [(h, (x_g if h == 470 else x_toe + (0.47 - h / 1000.0) / 0.05)) for h in MK_DEPTHS]
    control = """
C 11 6
C 12 3
C 13 3
C 14 21
C 15 21
C 16 3

B 1 0.02
B 10 0.0 %(xe)g 0.0 0.02 0.0 0.7
S 61 %(xt)g %(xe)g 0.0 0.02 0.0 %(ze).3f

M 10 2
M 20 2
""" % {"xt": x_toe, "xe": x_end, "ze": (x_end - x_toe) * 0.05}
    if model == "nhflow":
        mk = "A 10 5\nA 550 1\nB 89 1\nN 47 0.5"
        control = control.replace("B 1 0.02", "B 2 725 1 8")
    elif model == "fnpf":
        mk = "A 10 3\nA 341 2.0\nA 343 1\nA 350 1\nA 351 3\nA 352 3\nA 365 0.0025\nB 89 1\nN 47 1.0"
        control = control.replace("B 1 0.02", "B 2 725 1 10")
    else:
        mk = "A 10 2\nA 220 3\nA 246 1\nN 47 0.3"
    ctrl = """
%(mk)s
B 90 1
B 92 51
B 96 2.5 0.0
B 98 2
B 99 0
B 101 2
B 102 5.0
F 60 0.47
N 41 420.0
M 10 2
P 12 1000
W 22 -9.81
%(g)s
""" % {"mk": mk, "g": "\n".join("P 51 %.3f 0.01" % x for _, x in gx)}
    gauges = [{"name": "h = %g cm" % (h / 10.0), "index": k + 1,
               "data": "waves/mase_kirby/mase_kirby_eta_h%03dmm.dat" % h} for k, (h, _) in enumerate(gx)]
    night_cs = {"M 10": "1"}
    if model == "nhflow":
        night_cs["B 2"] = "362 1 5"
    elif model == "fnpf":
        night_cs["B 2"] = "362 1 8"
    else:
        night_cs.update({"B 1": "0.04", "B 10": "0.0 %g 0.0 0.04 0.0 0.7" % x_end,
                         "S 61": "%g %g 0.0 0.04 0.0 %.3f" % (x_toe, x_end, (x_end - x_toe) * 0.05)})
    case = {
        "description": "Irregular waves breaking on a 1:20 slope (Mase & Kirby 1992), %s. The waves are linear "
                       "components (0.2-3 Hz, 1147 components; SFLOW 0.2-1.6 Hz with the improved-dispersion pressure "
                       "A 220 3, since the depth-averaged model cannot carry the short components) from the FFT of the first 409.6 s of the measured record at the slope "
                       "toe (h = 47 cm; waverecon.dat written by tools/mase_kirby_waverecon.py, the decomposition "
                       "of FUNWAVE-TVD's fft4wavemaker.m), generated in a 2.5 m relaxation zone (B 92 51; the record "
                       "is reproduced at the wave origin x = 0, the toe gauge sits at x = 2.75 m; B 105 cannot move "
                       "the origin because it also moves the relaxation zone). Significant wave height Hm0 = 4 sigma and skewness at the "
                       "gauges (h = 47 ... 2.5 cm) against the measured records over the same time window." % model.upper(),
        "reference": ref("mase"),
        "tags": [model, "experiment", "2dv", "irregular", "breaking"],
        "np": 2,
        "time": 420.0,
        "generate": {"script": "mase_kirby_waverecon.py", "args": [0.2, 1.6] if model == "sflow" else []},
        "check": {"type": "wave_stats", "gauges": gauges, "t_start": 30.0, "t_end": 409.0, "data_dt": 0.05,
                  "data_y_factor": 0.01, "tol_Hm0": 0.25, "tol_skew": 0.3,
                  "note": "the measured record contains frequencies below 0.2 Hz that the generated waves do not "
                          "(Hm0 at the toe gauge: 6.62 cm generated, 6.65 cm measured)"},
        "levels": {
            "nightly": {"np": 1, "time": 150.0 if model == "sflow" else 200.0, "ctrl_set": {"M 10": "1"},
                        "control_set": night_cs,
                        "check_set": {"t_start": 30.0, "t_end": 150.0 if model == "sflow" else 200.0,
                                      "tol_Hm0": 0.3, "tol_skew": 0.35}},
            "release": {}},
    }
    if model == "fnpf":
        case["levels"]["release"]["xfail"] = (
            "hans_dev ff1bb68, release grid (725 x 10): Hm0 too high in the inner surf zone, +18 % at h = 5 cm "
            "and +39 % at h = 2.5 cm (nightly grid: +9 % / +26 %); within 9 % for h >= 7.5 cm and skewness "
            "within 0.13 everywhere. Not the breaking model: in the wind-wave band (f > 0.3 Hz) FNPF is "
            "within 6 % at all gauges; the excess is infragravity energy (f < 0.3 Hz, about 2x measured at "
            "the shoreline, 1.6x already at the toe): free long waves from the linear (first-order) wave "
            "generation, reflected at the static coastline")
    write_case("%s_mase_kirby_irregular" % model, control, ctrl, case)


def cfd_porous_dambreak():
    """Lin (1998) dam break through a crushed-rock dam (VRANS), set-up of the olaFlow tutorial CR35_dambreak"""
    times = ["0.35", "0.75", "1.15", "1.55", "1.95"]
    control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 0.005
B 10 0.0 0.89 0.0 0.005 0.0 0.58

M 10 2
M 20 2
"""
    ctrl = """
B 270 0.30 0.59 -1.0 1.0 0.0 0.58 0.49 0.0159 500.0 2.0
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
F 70 0.0 0.29 -1.0 1.0 0.0 0.35
F 70 0.0 0.89 -1.0 1.0 0.0 0.025
N 40 3
N 41 2.0
N 47 0.2
M 10 2
P 12 100
P 52 0.0025
P 55 0.05
T 10 0
W 22 -9.81
"""
    case = {
        "description": "Dam break through a porous dam of crushed rock (Lin 1998; Liu et al. 1999): reservoir 0.29 m "
                       "long and 0.35 m deep, 25 mm water downstream, porous dam at x = 0.30-0.59 m (n = 0.49, d50 = "
                       "15.9 mm), CFD with VRANS (B 270, alpha = 500, beta = 2.0, added mass 0.34), set-up of the "
                       "olaFlow tutorial CR35_dambreak. Free-surface profiles at t = 0.35 ... 1.95 s against the "
                       "measured profiles.",
        "reference": ref("lin1998"),
        "tags": ["cfd", "experiment", "2dv", "two-phase", "porous"],
        "np": 2,
        "time": 2.0,
        "check": {"type": "profiles_abs", "h_ref": 0.35, "t_tol": 0.03, "tol": 0.12,
                  "profiles": [[float(t), "porous/lin1998/lin1998_CR35_t%ss.dat" % t] for t in times],
                  "note": "VRANS coefficients alpha = 500, beta = 2.0 as in the REEF3D regression case; olaFlow's "
                          "tutorial uses a = 50, b = 2.0 in its own formulation"},
        "levels": {
            "nightly": {"np": 1, "control_set": {"B 1": "0.01", "B 10": "0.0 0.89 0.0 0.01 0.0 0.58", "M 10": "1"},
                        "ctrl_set": {"M 10": "1", "P 52": "0.005"}, "check_set": {"tol": 0.15}},
            "release": {}},
    }
    write_case("cfd_porous_dambreak_lin1998", control, ctrl, case)


def sflow_thacker():
    control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 0.005
B 10 0.0 4.0 0.0 0.005 0.0 2.5

G 10 1
G 15 2
G 31 0

M 10 1
M 20 2
"""
    ctrl = """
A 10 2
A 210 3
A 211 4
A 217 1
A 220 0
A 251 0.5
F 60 1.375
N 41 10.0303
N 47 0.25
M 10 1
P 12 100
W 22 -9.81
%s
""" % gauges([0.6, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.25, 3.4], 0.0025)
    case = {
        "description": "Thacker's planar surface oscillating in a parabolic bowl without friction (SWASHES 4.2.1: "
                       "a = 1 m, h0 = 0.5 m, L = 4 m, 5 periods of 2.006 s), hydrostatic SFLOW with wetting and "
                       "drying; bed from geo.dat (tools/thacker_geo.py), initial planar surface with the slope "
                       "option A 251. Water depth at 9 gauges (including points that fall dry and get wet again) "
                       "against the exact solution.",
        "reference": ref("thacker"),
        "tags": ["sflow", "analytical", "1d", "wetting-drying", "quick"],
        "time": 10.0303,
        "generate": {"script": "thacker_geo.py", "args": [0.005]},
        "check": {"type": "thacker1d", "a": 1.0, "h0": 0.5, "L": 4.0, "swl": 1.375, "t_start": 0.0,
                  "t_end": 10.0303, "h_dry": 0.002, "tol": 0.05},
        "levels": {
            "nightly": {"control_set": {"B 1": "0.02", "B 10": "0.0 4.0 0.0 0.02 0.0 2.5"},
                        "ctrl_set": {"P 51": ["%.4f 0.01" % x for x in [0.6, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.25, 3.4]]},
                        "generate": {"script": "thacker_geo.py", "args": [0.02]},
                        "check_set": {"tol": 0.08}},
            "release": {}},
    }
    write_case("sflow_thacker_parabola", control, ctrl, case)


def fnpf_jonswap():
    control = """
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 800 1 10
B 10 0.0 40.0 0.0 0.05 0.0 1.0
B 103 5
B 113 2.5
B 116 1.0

M 10 1
M 20 2
"""
    ctrl = """
A 10 3
A 310 3
A 311 4
A 343 0
B 84 2
B 85 2
B 86 512
B 88 3.3
B 90 1
B 92 31
B 93 0.1 1.5
B 96 4.0 10.0
B 98 2
B 99 1
B 139 1234
F 60 1.0
N 41 1100.0
N 47 1.0
M 10 1
P 12 1000
P 51 10.0 0.025
P 51 15.0 0.025
P 51 20.0 0.025
W 22 -9.81
"""
    case = {
        "description": "Irregular long-crested waves from a JONSWAP spectrum (Hs = 0.1 m, Tp = 1.5 s, gamma = 3.3, "
                       "d = 1 m, 512 components with the equal-energy discretisation, fixed seed), FNPF: spectral "
                       "significant wave height Hm0, peak period and spectral shape at x = 10, 15, 20 m against the "
                       "target spectrum (band-averaged periodogram of 1024 s / 512 s of record).",
        "reference": ref("jonswap"),
        "tags": ["fnpf", "analytical", "2dv", "irregular"],
        "time": 1100.0,
        "check": {"type": "spectrum", "Hs": 0.1, "Tp": 1.5, "gamma": 3.3, "t_start": 60.0, "t_end": 1100.0,
                  "dt": 0.05, "nband": 16, "tol_Hm0": 0.05, "tol_Tp": 0.08, "tol_shape": 0.2},
        "levels": {
            "nightly": {"time": 600.0, "control_set": {"B 2": "400 1 8"},
                        "check_set": {"t_start": 60.0, "t_end": 600.0, "nband": 8, "tol_Hm0": 0.08, "tol_Tp": 0.1,
                                      "tol_shape": 0.3}},
            "release": {}},
    }
    write_case("fnpf_jonswap_spectrum", control, ctrl, case)


def sloshing(model):
    L, h, a = 1.0, 0.5, 0.005
    def strips(dx):
        n = int(round(L / dx))
        out = []
        for i in range(n):
            xc = (i + 0.5) * dx
            z = h + a * math.cos(math.pi * xc / L)
            if model == "cfd":
                out.append("%.5f %.5f -1.0 1.0 0.0 %.6f" % (i * dx, (i + 1) * dx, z))
            else:
                out.append("%.5f %.5f -1.0 1.0 %.6f" % (i * dx, (i + 1) * dx, z))
        return out
    key = "F 70" if model == "cfd" else "F 72"
    if model == "cfd":
        control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 21

B 1 0.01
B 10 0.0 1.0 0.0 0.01 0.0 0.8

M 10 1
M 20 2
"""
        ctrl = """
D 10 4
D 20 2
D 30 1
F 30 3
F 40 3
N 40 3
N 41 12.0
N 47 0.2
M 10 1
P 12 100
T 10 0
W 22 -9.81
P 51 0.1667 0.005
"""
        night = {"control_set": {"B 1": "0.02", "B 10": "0.0 1.0 0.0 0.02 0.0 0.8"}}
    else:
        control = """
C 11 21
C 12 3
C 13 3
C 14 21
C 15 21
C 16 3

B 2 100 1 10
B 10 0.0 1.0 0.0 0.01 0.0 1.0

M 10 1
M 20 2
"""
        ctrl = """
A 10 5
F 60 0.5
N 41 12.0
N 47 0.5
M 10 1
P 12 100
W 22 -9.81
P 51 0.05 0.005
"""
        night = {"control_set": {"B 2": "50 1 5"}, "ctrl_set": {key: strips(0.02)}}
    if model == "cfd":
        # tilted plane z = h + a (1 - 2x/L): level set phi = z_s(x) - z from F 60 / F 62 / F 63; excites the odd modes, mode 1 dominates and
        # mode 3 has a node at the gauge x = L/6
        ctrl = ctrl + "F 60 %.4f\nF 62 %.4f\nF 63 0.0\n" % (h + a, h - a)
        night.pop("ctrl_set", None)
        night["ctrl_set"] = {"P 51": "%.4f 0.01" % (L / 6)}
    else:
        ctrl = ctrl + "\n".join("%s %s" % (key, s) for s in strips(0.01)) + "\n"
    case = {
        "description": "Free oscillation of the first sloshing mode in a closed rectangular tank (L = 1 m, h = 0.5 m, "
                       "kh = pi/2), %s: initial surface %s; period at a gauge against linear theory (T = 1.1816 s)." % (
                           model.upper(), "tilted plane with 5 mm at the walls (F 60 / F 62 / F 63), gauge at x = L/6 where the "
                           "third mode has a node" if model == "cfd" else "first-mode cosine with 5 mm amplitude, "
                           "set column by column (F 72), gauge at x = 0.05 m"),
        "reference": ref("sloshing"),
        "tags": [model, "analytical", "2dv", "quick"],
        "time": 12.0,
        "check": {"type": "oscillation_period", "gauge": {"index": 1}, "L": L, "h": h, "datum": 0.0 if model == "nhflow" else h,
                  "tol": 0.01},
        "levels": {"nightly": dict(night, check_set={"tol": 0.02}), "release": {}},
    }
    write_case("%s_sloshing_linear" % model, control, ctrl, case)


def main():
    if os.path.isdir(CASES):
        shutil.rmtree(CASES)
    sflow_dambreak("sflow_dambreak_ritter", 0.0, "ritter")
    sflow_dambreak("sflow_dambreak_stoker", 0.05, "stoker")
    sflow_solitary()
    sflow_beji()
    nhflow_beji()
    nhflow_stokes5()
    nhflow_synolakis()
    berkhoff("nhflow")
    berkhoff("fnpf")
    fnpf_beji()
    fnpf_stokes2()
    fnpf_ting_kirby()
    cfd_martin_moyce()
    cfd_kleefsman2()
    cfd_cylinder()
    cfd_beji()
    cfd_ting_kirby()
    cfd_chen()
    cfd_irschik()
    cfd_sphere()
    cfd_pier_scour()
    for m, c, hd, dxr, dxn in [("sflow", "A", 0.045, 0.05, 0.1), ("sflow", "C", 0.181, 0.05, 0.1),
                               ("nhflow", "C", 0.181, 0.05, 0.1)]:
        conical_geo_case(m, c, hd, dxr, dxn)
    for m in ("nhflow", "fnpf", "sflow"):
        mase_kirby(m)
    cfd_porous_dambreak()
    sflow_thacker()
    fnpf_jonswap()
    sloshing("cfd")
    sloshing("nhflow")
    print("\n".join(sorted(os.listdir(CASES))))


if __name__ == "__main__":
    main()
